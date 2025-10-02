import pickle
import os
import subprocess
import re
import shutil
from collections import deque
from dataclasses import dataclass
from IPython.display import display
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time

AORTA_PROJECT = r'D:\SimVascularProjects\0011_H_AO_H'
ABDOMINAL_PROJECT = r'D:\SimVascularProjects\0029_H_ABAO_H'
CORONARY_PROJECT = r'D:\SimVascularProjects\0073_H_CORO_H'
RESULTS_FOLDER = r'D:\SimVascularProjects\Results'
DEFAULT_RESULTS_DIR = os.path.join(os.getcwd(), 'default_results')

# SIMVASCULAR PATHS - Update these for your installation
MPI_PATH = r'C:\Program Files\Microsoft MPI\Bin\mpiexec.exe'
SVSOLVER_PATH = r'C:\Program Files\SimVascular\svSolver\2022-08-19'
SVSOLVER_EXE = 'svsolver-msmpi-bin.exe'
SVPOST_EXE = 'svpost-bin.exe'

# SIMULATION SETTINGS
N_PROCESSES = 8
SIM_SUBFOLDER = r'Simulations\Project'
SOLVER_FILE = 'solver.inp'
INFLOW_FILE = 'inflow.flow'
ITERATION_DURATION_MIN = 30
MAX_STATIONARY_WINDOW = 3

print("Configuration loaded successfully")
print(f"Aorta Project: {AORTA_PROJECT}")
print(f"SimVascular Path: {SVSOLVER_PATH}")
print(f"MPI Processes: {N_PROCESSES}")

# Vessel Configuration Classes
@dataclass
class PostureConfig:
    neck: float = 0
    torso: float = 0
    legs: float = 0

    def __post_init__(self):
        self.neck = max(-30, min(45, self.neck))
        self.torso = max(-15, min(60, self.torso))
        self.legs = max(-30, min(30, self.legs))

    def __str__(self):
        return f"Neck:{self.neck:+.0f}° Torso:{self.torso:+.0f}° Legs:{self.legs:+.0f}°"

    def as_tuple(self):
        return (self.neck, self.torso, self.legs)


def posture_distance(p1, p2):
    return abs(p1.neck - p2.neck) + abs(p1.torso - p2.torso) + abs(p1.legs - p2.legs)


class VesselConfig:
    # Aorta Configuration
    AORTA_SURFACES = {2: 'btrunk', 3: 'carotid', 5: 'outflow', 6: 'rt_carotid', 7: 'subclavian'}
    AORTA_RESISTANCE = {2: 12667, 3: 25333, 5: 2171, 6: 25333, 7: 25333}

    # Abdominal Aorta Configuration
    ABDOMINAL_SURFACES = {
        2: 'SMA', 3: 'hepatic', 5: 'left_external_iliac', 6: 'left_internal_iliac',
        7: 'left_renal', 8: 'right_external_iliac', 9: 'right_internal_iliac',
        10: 'right_renal', 11: 'splenic'
    }
    ABDOMINAL_RESISTANCE = {
        2: 50666, 3: 152000, 5: 20266, 6: 76000, 7: 50666,
        8: 20266, 9: 76000, 10: 50666, 11: 152000
    }

    # Coronary Configuration
    CORONARY_SURFACES = {
        2: 'aorta_outlet', 4: 'lca_br1', 5: 'lca_br2', 6: 'lca_br3',
        7: 'lca_br4', 8: 'lca_br5', 9: 'lca_br6', 10: 'rca_br1',
        11: 'rca_br2', 12: 'rca_br3'
    }
    CORONARY_RESISTANCE = {
        2: 2171, 4: 142857, 5: 200000, 6: 200000, 7: 333333,
        8: 500000, 9: 500000, 10: 150000, 11: 250000, 12: 375000
    }


def _fallback_distribution(surface_map, resistance_map):
    weights = {}
    for sid, name in surface_map.items():
        base_res = resistance_map.get(sid, 1.0)
        weights[name] = 1.0 / base_res if base_res else 1.0
    total = sum(weights.values())
    if not total:
        return {name: 1.0 / len(surface_map) for name in surface_map.values()}
    return {name: value / total for name, value in weights.items()}


def load_baseline_flow_profiles():
    mapping = {
        'aorta': ('Aorta', VesselConfig.AORTA_SURFACES, VesselConfig.AORTA_RESISTANCE),
        'abdominal': ('Abdominal_aorta', VesselConfig.ABDOMINAL_SURFACES, VesselConfig.ABDOMINAL_RESISTANCE),
        'coronary': ('Coronary', VesselConfig.CORONARY_SURFACES, VesselConfig.CORONARY_RESISTANCE)
    }
    profiles = {}
    for key, (folder, surface_map, resistance_map) in mapping.items():
        flows_path = os.path.join(DEFAULT_RESULTS_DIR, folder, 'all_results-flows.txt')
        distribution = {}
        baseline_inflow = 83.3
        try:
            df = pd.read_csv(flows_path, sep='\t')
            row = df.iloc[0].drop(labels='step')
            inflow = abs(row.get('inflow', baseline_inflow))
            if inflow > 0:
                baseline_inflow = inflow
            outlets = {col: abs(val) for col, val in row.items() if col != 'inflow'}
            total_out = sum(outlets.values())
            if total_out > 0:
                distribution = {col: val / total_out for col, val in outlets.items()}
        except FileNotFoundError:
            distribution = {}
        if not distribution:
            distribution = _fallback_distribution(surface_map, resistance_map)
        profiles[key] = {
            'baseline_inflow': baseline_inflow,
            'distribution': distribution
        }
    return profiles


BASELINE_PROFILES = load_baseline_flow_profiles()
FLOW_RATIO_BY_VESSEL = {'aorta': 0.94, 'abdominal': 0.90, 'coronary': 0.86}


@dataclass
class FlowParameters:
    vessel_key: str
    inflow: float
    resistance_multiplier: float
    expected_outflow: float
    outflow_ratio: float
    outlet_distribution: dict


def compute_expected_outflow(vessel_key, inflow, resistance_multiplier):
    base_ratio = FLOW_RATIO_BY_VESSEL.get(vessel_key, 0.9)
    posture_factor = 1.05 - 0.12 * (resistance_multiplier - 1.0)
    ratio = base_ratio * posture_factor
    ratio = max(0.55, min(0.97, ratio))
    return inflow * ratio, ratio


def scale_distribution(distribution, total_outflow):
    return {name: total_outflow * fraction for name, fraction in distribution.items()}


def calculate_flow_parameters(posture, vessel_type):
    vessel_key = vessel_type.lower()
    base_flow = BASELINE_PROFILES.get(vessel_key, {}).get('baseline_inflow', 83.3)

    co_factor = 1.0 + posture.torso * 0.002 + posture.legs * 0.001 - abs(posture.neck) * 0.0005
    inflow = base_flow * max(0.7, min(1.3, co_factor))

    if vessel_key == 'aorta':
        resistance_mult = 1.0 + posture.neck * 0.001 + posture.torso * 0.0005
    elif vessel_key == 'abdominal':
        resistance_mult = 1.0 + posture.torso * 0.001 + posture.legs * 0.002
    elif vessel_key == 'coronary':
        resistance_mult = 1.0 + posture.torso * 0.002
    else:
        resistance_mult = 1.0

    resistance_mult = max(0.7, min(1.3, resistance_mult))

    expected_outflow, ratio = compute_expected_outflow(vessel_key, inflow, resistance_mult)
    outlet_distribution = BASELINE_PROFILES.get(vessel_key, {}).get('distribution', {})
    if not outlet_distribution:
        if vessel_key == 'aorta':
            outlet_distribution = _fallback_distribution(VesselConfig.AORTA_SURFACES, VesselConfig.AORTA_RESISTANCE)
        elif vessel_key == 'abdominal':
            outlet_distribution = _fallback_distribution(VesselConfig.ABDOMINAL_SURFACES, VesselConfig.ABDOMINAL_RESISTANCE)
        else:
            outlet_distribution = _fallback_distribution(VesselConfig.CORONARY_SURFACES, VesselConfig.CORONARY_RESISTANCE)

    return FlowParameters(vessel_key, inflow, resistance_mult, expected_outflow, ratio, outlet_distribution)

# Load the variable from the pickle file
with open('run_summary.pkl', 'rb') as f:
    run_summary = pickle.load(f)

# Print the full variable
print(run_summary)
