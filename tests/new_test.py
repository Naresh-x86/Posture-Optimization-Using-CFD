import os
import subprocess

# Path to your case folder (where solver.inp lives)
case_dir = r"D:\SimVascularProjects\0007_H_AO_H\Simulations\Project_Simulation"

def run_simulation(case_dir):
    # Preprocess
    subprocess.run(["svpre-bin.exe"], cwd=case_dir, check=True)

    # Solve (no MPI version)
    subprocess.run(["svsolver-nompi-bin.exe"], cwd=case_dir, check=True)

    # Postprocess
    subprocess.run(["svpost-bin.exe"], cwd=case_dir, check=True)

run_simulation(case_dir)
