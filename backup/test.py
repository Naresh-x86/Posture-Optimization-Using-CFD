import subprocess
import os
from pathlib import Path
import pyvista as pv  # pip install pyvista

case_dir = Path(r"D:\SimVascularProjects\0007_H_AO_H\Simulations\Project_Simulation\mesh-complete")

def run_simulation(case_dir):
    # Step 1: prepare data
    subprocess.run([r"C:\Program Files\SimVascular\svSolver\2022-08-19\svpre-bin.exe"], cwd=case_dir, check=True)
    # Step 2: run solver
    subprocess.run([r"C:\Program Files\SimVascular\svSolver\2022-08-19\svsolver-nompi-bin.exe"], cwd=case_dir, check=True)
    # Step 3: run post
    subprocess.run([r"C:\Program Files\SimVascular\svSolver\2022-08-19\svpost-bin.exe"], cwd=case_dir, check=True)

def extract_outflow_results(case_dir, outlet_names):
    """
    Read .vtu results from all_results/ and compute average velocity magnitude
    across each outlet surface.
    """
    results_dir = case_dir / "all_results"
    last_step = sorted(results_dir.glob("*.vtu"))[-1]  # pick last timestep
    
    mesh = pv.read(last_step)
    print(f"[OK] Loaded {last_step}")

    outlet_data = {}
    for outlet in outlet_names:
        # Each outlet surface is usually a group named in SimVascular
        if outlet not in mesh.cell_data:
            print(f"[WARN] No surface data found for {outlet}")
            continue

        # Average velocity magnitude for that surface
        vel = mesh.point_data["velocity"]  # velocity field
        mag = (vel**2).sum(axis=1)**0.5
        outlet_data[outlet] = mag.mean()
    
    return outlet_data

if __name__ == "__main__":
    run_simulation(case_dir)
    outlets = ["btrunk", "rt_carotid", "carotid", "subclavian", "outflow"]
    data = extract_outflow_results(case_dir, outlets)
    print("\n=== Outflow results ===")
    for k, v in data.items():
        print(f"{k}: {v:.3f} cm/s (avg velocity magnitude)")
