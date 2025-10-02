import pyvista as pv
import numpy as np
import csv
from pathlib import Path

# --- SETTINGS ---
input_files = [
    r"D:\SimVascularProjects\Results\average_result.vtp",
    r"D:\SimVascularProjects\Results\sim1_00010.vtp",
    r"D:\SimVascularProjects\Results\sim1_00010.vtu"
]
output_csv = r"D:\SimVascularProjects\Results\placeholder_flows.csv"

# Branch names (fake / placeholder)
branches = ["btrunk", "carotid", "inflow", "outflow", "rt_carotid", "subclavian"]

# --- HELPER FUNCTIONS ---
def avg_velocity_magnitude(mesh):
    """Compute average velocity magnitude if velocity field exists."""
    if "velocity" in mesh.array_names:
        vel = mesh["velocity"]
        return np.mean(np.linalg.norm(vel, axis=1))
    elif "Velocity" in mesh.array_names:
        vel = mesh["Velocity"]
        return np.mean(np.linalg.norm(vel, axis=1))
    else:
        return 0.0

# --- MAIN ---
with open(output_csv, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    # Header
    writer.writerow(["step"] + branches)
    
    # Placeholder step number
    step = 10
    
    for file in input_files:
        mesh = pv.read(file)
        avg_vel = avg_velocity_magnitude(mesh)
        
        # Create fake branch flows by distributing average velocity
        # (sum to avg_vel for simplicity)
        flows = np.linspace(avg_vel/2, avg_vel/20, len(branches))
        
        # Write row
        writer.writerow([step] + list(flows))

print(f"Placeholder flow table written to {output_csv}")
