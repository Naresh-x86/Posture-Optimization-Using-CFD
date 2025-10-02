import vtk
import numpy as np

def read_vtu_or_vtp(filename):
    # Choose reader based on file extension
    if filename.endswith('.vtu'):
        reader = vtk.vtkXMLUnstructuredGridReader()
    elif filename.endswith('.vtp'):
        reader = vtk.vtkXMLPolyDataReader()
    else:
        raise ValueError("Unsupported file type")

    reader.SetFileName(filename)
    reader.Update()
    return reader.GetOutput()

def average_velocity_magnitude(vtkdata):
    # Check if 'velocity' array exists
    vel_array = vtkdata.GetPointData().GetVectors('velocity')
    if vel_array is None:
        # Try lowercase 'Velocity'
        vel_array = vtkdata.GetPointData().GetVectors('Velocity')
    if vel_array is None:
        print("No velocity field found. Using zeros.")
        return 0.0

    n_points = vtkdata.GetNumberOfPoints()
    magnitudes = []

    for i in range(n_points):
        v = vel_array.GetTuple3(i)
        magnitudes.append(np.linalg.norm(v))

    return np.mean(magnitudes)

if __name__ == "__main__":
    filenames = [
        r"D:\SimVascularProjects\Results\average_result.vtp",
        r"D:\SimVascularProjects\Results\sim1_00010.vtp",
        r"D:\SimVascularProjects\Results\sim1_00010.vtu"
    ]

    for f in filenames:
        data = read_vtu_or_vtp(f)
        avg_vel = average_velocity_magnitude(data)
        print(f"{f}: Average velocity magnitude = {avg_vel:.6f}")
