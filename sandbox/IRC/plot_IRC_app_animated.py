import os
import numpy as np
import streamlit as st
import matplotlib.pyplot as plt
from rdkit import Chem
import py3Dmol
import time

def atomic_number_to_symbol(atomic_number):
    """Convert atomic number to element symbol."""
    element = Chem.rdchem.GetPeriodicTable().GetElementSymbol(atomic_number)
    return element

# Function to parse the Gaussian log file
def parse_gaussian_log(file_content):
    lines = file_content.splitlines()

    ts_geometry = None
    forward_geometries = []
    reverse_geometries = []
    energies = []

    temp_geom = []
    capture_geom = False
    skip_lines_after_marker = 4  # Assuming there are 4 lines to skip after the geometry section marker
    reading_reverse = False

    for line in lines:
        if "Point Number  1 in FORWARD path direction." in line:
            capture_geom = True
            reading_reverse = False
            skipped_lines = 0
            continue
        elif "Point Number  1 in REVERSE path direction." in line:
            capture_geom = True
            reading_reverse = True
            skipped_lines = 0
            continue
        elif "Input orientation:" in line or "Standard orientation:" in line:
            capture_geom = True
            skipped_lines = 0
            if temp_geom:
                if reading_reverse:
                    reverse_geometries.append("\n".join(temp_geom))
                else:
                    forward_geometries.append("\n".join(temp_geom))
                temp_geom = []
            continue

        if capture_geom:
            if skipped_lines < skip_lines_after_marker:
                # Skip header lines immediately following the marker
                skipped_lines += 1
                continue
            elif "---------------------------------------------------------------------" in line:
                # End of geometry section
                if temp_geom:
                    if reading_reverse:
                        reverse_geometries.append("\n".join(temp_geom))
                    else:
                        forward_geometries.append("\n".join(temp_geom))
                    temp_geom = []
                capture_geom = False
            else:
                # This line is part of the geometry data; append it to the current geometry list
                temp_geom.append(line.strip())

        if "SCF Done:" in line:
            parts = line.split()
            energy = parts[4]  # Extract the energy value
            if reading_reverse:
                energies.insert(0, float(energy))  # Insert energy at the beginning for reverse
            else:
                energies.append(float(energy))

    # Check for any remaining geometry data that hasn't been appended
    if temp_geom:
        if reading_reverse:
            reverse_geometries.append("\n".join(temp_geom))
        else:
            forward_geometries.append("\n".join(temp_geom))

    # Reverse the list of reverse geometries to have them in the correct order
    reverse_geometries.reverse()

    # Combine all geometries: reverse (in reversed order) -> TS -> forward
    combined_geometries = reverse_geometries + forward_geometries

    return combined_geometries, energies

def compute_centroid(geometry):
    """Compute the centroid of a geometry."""
    lines = geometry.splitlines()
    coordinates = np.array([[float(parts.split()[1]), float(parts.split()[2]), float(parts.split()[3])] for parts in lines])
    centroid = np.mean(coordinates, axis=0)
    return centroid

def align_geometries(geometries):
    """Align all geometries to a common reference frame."""
    aligned_geometries = []

    for geometry in geometries:
        lines = geometry.splitlines()
        # Compute the centroid
        centroid = compute_centroid(geometry)
        coordinates = []

        # Translate all atoms to center at the origin
        for line in lines:
            parts = line.split()
            element_symbol = parts[0]
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
            x -= centroid[0]
            y -= centroid[1]
            z -= centroid[2]
            coordinates.append([element_symbol, x, y, z])

        # Convert coordinates to numpy array for matrix operations
        coord_array = np.array([[x, y, z] for _, x, y, z in coordinates])

        # Compute the moment of inertia tensor
        inertia_tensor = np.zeros((3, 3))
        for _, x, y, z in coordinates:
            inertia_tensor += np.array([[y**2 + z**2, -x * y, -x * z],
                                        [-y * x, x**2 + z**2, -y * z],
                                        [-z * x, -z * y, x**2 + y**2]])

        # Eigen decomposition for principal axes
        eigenvalues, eigenvectors = np.linalg.eigh(inertia_tensor)

        # Align the molecule to the principal axes
        aligned_coords = np.dot(coord_array, eigenvectors)
        aligned_geometry = "\n".join([f"{coordinates[i][0]} {aligned_coords[i][0]:.6f} {aligned_coords[i][1]:.6f} {aligned_coords[i][2]:.6f}" for i in range(len(coordinates))])
        aligned_geometries.append(aligned_geometry)

    return aligned_geometries

def clean_geometries_for_py3dmol(geometries):
    cleaned_geometries = []

    for geometry in geometries:
        lines = geometry.split('\n')
        cleaned_geometry = []
        for line in lines:
            parts = line.split()
            if len(parts) < 6:
                continue  # Skip lines that don't have enough data
            atomic_num = int(parts[1])
            x, y, z = parts[3], parts[4], parts[5]
            element_symbol = atomic_number_to_symbol(atomic_num)
            cleaned_geometry.append(f"{element_symbol} {x} {y} {z}")
        cleaned_geometries.append("\n".join(cleaned_geometry))
    return cleaned_geometries

def plot_energies(energies):
    """Plot the energy profile of the geometries with optional normalization and unit conversion."""
    
    # Normalize energies if the user selects the option
    normalize = st.checkbox('Normalize energies relative to the lowest value')
    normalize = True
    if normalize:
        # Convert energies to be relative to the lowest value
        min_energy = min(energies)
        energies = [energy - min_energy for energy in energies]
        
        # Unit conversion options
        unit = st.radio('Select energy units:', ('Hartree', 'eV', 'Kcal/mol'))
        
        if unit == 'eV':
            energies = [energy * 27.21 for energy in energies]  # Convert Hartree to eV
        elif unit == 'Kcal/mol':
            energies = [energy * 627.51 for energy in energies]  # Convert Hartree to Kcal/mol

    # Plotting the energy profile
    fig, ax = plt.subplots()
    ax.plot(range(len(energies)), energies, marker='o', linestyle='-')
    ax.set_xlabel('Step')
    ax.set_ylabel(f'Energy ({unit})')
    ax.set_title('Energy Profile along IRC')
    
    st.pyplot(fig)

def visualize_molecule(geometry, style='stick'):
    """Visualize the molecule using py3Dmol."""
    # Calculate number of atoms
    lines = geometry.splitlines()
    num_atoms = len(lines)

    # Reformat geometry to standard XYZ format
    xyz_format = f"{num_atoms}\n\n" + geometry

    # Style options for visualization
    style_options = {
        'Ball and Stick': {'stick': {}, 'sphere': {'radius': 0.5}},
        'Stick': {'stick': {}},
        'Spacefill': {'sphere': {}}
    }
    selected_style = st.radio('Select visualization style', list(style_options.keys()))

    # Visualize the molecule using py3Dmol
    scale = 1
    width = int(640.0 * scale)
    height = int(480.0 * scale)

    xyzview = py3Dmol.view(width=width, height=height)
    xyzview.addModel(xyz_format, 'xyz')
    xyzview.setStyle(style_options[selected_style])  # Use the selected style
    xyzview.zoomTo()

    # Display the visualization in Streamlit
    st.components.v1.html(xyzview._make_html(), width=width, height=height, scrolling=False)

def main():
    st.title('Gaussian IRC Visualization App')
    st.markdown("### Beta version")

    uploaded_file = st.file_uploader("Upload a Gaussian IRC .log file", type="log")

    if uploaded_file is not None:
        file_content = uploaded_file.read().decode('utf-8')
        geometries, energies = parse_gaussian_log(file_content)

        if not geometries or not energies:
            st.error("Failed to parse geometries or energies. Please upload a valid Gaussian IRC log file.")
            return

        cleaned_geometries = clean_geometries_for_py3dmol(geometries)

        # Align all geometries to the same reference frame
        aligned_geometries = align_geometries(cleaned_geometries)

        st.header("Energy Profile")
        plot_energies(energies)

        st.header("Geometry Viewer")

        if 'step' not in st.session_state:
            st.session_state.step = 0
            st.session_state.play = False

        # Play/Pause button
        if st.button('Play Animation'):
            st.session_state.play = True

        # Stop button
        if st.button('Stop Animation'):
            st.session_state.play = False

        # Animation loop
        if st.session_state.play:
            st.session_state.step += 1
            if st.session_state.step >= len(aligned_geometries):
                st.session_state.step = 0  # Loop back to start
            time.sleep(0.2)
            st.experimental_rerun()

        # Slider for manual selection
        step = st.slider('Select Geometry Step', 0, len(aligned_geometries) - 1, st.session_state.step)
        st.session_state.step = step

        # 3D Visualization of selected geometry
        selected_geometry = aligned_geometries[step]
        st.write(f'### Geometry at step {step}')
        visualize_molecule(selected_geometry)

if __name__ == "__main__":
    main()


# Animation button still doesnt work
# 
