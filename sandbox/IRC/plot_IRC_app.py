import os
import streamlit as st
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import py3Dmol

def atomic_number_to_symbol(atomic_number):
    """Convert atomic number to element symbol."""
    element = Chem.rdchem.GetPeriodicTable().GetElementSymbol(atomic_number)
    return element

# Function to parse the Gaussian log file
def parse_gaussian_log(file_content):
    lines = file_content.splitlines()

    geometries = []
    energies = []
    temp_geom = []
    capture_geom = False
    skip_lines_after_marker = 4  # Assuming there are 4 lines to skip after the geometry section marker

    for line in lines:
        if "Input orientation:" in line or "Standard orientation:" in line:
            capture_geom = True
            skipped_lines = 0
            if temp_geom:  # If there's an existing geometry block, store it before starting a new one
                geometries.append("\n".join(temp_geom))
                temp_geom = []
            continue

        if capture_geom:
            if skipped_lines < skip_lines_after_marker:
                # Skip header lines immediately following the marker
                skipped_lines += 1
                continue
            elif "---------------------------------------------------------------------" in line:
                # End of geometry section
                if temp_geom:  # Check to ensure we don't capture empty geometries
                    geometries.append("\n".join(temp_geom))
                    temp_geom = []
                capture_geom = False
            else:
                # This line is part of the geometry data; append it to the current geometry list
                temp_geom.append(line.strip())

        if "SCF Done:" in line:
            parts = line.split()
            energy = parts[4]  # Extract the energy value
            energies.append(float(energy))

    # Check for any remaining geometry data that hasn't been appended
    if temp_geom:
        geometries.append("\n".join(temp_geom))

    return geometries, energies

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
    """Plot the energy profile of the geometries."""
    fig, ax = plt.subplots()
    ax.plot(range(len(energies)), energies, marker='o', linestyle='-')
    ax.set_xlabel('Step')
    ax.set_ylabel('Energy (a.u.)')
    ax.set_title('Energy Profile along IRC')
    st.pyplot(fig)

def visualize_molecule(geometry, style='stick'):
    """Visualize the molecule using py3Dmol."""
    view = py3Dmol.view(width=400, height=300)
    view.addModel(geometry, 'xyz')
    view.setStyle({style: {}})
    view.zoomTo()
    return view

def main():
    st.title('Gaussian IRC Visualization App')

    uploaded_file = st.file_uploader("Upload a Gaussian IRC .log file", type="log")

    if uploaded_file is not None:
        file_content = uploaded_file.read().decode('utf-8')
        geometries, energies = parse_gaussian_log(file_content)

        if not geometries or not energies:
            st.error("Failed to parse geometries or energies. Please upload a valid Gaussian IRC log file.")
            return

        cleaned_geometries = clean_geometries_for_py3dmol(geometries)

        st.header("Energy Profile")
        plot_energies(energies)

        st.header("Geometry Viewer")

        step = st.slider('Select Geometry Step', 0, len(cleaned_geometries) - 1, 0)

        # Display selected geometry
        selected_geometry = cleaned_geometries[step]
        st.write(f'### Geometry at step {step}')
        st.text(selected_geometry)

        # 3D Visualization
        xyzview = visualize_molecule(selected_geometry)
        st.components.v1.html(xyzview._make_html(), width=400, height=300, scrolling=False)

if __name__ == "__main__":
    main()

