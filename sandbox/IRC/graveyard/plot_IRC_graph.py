import streamlit as st
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from io import BytesIO
import base64
import py3Dmol

def atomic_number_to_symbol(atomic_number):
    element = Chem.rdchem.GetPeriodicTable().GetElementSymbol(atomic_number)
    return element

# Function to parse the Gaussian log file
def parse_gaussian_log(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

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

    print("Geometries:",geometries)
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
    print("cleaned geometries:", cleaned_geometries) 
    return cleaned_geometries

def show_structure(geometry):
    # Initialize a Py3Dmol view
    view = py3Dmol.view(width=400, height=300)
    
    # Add the molecule in XYZ format
    view.addModel(geometry, 'xyz')
    
    # Style the molecule. You can customize this according to your needs.
    view.setStyle({'stick': {}})
    
    # Zoom to fit the molecule in the view
    view.zoomTo()
    
    # Return the HTML needed to render the view. Streamlit's `components.html` can be used to display this.
    return view.show()

def convert_atomic_numbers_to_symbols(geometry_lines):
    from rdkit.Chem.PeriodicTable import GetElementSymbol
    converted_lines = []
    for line in geometry_lines:
        parts = line.split()
        if len(parts) == 5:  # Atomic number and coordinates
            atom_symbol = GetElementSymbol(int(parts[0]))
            converted_line = f"{atom_symbol} {parts[2]} {parts[3]} {parts[4]}"
            converted_lines.append(converted_line)
    return converted_lines

# Function to plot energies
def plot_energies(energies):
    fig, ax = plt.subplots()
    ax.plot(energies, marker='o', linestyle='-')
    ax.set_xlabel('Step')
    ax.set_ylabel('Energy (a.u.)')
    return fig

# Streamlit app main function
def main():
    st.title('Gaussian IRC Visualization')

    file_path = "H-N-IRC.log"
    geometries, energies = parse_gaussian_log(file_path)
    geometries = clean_geometries_for_py3dmol(geometries)
    
    # Plotting energies
    if energies:
        st.write("Energy Profile")
        fig = plot_energies(energies)
        st.pyplot(fig)
    
    # Visualizing molecular structures
    if geometries:
        point_to_view = st.slider("Select a point to view", 0, len(geometries)-1, 0)
        st.write("Molecular Structure at Selected Point")
        html_str = show_structure(geometries[point_to_view])
        st.markdown(html_str, unsafe_allow_html=True)
    else:
        st.write("No geometries to display. Please check the log file and parsing logic.")

if __name__ == "__main__":
    main()

