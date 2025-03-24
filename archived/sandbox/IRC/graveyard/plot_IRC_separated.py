import os
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

def create_xyz_files(cleaned_geometries, directory="xyz_files"):
    # Create the directory if it doesn't exist
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    xyz_file_paths = []
    
    # Write each geometry to a separate .xyz file
    for i, geometry in enumerate(cleaned_geometries):
        file_path = os.path.join(directory, f"molecule_{i}.xyz")
        with open(file_path, 'w') as file:
            # Assuming the first line needs to be the atom count and the second line is a comment
            atom_count = geometry.count('\n') + 1  # +1 to count the last line as well
            file.write(f"{atom_count}\n")
            file.write(f"Molecule {i}\n")
            file.write(geometry)
        
        xyz_file_paths.append(file_path)
    
    return xyz_file_paths


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
import streamlit as st
import os
import py3Dmol

def list_xyz_files(directory):
    """List all .xyz files in the specified directory."""
    return [f for f in os.listdir(directory) if f.endswith('.xyz')]

def read_xyz_file(xyz_file_path):
    """Read and return the content of the specified XYZ file."""
    with open(xyz_file_path, 'r') as file:
        return file.read()

def main():
    st.title('Gaussian IRC Visualization')

    # Assuming your xyz files are in the "xyz_files" directory
    xyz_files_directory = 'xyz_files'
    xyz_files = list_xyz_files(xyz_files_directory)

    if not xyz_files:
        st.write("No .xyz files found in the xyz_files directory.")
        return

    # Dropdown to select an XYZ file
    selected_file = st.selectbox('Select an XYZ file', xyz_files)

    # Style options as radio buttons
    style_options = {
        'Stick': {'stick': {}},
        'Ball and Stick': {'stick': {}, 'sphere': {'radius': 0.5}},
        'Spacefill': {'sphere': {}}
    }
    selected_style = st.radio('Select visualization style', list(style_options.keys()))

    # Button to visualize the selected XYZ file
    if st.button('Visualize'):
        # Full path to the selected file
        full_path_to_file = os.path.join(xyz_files_directory, selected_file)

        # Read the selected XYZ file
        xyz_content = read_xyz_file(full_path_to_file)
        print(xyz_content)
        width = 640
        height = 480

        # Visualize the molecule using py3Dmol
        xyzview = py3Dmol.view(width=width, height=height)
        xyzview.addModel(xyz_content, 'xyz')
        xyzview.setStyle(style_options[selected_style])  # Use the selected style
        xyzview.zoomTo()

        # Display the visualization in Streamlit
        xyzview.show()
        st.components.v1.html(xyzview._make_html(), width=width, height=height, scrolling=False)

        # Display the visualization in Streamlit
       # html_str = xyzview.show()
       # st.components.v1.html(html_str, width=width, height=height, scrolling=False)

if __name__ == "__main__":
    main()

