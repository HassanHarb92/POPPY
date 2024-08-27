import os
import streamlit as st
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from io import BytesIO
import base64
import py3Dmol

import pandas as pd

# Assuming energies is a list of energies extracted from the log file
# and xyz_file_paths is a list of paths to the generated .xyz files

def create_energy_dataframe(xyz_file_paths, energies):
    """Create a DataFrame with XYZ filenames and their corresponding energies."""
    # Extract filenames from paths
    filenames = [os.path.basename(path) for path in xyz_file_paths]
    
    # Create the DataFrame
    df = pd.DataFrame({
        'Filename': filenames,
        'Energy': energies
    })
    
    return df


def atomic_number_to_symbol(atomic_number):
    element = Chem.rdchem.GetPeriodicTable().GetElementSymbol(atomic_number)
    return element

def parse_gaussian_log(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Initialize containers for parsed data
    geometries = {'initial': [], 'forward': [], 'reverse': []}
    energies = {'initial': None, 'forward': [], 'reverse': []}
    temp_geom = []
    capture_geom = False
    current_section = 'initial'  # Start with 'initial'

    for line in lines:
        if "Input orientation:" in line or "Standard orientation:" in line:
            capture_geom = True
            skipped_lines = 0  # Reset skipped lines counter for header lines
            if temp_geom:  # If there's an existing geometry block
                geometries[current_section].append("\n".join(temp_geom))
                temp_geom = []
            continue

        if capture_geom:
            if skipped_lines < 4:  # Assuming there are 4 header lines to skip
                skipped_lines += 1
                continue
            elif "---------------------------------------------------------------------" in line:
                capture_geom = False  # Stop capturing geometry
                if temp_geom:  # Ensure we don't append empty geometries
                    geometries[current_section].append("\n".join(temp_geom))
                    temp_geom = []
            else:
                temp_geom.append(line.strip())  # Append geometry line

        if "SCF Done:" in line:
            energy = float(line.split()[4])
            if current_section == 'initial' and energies['initial'] is None:
                energies['initial'] = energy
            else:
                energies[current_section].append(energy)

        # Identifying section markers for initial, forward, and reverse directions
        # Adjust these conditions based on your log file's specific markers
        if "Start of forward direction" in line:
            current_section = 'forward'
        elif "Start of reverse direction" in line:
            current_section = 'reverse'

    # Check if there are any remaining geometries not appended due to missing end line
    if temp_geom:
        geometries[current_section].append("\n".join(temp_geom))

    # The function now returns dictionaries for geometries and energies,
    # with keys for 'initial', 'forward', and 'reverse' containing the respective data.
    return geometries, energies

# Example usage
file_path = "H-N-IRC.log"  # Update this with your actual file path
geometries, energies = parse_gaussian_log(file_path)

# Accessing the parsed data
initial_structure = geometries['initial'][0] if geometries['initial'] else None
forward_structures = geometries['forward']
reverse_structures = geometries['reverse']
initial_energy = energies['initial']
forward_energies = energies['forward']
reverse_energies = energies['reverse']


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
    # Example code within the main function
    geometries, energies = parse_gaussian_log("H-N-IRC.log")  # Update path as needed
    cleaned_geometries = clean_geometries_for_py3dmol(geometries)
    xyz_file_paths = create_xyz_files(cleaned_geometries)
    energy_df = create_energy_dataframe(xyz_file_paths, energies)

    fig, ax = plt.subplots()
    ax.plot(energy_df['Energy'], marker='o', linestyle='-')
    ax.set_xlabel('Structure Index')
    ax.set_ylabel('Energy (a.u.)')
    ax.set_title('Energy Profile')
    st.pyplot(fig)

    # Replace the direct file selection with DataFrame-based selection
    file_index = st.slider('Select a structure to view', 0, len(energy_df) - 1, 0)
    selected_file = energy_df.loc[file_index, 'Filename']
    selected_energy = energy_df.loc[file_index, 'Energy']
    
    # Display selected energy
    st.write(f"Selected Energy: {selected_energy} a.u.")
    
    # Proceed with visualization as before
    full_path_to_file = os.path.join(xyz_files_directory, selected_file)
    xyz_content = read_xyz_file(full_path_to_file)
    # Visualization code remains the same...


    st.title('Gaussian IRC Visualization')

    xyz_files_directory = 'xyz_files'
    xyz_files = list_xyz_files(xyz_files_directory)

    if not xyz_files:
        st.write("No .xyz files found in the xyz_files directory.")
        return

    # Create a slider for selecting which XYZ file to view
    file_index = st.slider('Select a structure to view', 0, len(xyz_files) - 1, 0)

    # Get the selected file based on slider's position
    selected_file = xyz_files[file_index]

    style_options = {
        'Stick': {'stick': {}},
        'Ball and Stick': {'stick': {}, 'sphere': {'radius': 0.5}},
        'Spacefill': {'sphere': {}}
    }
    selected_style = st.radio('Select visualization style', list(style_options.keys()))

    # Full path to the selected file
    full_path_to_file = os.path.join(xyz_files_directory, selected_file)

    # Read the selected XYZ file
    xyz_content = read_xyz_file(full_path_to_file)
    width = 640
    height = 480

    # Visualize the molecule using py3Dmol
    xyzview = py3Dmol.view(query=f'{full_path_to_file}', width=width, height=height)
    xyzview.addModel(xyz_content, 'xyz')
    xyzview.setStyle(style_options[selected_style])  # Apply the selected style
    xyzview.zoomTo()

    # Display the visualization in Streamlit
    st.components.v1.html(xyzview._make_html(), width=width, height=height, scrolling=False)

if __name__ == "__main__":
    main()

