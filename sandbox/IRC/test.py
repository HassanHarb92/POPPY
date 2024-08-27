import streamlit as st
import matplotlib.pyplot as plt
import py3Dmol
from rdkit import Chem

def atomic_number_to_symbol(atomic_number):
    element = Chem.rdchem.GetPeriodicTable().GetElementSymbol(atomic_number)
    return element

def parse_gaussian_log(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    geometries = []
    energies = []
    temp_geom = []
    capture_geom = False
    skip_lines_after_marker = 4

    for line in lines:
        if "Input orientation:" in line or "Standard orientation:" in line:
            capture_geom = True
            skipped_lines = 0
            if temp_geom:
                geometries.append("\n".join(temp_geom))
                temp_geom = []
            continue

        if capture_geom:
            if skipped_lines < skip_lines_after_marker:
                skipped_lines += 1
                continue
            elif "---------------------------------------------------------------------" in line:
                if temp_geom:
                    geometries.append("\n".join(temp_geom))
                    temp_geom = []
                capture_geom = False
            else:
                temp_geom.append(line.strip())

        if "SCF Done:" in line:
            parts = line.split()
            energy = parts[4]
            energies.append(float(energy))

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
                continue
            atomic_num = int(parts[1])
            x, y, z = parts[3], parts[4], parts[5]
            element_symbol = atomic_number_to_symbol(atomic_num)
            cleaned_geometry.append(f"{element_symbol} {x} {y} {z}")
        cleaned_geometries.append("\n".join(cleaned_geometry))

    return cleaned_geometries

def show_structure(geometry):
    # This remains unchanged, setting up the molecular visualization
    view = py3Dmol.view(width=400, height=300)
    view.addModel(geometry, 'xyz')
    view.setStyle({'stick': {}})
    view.zoomTo()
    # Instead of returning view.show(), just prepare the view for rendering
    return view

def main():
    st.title('Gaussian IRC and XYZ Visualization')

    # UI to select mode: Parse log or view XYZ
    mode = st.sidebar.selectbox("Mode", ["View XYZ Files", "Parse Gaussian Log"])

    if mode == "View XYZ Files":
        # Assuming your xyz files are in the "xyz_files" directory
        xyz_files_directory = 'xyz_files'
        xyz_files = list_xyz_files(xyz_files_directory)
        if not xyz_files:
            st.write("No .xyz files found in the xyz_files directory.")
            return
        selected_file = st.selectbox('Select an XYZ file', xyz_files)
        selected_style = st.radio('Select visualization style', ['Stick', 'Ball and Stick', 'Spacefill'])
        if st.button('Visualize'):
            full_path_to_file = os.path.join(xyz_files_directory, selected_file)
            xyz_content = read_xyz_file(full_path_to_file)
            show_molecule(xyz_content, selected_style)

    elif mode == "Parse Gaussian Log":
        file_path = st.text_input("Enter the path to the Gaussian log file", "H-N-IRC.log")
        if st.button('Parse and Visualize'):
            geometries, energies = parse_gaussian_log(file_path)
            cleaned_geometries = clean_geometries_for_py3dmol(geometries)
            if cleaned_geometries:
                create_xyz_files(cleaned_geometries)  # Optional: Save to XYZ files
                point_to_view = st.slider("Select a point to view", 0, len(cleaned_geometries)-1, 0)
                html_str = show_structure(cleaned_geometries[point_to_view])
                st.markdown(html_str, unsafe_allow_html=True)
                if energies:
                    st.write("Energy Profile")
                    plot_energies(energies)

if __name__ == "__main__":
    main()

## Need to debug manually
## check separate scripts
## make sure to first of all order the molecules in the way we need them (TS) then (reverse) then (forward)
## have test cases to pick one structure to visualize with the plot then change it to incorporate multiple structures
