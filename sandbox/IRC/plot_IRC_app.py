import os
import streamlit as st
import matplotlib.pyplot as plt
from rdkit import Chem
import py3Dmol

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

def old_visualize_molecule(geometry, style='stick'):
    print ("geometry debug")
    print (geometry)
    style_options = {
        'Stick': {'stick': {}},
        'Ball and Stick': {'stick': {}, 'sphere': {'radius': 0.5}},
        'Spacefill': {'sphere': {}}
    }
    selected_style = st.radio('Select visualization style', list(style_options.keys()))

    """Visualize the molecule using py3Dmol."""
    scale = 1
    width = int(640.0 * scale)
    height = int(480.0 * scale)

    xyzview = py3Dmol.view(width=width, height=height)
    xyzview.addModel(geometry, 'xyz')
#    xyzview.setStyle({style: {}})
    xyzview.setStyle(style_options[selected_style])  # Use the selected style
    xyzview.zoomTo()

    xyzview.show()
    # Display the visualization in Streamlit
    st.components.v1.html(xyzview._make_html(), width=width, height=height, scrolling=False)


def visualize_molecule(geometry, style='stick'):
    """Visualize the molecule using py3Dmol."""
    # Calculate number of atoms
    lines = geometry.splitlines()
    num_atoms = len(lines)

    # Reformat geometry to standard XYZ format
    xyz_format = f"{num_atoms}\n\n" + geometry

#    print("geometry debug")
#    print(xyz_format)  # Debugging output

    # Style options for visualization
    style_options = {
        'Stick': {'stick': {}},
        'Ball and Stick': {'stick': {}, 'sphere': {'radius': 0.5}},
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

        # 3D Visualization of selected geometry
        selected_geometry = cleaned_geometries[step]
        st.write(f'### Geometry at step {step}')
        visualize_molecule(selected_geometry)

if __name__ == "__main__":
    main()

