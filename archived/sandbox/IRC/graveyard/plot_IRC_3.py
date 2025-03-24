import streamlit as st
import matplotlib.pyplot as plt
import py3Dmol

# Function to parse the Gaussian log file
def parse_gaussian_log(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    geometries = []
    energies = []
    temp_geom = []
    capture_geom = False

    for line in lines:
        if "Input orientation:" in line or "Standard orientation:" in line:
            capture_geom = True
            # Skip the header lines of the geometry section
            header_lines_skipped = 0
            continue

        if capture_geom:
            if "---------------------------------------------------------------------" in line:
                if temp_geom:  # Ensure we don't add empty geometries
                    # Join the geometry data into a single string and add to geometries
                    geometries.append("\n".join(temp_geom))
                    temp_geom = []  # Reset for the next geometry
                capture_geom = False
                continue
            # Ensure we skip the header lines of the geometry section
            if header_lines_skipped < 5:
                header_lines_skipped += 1
                continue
            # Capture the geometry data
            parts = line.strip().split()
            if len(parts) == 6:  # Ensure this line is part of the geometry data
                atom_index, atom_type, x, y, z = parts[1], parts[2], parts[3], parts[4], parts[5]
                temp_geom.append(f"{atom_type} {x} {y} {z}")

        if "SCF Done:" in line:
            parts = line.split()
            energy = parts[4]  # Assuming the energy value is always in this position
            energies.append(float(energy))

    return geometries, energies



# Function to plot energies
def plot_energies(energies):
    fig, ax = plt.subplots()
    ax.plot(energies, marker='o', linestyle='-')
    ax.set_xlabel('Step')
    ax.set_ylabel('Energy (a.u.)')
    return fig

# Function to display molecular structure
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import base64
from io import BytesIO

def show_structure(geometry):
    # Convert the geometry string to an RDKit molecule
    mol = Chem.MolFromXYZBlock(geometry, sanitize=False)
    if mol:
        # Compute 2D coordinates for better visualization
        AllChem.Compute2DCoords(mol)
        # Use RDKit to draw the molecule to an image
        img = Draw.MolToImage(mol, size=(300, 300))
        # Convert the image to a format that Streamlit can display
        buffered = BytesIO()
        img.save(buffered, format="PNG")
        img_str = base64.b64encode(buffered.getvalue()).decode()
        return f'<img src="data:image/png;base64,{img_str}"/>'
    else:
        return "Unable to generate structure."

# Update the part in the Streamlit app that displays the structure:
if geometries:  # Ensure there are geometries to work with
    point_to_view = st.slider("Select a point to view", 0, len(geometries)-1)
    st.write("Molecular Structure at Selected Point")
    html_str = show_structure(geometries[point_to_view])
    st.markdown(html_str, unsafe_allow_html=True)
else:
    st.write("No geometries to display. Please check the log file and parsing logic.")

# Streamlit app
def main():
    st.title('Gaussian IRC Visualization')

    file_path = "H-N-IRC.log"  # Directly specifying the file path
    geometries, energies = parse_gaussian_log(file_path)
    
    # Plotting energies
    st.write("Energy Profile")
    fig = plot_energies(energies)
    st.pyplot(fig)

    # Selecting a point to visualize - make sure to adjust the slider's range dynamically
    if geometries:  # Ensure there are geometries to work with
        point_to_view = st.slider("Select a point to view", 0, len(geometries)-1)
        st.write("Molecular Structure at Selected Point")
        view = show_structure(geometries[point_to_view])
        st.write(view.show())
    else:
        st.write("No geometries to display. Please check the log file and parsing logic.")


    
    # Selecting a point to visualize
    point_to_view = st.slider("Select a point to view", 0, len(geometries)-1, 0)
    
    # Visualizing the structure
    st.write("Molecular Structure at Selected Point")
    view = show_structure(geometries[point_to_view])
    st.write(view.show())

if __name__ == "__main__":
    main()

