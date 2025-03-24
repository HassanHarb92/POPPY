import streamlit as st
import os

# Function to list all .xyz files in a directory
def list_xyz_files(directory):
    return [f for f in os.listdir(directory) if f.endswith('.xyz')]

# Function to read and return the content of an XYZ file
def read_xyz_file(xyz_file_path):
    with open(xyz_file_path, 'r') as file:
        return file.read()

# Streamlit app starts here
st.title('Molecular Visualizer')
st.markdown('## Visualize Molecules from XYZ Files')

# Directory containing the .xyz files
xyz_files_directory = 'xyz_files'

# List all .xyz files in the directory
xyz_files = list_xyz_files(xyz_files_directory)

if xyz_files:
    # Dropdown to select an XYZ file
    selected_file = st.selectbox('Select an XYZ file:', xyz_files)

    # Style options as radio buttons
    style_options = {
        'Stick': {'stick': {}},
        'Ball and Stick': {'stick': {}, 'sphere': {'radius': 0.5}},
        'Spacefill': {'sphere': {}}
    }
    selected_style = st.radio('Select visualization style:', list(style_options.keys()))

    # Button to visualize the selected XYZ file
    if st.button('Visualize'):
        # Full path to the selected file
        full_path_to_file = os.path.join(xyz_files_directory, selected_file)

        # Read the selected XYZ file
        xyz_content = read_xyz_file(full_path_to_file)

        # Convert the style selection to a JavaScript object
        style_js = str(style_options[selected_style]).replace("'", "\"")

        # Define the JavaScript and HTML for embedding the 3D visualization
        html_template = f"""
        <div id="moldiv" style="height: 400px; width: 600px;"></div>
        <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
        <script>
            var viewer = $3Dmol.createViewer("moldiv", {{backgroundColor: "white"}});
            viewer.addModel(`{xyz_content}`, "xyz");
            viewer.setStyle({{ }}, {style_js});
            viewer.zoomTo();
            viewer.render();
        </script>
        """

        # Display the visualization in Streamlit
        st.components.v1.html(html_template, height=400)
else:
    st.error("No XYZ files found in the specified directory.")

