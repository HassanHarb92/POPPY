import subprocess
import os
import glob


def convert_log_to_xyz(log_file_path, xyz_file_path):
    """
    Uses Open Babel to convert a log file to XYZ format.
    Differentiates between Gaussian and ORCA output files based on their content.
    """
    # Read the content of the log file
    with open(log_file_path, 'r') as file:
        content = file.read()
    
    # Check for Gaussian output file
    if "Gaussian, Inc." in content:
        command = ["obabel", log_file_path, "-O", xyz_file_path]
    # Check for ORCA output file
    elif "Max Planck Institute fuer Kohlenforschung" in content:
        command = ["obabel", "-i", "orca", log_file_path, "-o", "xyz", "-O", xyz_file_path]

    elif "POSCAR" or "CONTCAR" in log_file_path:
        command = ['obabel','-i','vasp',log_file_path,'-o','xyz','-O',xyz_file_path]
    else:
        raise ValueError("Unknown file type. Cannot determine if the file is Gaussian or ORCA output.")

    # Execute the command
    subprocess.run(command, check=True)

# Ensure the xyz_files directory exists
xyz_dir = 'xyz_files'
if not os.path.exists(xyz_dir):
    os.makedirs(xyz_dir)

log_dir = 'log_files'
for log_file in glob.glob(os.path.join(log_dir, '*')):
    # Generate the XYZ file name by changing the extension to .xyz
    file_root, _ = os.path.splitext(os.path.basename(log_file))
    xyz_file_name = f"{file_root}.xyz"
    xyz_file_path = os.path.join(xyz_dir, xyz_file_name)

    # Convert the log file to XYZ format
    convert_log_to_xyz(log_file, xyz_file_path)

print("XYZ files have been generated in the 'xyz_files' directory.")

# Ask the user for the title of the project/manuscript
project_title = input("Enter the title of the project/manuscript: ")

# Content of the stream_app.py file with a placeholder for the project title
stream_app_content = f"""
import streamlit as st
import os
import py3Dmol

# Function to list all .xyz files in a directory
def list_xyz_files(directory):
    return [f for f in os.listdir(directory) if f.endswith('.xyz')]

# Function to read and return the content of an XYZ file
def read_xyz_file(xyz_file_path):
    with open(xyz_file_path, 'r') as file:
        return file.read()

# Set the page to wide mode
#st.set_page_config(layout="wide")

# Streamlit app starts here
st.title('{project_title}')
st.markdown('## Supplementary Information Visualizer')

# Directory containing the .xyz files
xyz_files_directory = 'xyz_files'

# List all .xyz files in the directory
xyz_files = list_xyz_files(xyz_files_directory)

# Dropdown to select an XYZ file
selected_file = st.selectbox('Select an XYZ file', xyz_files)

# Style options as radio buttons
style_options = {{
    'Stick': {{'stick': {{}}}},
    'Ball and Stick': {{'stick': {{}}, 'sphere': {{'radius': 0.5}}}},
    'Spacefill': {{'sphere': {{}}}}
}}
selected_style = st.radio('Select visualization style', list(style_options.keys()))

# Button to visualize the selected XYZ file
if st.button('Visualize'):
    # Full path to the selected file
    full_path_to_file = os.path.join(xyz_files_directory, selected_file)

    # Read the selected XYZ file
    xyz_content = read_xyz_file(full_path_to_file)
    scale = 1

    width = int(640.0*scale)
    height = int(480.0*scale)

    # Visualize the molecule using py3Dmol
    xyzview = py3Dmol.view(width=width, height=height)
    xyzview.addModel(xyz_content, 'xyz')
    xyzview.setStyle(style_options[selected_style])  # Use the selected style
    xyzview.zoomTo()

    # Display the visualization in Streamlit
    xyzview.show()
    st.components.v1.html(xyzview._make_html(), width=width, height=height, scrolling=False)
"""

# Write the content to stream_app.py file
with open("molecule_viewer_app.py", "w") as file:
    file.write(stream_app_content)

print("stream_app.py has been generated.")

# Define the content for packages.txt
packages_content = "libxrender1"

# Define the content for requirements.txt
requirements_content = """\
# add the required pypi packages here that your app actually needs
numpy
pandas
pyarrow
matplotlib
plotly
pillow
streamlit
rdkit
# due to changes of py3Dmol, we need to install a specific version:
py3Dmol==2.0.0.post2
# required for stmol:
ipython_genutils
stmol
"""

# Write to packages.txt
with open("packages.txt", "w") as file:
    file.write(packages_content)

# Write to requirements.txt
with open("requirements.txt", "w") as file:
    file.write(requirements_content)

print("Environment files have been created.")

print("Your app is ready to run!!")
print("Make sure you have Streamlit installed. if streamlit isnt found run: pip install streamlit")
print("Run the app: streamlit run stream_app.py ")

print("\n\n")
print("App creation complete!\n\n")
print("To deploy your web app:") 
print("Upload the directory to a new Github repository")
print("Sign in to https://streamlit.io/ with your github account and locate your app's repo")
print("Name your app, this will be the URL!")
print("Click Deploy! and wait for the app to build")
print("Enyoy the app and share the link in your manuscripts and presentations.\n\n")

# Prompt the user to decide whether to run the web app now
user_decision = input("Do you want to run the web app now? (y/n): ").strip().lower()

if user_decision == 'y':
    # Run the Streamlit app
    try:
        subprocess.run(["streamlit", "run", "stream_app.py"])
    except Exception as e:
        print(f"An error occurred while trying to run the app: {e}")
else:
    # Provide information on how to run the app later
    print("You can always start the app by running 'streamlit run stream_app.py'")

