import os
import glob
import subprocess

# Ask the user for the title of the Streamlit app
app_title = input("Enter the title of your Streamlit app: ")

# Streamlit app code with dynamic app title
streamlit_app_code = f"""
import streamlit as st
import py3Dmol
import os
import glob

st.title('{app_title}')

# Function to visualize the cube file using py3Dmol with adjustable alpha value
def visualize_cube(cube_file_path, alpha):
    width = 400
    height = 480

    viewer = py3Dmol.view(width=width, height=height)

    with open(cube_file_path, 'r') as cube_file:
        cube_data = cube_file.read()
    viewer.addModel(cube_data, "cube")

    viewer.setStyle({{'stick': {{}}, 'sphere': {{'radius': 0.5}}}})
    viewer.addVolumetricData(cube_data, "cube", {{'isoval': 0.02, 'color': "blue", 'alpha': alpha}})
    viewer.addVolumetricData(cube_data, "cube", {{'isoval': -0.02, 'color': "red", 'alpha': alpha}})

    viewer.zoomTo()

    viewer_html = viewer._make_html()
    st.components.v1.html(viewer_html, width=width, height=height, scrolling=False)

# Streamlit interface setup
alpha = st.slider('Adjust the opacity for visualization:', min_value=0.0, max_value=1.0, value=0.65, step=0.05)

cubes_directory = 'cubes'
cube_files = glob.glob(os.path.join(cubes_directory, '*.cube'))
cube_files_names = [os.path.basename(cube_file) for cube_file in cube_files]

st.write("Visualizer")
selected_cube_file_name = st.selectbox('Select a cube file for the visualizer:', cube_files_names)

if st.button('Visualize Selected Orbital'):
    if selected_cube_file_name:
        selected_cube_file_path = os.path.join(cubes_directory, selected_cube_file_name)
        visualize_cube(selected_cube_file_path, alpha)
""".strip()

# Write the Streamlit app code to a file
streamlit_app_filename = "orbital_viewer_app.py"
with open(streamlit_app_filename, "w") as file:
    file.write(streamlit_app_code)

# Environment files content
packages_content = "libxrender1"
requirements_content = """
streamlit
py3Dmol==0.9.1
"""

# Write to packages.txt and requirements.txt
with open("packages.txt", "w") as file:
    file.write(packages_content)
with open("requirements.txt", "w") as file:
    file.write(requirements_content)

# Inform the user about the app and environment setup
print("Environment files have been created.")
print("\nYour app is ready to run!!\n")
print(f"Make sure you have Streamlit installed. If Streamlit isn't found, run: pip install streamlit")
print(f"Run the app: streamlit run {streamlit_app_filename}\n")

print("App creation complete!\n\nTo deploy your web app:")
print("Upload the directory to a new GitHub repository")
print("Sign in to https://streamlit.io/ with your GitHub account and locate your app's repo")
print("Name your app, this will be the URL!")
print("Click Deploy! and wait for the app to build")
print("Enjoy the app and share the link in your manuscripts and presentations.\n\n")

# Prompt the user to decide whether to run the web app now
user_decision = input("Do you want to run the web app now? (y/n): ").strip().lower()

if user_decision == 'y':
    # Attempt to run the Streamlit app
    try:
        subprocess.run(["streamlit", "run", streamlit_app_filename])
    except Exception as e:
        print(f"An error occurred while trying to run the app: {e}")
else:
    print(f"You can always start the app by running 'streamlit run {streamlit_app_filename}'")

