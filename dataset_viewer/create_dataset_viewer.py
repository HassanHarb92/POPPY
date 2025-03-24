import subprocess
import os
import glob


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

def create_streamlit_app_script(title, df_filename):
    app_script_content = f"""import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol

# Load the CSV file into a DataFrame
df = pd.read_csv('{df_filename}')

def SmilesTo3DViewer(molecule_smiles):
    mol = Chem.MolFromSmiles(molecule_smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.MMFFOptimizeMolecule(mol)
    xyz = Chem.MolToXYZBlock(mol)
    return xyz

def display_molecule(xyz, style_options, selected_style):
    scale = 1
    width = int(640.0 * scale)
    height = int(480.0 * scale)
    xyzview = py3Dmol.view(width=width, height=height)
    xyzview.addModel(xyz, 'xyz')
    xyzview.setStyle(style_options[selected_style])
    xyzview.zoomTo()
    xyzview.show()
    st.components.v1.html(xyzview._make_html(), width=width, height=height, scrolling=False)

# Streamlit UI
st.title("{title}")
# Dropdown to select a molecule name
selected_name = st.selectbox('Select a molecule', df['Names'])

# Get the SMILES string for the selected molecule
selected_smiles = df.loc[df['Names'] == selected_name, 'SMILES'].iloc[0]
# Style options
style_options = {{
    'Stick': {{'stick': {{}}}},
    'Ball and Stick': {{'stick': {{}}, 'sphere': {{'radius': 0.5}}}},
    'Spacefill': {{'sphere': {{}}}},
}}
selected_style = st.radio('Select visualization style', list(style_options.keys()))

# Visualization button
if st.button("Visualize"):
    # Display additional properties if they exist
    st.markdown("### SMILES String:")
    st.markdown(selected_smiles)
    if len(df.columns) > 2:
        st.write("### Molecule Properties")
        # Retrieve and display the properties for the selected molecule
        properties_df = df.loc[df['Names'] == selected_name, df.columns[2:]]
        properties = properties_df.iloc[0].to_dict()
        # Create a two-column layout for a neat display of properties
        col1, col2 = st.columns(2)
        for index, (prop, value) in enumerate(properties.items()):
            with col1 if index % 2 == 0 else col2:
                st.metric(label=prop, value=value)

    xyz_content = SmilesTo3DViewer(selected_smiles)
    display_molecule(xyz_content, style_options, selected_style)
"""

    return app_script_content

def main():
    title = input("Enter the title of the Streamlit app: ")
    df_filename = input("Enter the filename of the DataFrame (CSV): ")

    script_content = create_streamlit_app_script(title, df_filename)

    output_filename = "dataset_viewer_app.py"
    with open(output_filename, "w") as file:
        file.write(script_content)
    
    print(f"Streamlit app script has been created: {output_filename}")


# Write to packages.txt
    with open("packages.txt", "w") as file:
        file.write(packages_content)
    
    # Write to requirements.txt
    with open("requirements.txt", "w") as file:
        file.write(requirements_content)
    
    print("Environment files have been created.")
    
    print("Your app is ready to run!!")
    print("Make sure you have Streamlit installed. if streamlit isnt found run: pip install streamlit")
    print("Run the app: streamlit run streamlit_app_db.py ")
    
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
            subprocess.run(["streamlit", "run", "streamlit_app_db.py"])
        except Exception as e:
            print(f"An error occurred while trying to run the app: {e}")
    else:
        # Provide information on how to run the app later
        print("You can always start the app by running 'streamlit run stream_app.py'")
    
    
if __name__ == "__main__":
    main()







