import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol

# Load the CSV file into a DataFrame
df = pd.read_csv('')

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
st.title("test DB")
# Dropdown to select a molecule name
selected_name = st.selectbox('Select a molecule', df['Names'])

# Get the SMILES string for the selected molecule
selected_smiles = df.loc[df['Names'] == selected_name, 'SMILES'].iloc[0]
# Style options
style_options = {
    'Stick': {'stick': {}},
    'Ball and Stick': {'stick': {}, 'sphere': {'radius': 0.5}},
    'Spacefill': {'sphere': {}},
}
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
