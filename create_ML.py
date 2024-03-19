def generate_streamlit_app(model_file_name, app_title):
    streamlit_code = f"""
\"\"\"
{app_title}
\"\"\"
import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem
import joblib

# Function to convert SMILES to Morgan fingerprints
def smiles_to_morgan(smiles, n_bits=2048):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=n_bits)
    return [int(b) for b in fp]

# Load the model
model = joblib.load("{model_file_name}")

st.title("{app_title}")
smiles_input = st.text_input("Enter a SMILES string:", "")

if st.button("Predict"):
    if smiles_input:
        fp = smiles_to_morgan(smiles_input)
        if fp is not None:
            prediction = model.predict([fp])[0]
            st.write(f"Predicted property: {{prediction}}")
        else:
            st.write("Invalid SMILES string. Please check your input.")
    else:
        st.write("Please enter a SMILES string.")
"""
    with open("streamlit_ML.py", "w") as file:
        file.write(streamlit_code)

def generate_requirements_txt():
    requirements = """
numpy
pandas
pyarrow
matplotlib
plotly
pillow
streamlit
rdkit-pypi
py3Dmol==2.0.0.post2
ipython_genutils
stmol
joblib
scikit-learn==1.2.2
"""
    with open("requirements.txt", "w") as file:
        file.write(requirements.strip())

def generate_packages_txt():
    packages = "libxrender1"
    with open("packages.txt", "w") as file:
        file.write(packages.strip())

if __name__ == "__main__":
    model_file_name = input("Enter the name of the model file (including extension): ")
    app_title = input("Enter a title for the App: ")
    generate_streamlit_app(model_file_name, app_title)
    generate_requirements_txt()
    generate_packages_txt()
    print("Streamlit app script and environment setup files generated successfully.")

