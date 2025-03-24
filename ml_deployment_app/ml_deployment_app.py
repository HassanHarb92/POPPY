
"""
ML Hydrogenation test
"""
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
model = joblib.load("Hydrogenation_random_forest_regressor_model.joblib")

st.title("ML Hydrogenation test")
smiles_input = st.text_input("Enter a SMILES string:", "")

if st.button("Predict"):
    if smiles_input:
        fp = smiles_to_morgan(smiles_input)
        if fp is not None:
            prediction = model.predict([fp])[0]
            st.write(f"Predicted property: {prediction}")
        else:
            st.write("Invalid SMILES string. Please check your input.")
    else:
        st.write("Please enter a SMILES string.")
