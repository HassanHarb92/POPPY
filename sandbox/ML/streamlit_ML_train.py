import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
import numpy as np
import matplotlib.pyplot as plt

# Set up Streamlit app
st.title("Molecular Property Prediction App")

# Upload CSV file
uploaded_file = st.file_uploader("Upload your CSV file", type=["csv"])

if uploaded_file is not None:
    # Read the uploaded CSV
    data = pd.read_csv(uploaded_file)

    # Display a preview of the uploaded data
    st.write("Data preview:")
    st.write(data.head())

    # Ask for user input: SMILES column and the target column
    smiles_column = st.selectbox("Select the SMILES column", data.columns)
    target_column = st.selectbox("Select the target column for prediction (e.g., delta_H)", data.columns)

    # Function to compute Morgan fingerprints
    def compute_morgan_fingerprint(smiles, radius=2, nBits=2048):
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits)
            return np.array(fp)
        else:
            return np.zeros((nBits,))

    # Compute Morgan fingerprints for each SMILES string
    st.write("Computing Morgan fingerprints...")
    data['MorganFP'] = data[smiles_column].apply(compute_morgan_fingerprint)

    # Split the dataset into training and testing sets
    X = np.array(list(data['MorganFP']))
    y = data[target_column].values
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    # Train the regression model
    st.write("Training RandomForestRegressor model...")
    model = RandomForestRegressor(random_state=42)
    model.fit(X_train, y_train)

    # Predict on the test set
    y_pred = model.predict(X_test)

    # Calculate evaluation metrics
    mae = mean_absolute_error(y_test, y_pred)
    rmsd = np.sqrt(mean_squared_error(y_test, y_pred))
    r_squared = r2_score(y_test, y_pred)

    # Display the metrics
    st.write(f"Mean Absolute Error (MAE): {mae}")
    st.write(f"Root Mean Square Deviation (RMSD): {rmsd}")
    st.write(f"R-squared: {r_squared}")

    # Plot actual vs predicted delta_H
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(y_test, y_pred, alpha=0.5)
    ax.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'k--', lw=2)
    ax.set_title("Actual vs Predicted")
    ax.set_xlabel("Actual delta_H")
    ax.set_ylabel("Predicted delta_H")

    # Show the plot in Streamlit
    st.pyplot(fig)

