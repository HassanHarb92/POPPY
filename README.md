# POPPY: Code Documentation

Last updated March 24, 2025.

## Introduction

POPPY (Python-based Output Processing and Presentation for chemistrY) is an open-source tool designed to facilitate the creation of interactive web applications for computational chemistry data. Leveraging Python and popular libraries like Streamlit, Open Babel, and Py3Dmol, POPPY enables users to visualize and share molecular structures, datasets, and properties seamlessly across devices using just a web browser. This documentation provides a step-by-step guide to set up, customize, and deploy applications created with POPPY.

## System Requirements

### Software

- **Python**: Version 3.8 or later.
- Required Libraries:
  - `streamlit`
  - `openbabel`
  - `py3Dmol`
  - `matplotlib`
  - `rdkit`
  - `joblib`
  - `numpy`
  - `pandas`
  - `pyarrow`
  - `plotly`
  - `pillow`
  - `ipython_genutils`
  - `stmol`
  - `scikit-learn==1.2.2`
  - `pubchempy`

Install dependencies by running:

```
pip install -r requirements.txt
```

### Hardware

POPPY applications are browser-based and do not require specific hardware. Any modern device with an internet connection and a web browser can be used to access the deployed apps.

## Features Overview

POPPY supports three primary application types:

1. **Computational Chemistry Output Visualization**: Processes output files (e.g., Gaussian, Orca) to render molecular structures and computed data.
2. **Molecular Orbital Visualization**: Displays 3D visualizations of molecular orbitals generated as cube files.
3. **Dataset Visualization**: Enables interactive exploration of molecular datasets in CSV format.

**Details on the scripts are presented in this table:**

| **Folder Name**      | **Contents**                                          | **Purpose**                                                  |
| -------------------- | ----------------------------------------------------- | ------------------------------------------------------------ |
| `molecule_viewer/`   | `create_molecule_viewer.py`, `molecule_viewer_app.py` | Visualize molecular structures                               |
| `orbital_viewer/`    | `create_orbital_viewer.py`, `orbital_viewer_app.py`   | Visualize molecular orbitals                                 |
| `dataset_viewer/`    | `create_dataset_viewer.py`, `dataset_viewer_app.py`   | Browse and explore datasets                                  |
| `ml_training_app/`   | `ml_training_app.py`                                  | Train machine learning models                                |
| `ml_deployment_app/` | `ml_deployment_app.py`                                | Deploy trained machine learning models                       |
| `particle_in_box/`   | `particle_in_a_box.py`                                | Solve and visualize particle-in-a-box                        |
| `irc_visualizer/`    | `IRC_visualizer_app.py`                               | Visualize Intrinsic Reaction Coordinates (from a Gaussian output file) |
| `xtb_wrapper/`       | `xtb2Go.py`                                           | xTB calculation launcher                                     |



## Setup and Usage

### 1. Computational Chemistry Output Visualization

#### Step 1: Prepare Input Files

- Collect Gaussian or Orca output files.
- Store these in a directory named `log_files`.
- Ensure the files are error-free and properly formatted.

#### Step 2: Generate the Application

Run the script to create the app:

```
python create_molecule_viewer.py
```

You will be prompted to enter the title of your app. Upon execution, the following files are generated:

- `molecule_viewer_app.py`: Core script for the Streamlit app.
- `packages.txt` and `requirements.txt`: Dependencies for the app.

#### Step 3: Preview Locally

Test the app locally by running:

```
streamlit run molecule_viewer_app.py
```

Access the app in your default web browser to ensure proper functionality.

### 2. Molecular Orbital Visualization

#### Step 1: Prepare Cube Files

- Generate molecular orbital cube files using your preferred computational chemistry software (e.g., Gaussian, Orca).
- Place these files in a directory named `cubes`.

#### Step 2: Create the App

Run the script to generate the app:

```
python create_orbital_viewer.py
```

Similar to the output visualization, this creates the necessary Streamlit scripts for deploying the app.

#### Step 3: Preview Locally

Run the app locally to confirm functionality:

```
streamlit run orbital_viewer_app.py
```

### 3. Dataset Visualization

#### Step 1: Format the Dataset

Create a CSV file with the following columns:

- **Name**: Descriptive name for each molecule.
- **SMILES**: Simplified molecular input line entry system representation.
- Additional columns for molecular properties as needed.

#### Step 2: Create the App

Run the script to create the dataset visualization app:

```
python create_dataset_viewer.py
```

Provide the title of the app and the path to the CSV file when prompted.

#### Step 3: Preview Locally

Test the dataset visualization app:

```
streamlit run dataset_viewer_app.py
```

## Deployment (optional)

POPPY apps can be deployed on platforms like Streamlit Cloud for easy sharing and accessibility.

### Steps for Deployment

1. **Prepare Repository**:
   - Create a GitHub repository and upload the necessary files (e.g., `stream_app.py`, `requirements.txt`).
   - Exclude large directories like `log_files` or `cubes` to keep the repository lightweight.
2. **Deploy on Streamlit Cloud**:
   - Sign in to [Streamlit Cloud](https://streamlit.io/).
   - Click on "New App" and select your GitHub repository.
   - Specify the main script file (e.g., `stream_app.py`).
3. **Finalize Deployment**:
   - Assign a custom URL to the app.
   - Share the URL with collaborators or include it in publications and presentations.

