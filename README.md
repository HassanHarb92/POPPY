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

## Setup and Usage

### 1. Computational Chemistry Output Visualization

#### Step 1: Prepare Input Files

- Collect Gaussian or Orca output files.
- Store these in a directory named `log_files`.
- Ensure the files are error-free and properly formatted.

#### Step 2: Generate the Application

Run the script to create the app:

```
python create.py
```

You will be prompted to enter the title of your app. Upon execution, the following files are generated:

- `stream_app.py`: Core script for the Streamlit app.
- `packages.txt` and `requirements.txt`: Dependencies for the app.

#### Step 3: Preview Locally

Test the app locally by running:

```
streamlit run stream_app.py
```

Access the app in your default web browser to ensure proper functionality.

### 2. Molecular Orbital Visualization

#### Step 1: Prepare Cube Files

- Generate molecular orbital cube files using your preferred computational chemistry software (e.g., Gaussian, Orca).
- Place these files in a directory named `cubes`.

#### Step 2: Create the App

Run the script to generate the app:

```
python create_MO.py
```

Similar to the output visualization, this creates the necessary Streamlit scripts for deploying the app.

#### Step 3: Preview Locally

Run the app locally to confirm functionality:

```
streamlit run stream_MO.py
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
python create_db.py
```

Provide the title of the app and the path to the CSV file when prompted.

#### Step 3: Preview Locally

Test the dataset visualization app:

```
streamlit run streamlit_app_db.py
```

## Deployment

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

## Best Practices

### Code Quality

- Use tools like `pylint` or `flake8` to ensure clean and maintainable code.

- Example:

  ```
  pylint create.py
  ```

### Testing

- Test all scripts with a variety of sample datasets and file formats to identify bugs or compatibility issues.
- Automate testing with `pytest`.

### Error Handling

- Add robust error-handling mechanisms to manage missing files or unsupported formats gracefully.

### Documentation

- Provide detailed comments in all scripts.
- Maintain an updated `README.md` file in your GitHub repository with clear instructions and examples.

## Example Applications

### Visualizing Gaussian Output Files

This app showcases computed molecular structures, allowing users to:

- Select output files from a dropdown menu.
- Visualize 3D structures in various styles (e.g., stick, ball-and-stick).

### Molecular Orbital Viewer

Interactive visualization of molecular orbitals with adjustable opacity sliders for better analysis of electron density distributions.

### Dataset Explorer

Query and filter molecular datasets interactively based on user-defined parameters (e.g., property thresholds).

## Troubleshooting

### Common Issues

1. **App crashes on deployment**:
   - Ensure all dependencies are correctly listed in `requirements.txt`.
   - Check compatibility of Python libraries.
2. **Cube files not rendering**:
   - Confirm that cube files are properly formatted.
   - Ensure the `py3Dmol` library is correctly installed.
3. **Dataset app not loading CSV**:
   - Check that the CSV file follows the expected structure.
   - Validate data integrity to avoid parsing errors.

## Frequently Asked Questions (FAQ)

### Q1: Can POPPY handle non-Gaussian output files?

Yes. POPPY relies on Open Babel for file format conversion and can be extended to support various computational chemistry outputs.

### Q2: How do I customize the app layout?

Modify the `stream_app.py` script. Streamlit’s Markdown and widgets allow easy customization of titles, colors, and layouts.

### Q3: Can I integrate machine learning models with POPPY?

Yes. POPPY supports Python libraries like `scikit-learn` and `joblib` for training and deploying ML models within the app.

