# Web App Generator for Computational Chemistry Data Visualization

This guide provides a streamlined process for creating and deploying web applications. These apps are designed to share and visualize and share computational chemistry data, including Gaussian and Orca output files, molecular orbital data from cube files, and molecular datasets via CSV files. The guide facilitates the visualization of computational chemistry output files and molecular orbitals, enhancing the accessibility and impact of your research findings.

## Dependencies

- `subprocess`
- `os`
- `glob`
- `openbabel` 
- `streamlit`

## Creating Your Web App

### For Computational Chemistry Output Files

1. **Prepare Data**: Collect all your output files in a directory named `log_files`.
2. **Generate App Scripts**: Run `create.py`. This script prompts you for a title for the web app and generates `stream_app.py`, `packages.txt`, and `requirements.txt`. `stream_app.py` is the core of your computational chemistry data visualization web app.

### For Molecular Orbital Visualization

1. **Generate Cube Files**: Use your computational chemistry software to generate cube files containing molecular orbital data. To create cube files with Gaussian, check this [Gaussian CubeGen documentation](https://gaussian.com/cubegen/). 
2. **Organize Cube Files**: Place all cube files into a directory named `cubes`.
3. **Generate App Scripts**: Execute `create_MO.py`. Similar to `create.py`, it asks for a title, then generates `stream_MO.py`, `packages.txt`, and `requirements.txt`, with `stream_MO.py` serving as the core of your molecular orbital visualization web app.

### For Visualizing Molecular Datasets

1. **Prepare Your CSV File**: Ensure your CSV is formatted correctly: the first column should be named "Name" and contain the molecule names; the second column, "SMILES", should contain the SMILES strings of the molecules; subsequent columns can include any properties of interest related to the molecules.
2. **Generate App Scripts**: Run `create_db.py`. This script prompts you for a title for the web app and the location of the dataset. Following these inputs, it automatically generates `streamlit_app_db.py`, the primary file for your web app, alongside `requirements.txt` and `packages.txt`, ensuring your app has all it needs to run effectively.

## Running Your Web App Locally

- **Preview Locally**: To preview either web app locally, execute the appropriate Streamlit command:
  - For computational chemistry output files: `streamlit run stream_app.py`
  - For molecular orbital data: `streamlit run stream_MO.py`
  

This command launches the selected app in your default web browser for interactive visualization.

## Deploying Your Web App

The deployment process is identical for both types of web apps and can be summarized in the following steps:

1. **GitHub Repository**: Create a new GitHub repository and upload the necessary files. Exclude the `log_files` or `cubes` directory to keep your repository streamlined.
2. **Streamlit Signup**: Sign in to [Streamlit](https://streamlit.io/) using your GitHub credentials.
3. **Create New App**: Navigate to "New app", then select your repository, branch, and specify the main file (`stream_app.py` or `stream_MO.py`).
4. **Configure URL**: Choose a URL for your web app.
5. **Deploy**: Initiate the deployment process by clicking "Deploy". Wait for the build to complete, which may take a few minutes.
6. **Share Your App**: Once deployed, share the web app link with colleagues and include it in your manuscripts and presentations to broaden the accessibility and impact of your research.

By following these steps, you can efficiently create and deploy web applications tailored to the visualization needs of your computational chemistry research, making your supplementary information more accessible and interactive.

*Note*: To create a supercell from a POSCAR file, check extra/ase_poscar_larger.py 

