# SI Web App Generator: Create Web Apps for Your Computational SI Data

This guide provides a streamlined process to create and deploy web applications for sharing and visualizing computational Supplementary Information (SI) data, particularly useful for computational chemistry output files (current version supports Gaussian and Orca output files).

## Dependencies

- `subprocess`
- `os`
- `glob`
- `openbabel`
- `streamlit`

## Creating Your Web App

Follow these steps to create your web app:

1. **Prepare Data**: Collect all your Gaussian output files in one directory named `log_files`.
2. **Generate App Scripts**: Run `create.py`. This script will prompt you for a title for the web app. Upon completion, it generates new scripts: `stream_app.py`, `packages.txt`, and `requirements.txt`. The `stream_app.py` script is the core of your web app.
3. **Run Locally**: To preview your app locally, execute `streamlit run stream_app.py`. This command launches the app in your default web browser.

## Deploying Your Web App

To make your web app accessible online, follow these deployment steps:

1. **GitHub Repository**: Create a new repository on [GitHub](https://www.github.com) and move all necessary files there, excluding the `log_files` directory.
2. **Streamlit Signup**: Sign in to [Streamlit](https://streamlit.io/) using your GitHub credentials.
3. **Create New App**: Click on "New app", then select your app's repository, branch, and specify `stream_app.py` as the app's main file.
4. **Configure URL**: Assign a URL to your app.
5. **Deploy**: Click "Deploy" and wait for the app to build. This process may take several minutes.
6. **Share Your App**: Once deployment is complete, share the web app link with your coworkers and include it in your manuscripts and presentations for broader accessibility and impact.

By following these steps, you can efficiently create and deploy web applications to enhance the accessibility and visualization of your computational SI data.

