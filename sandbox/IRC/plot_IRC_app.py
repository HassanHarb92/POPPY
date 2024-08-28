import streamlit as st
import re
import matplotlib.pyplot as plt

def parse_log_file(file):
    """
    Parse the Gaussian .log file to extract geometries and energies.
    """
    with open(file, 'r') as f:
        lines = f.readlines()
    
    geometries = []
    energies = []

    current_geometry = []
    reading_geometry = False
    energy_pattern = re.compile(r'SCF Done:.*=(-?\d+\.\d+)')

    for line in lines:
        if 'Input orientation:' in line:
            reading_geometry = True
            current_geometry = []
            continue
        
        if reading_geometry:
            if '----' in line or 'Distance matrix' in line:
                reading_geometry = False
                geometries.append(current_geometry)
                continue
            if line.strip() and len(line.split()) == 6:  # Check for geometry line
                current_geometry.append(line.split()[3:6])

        energy_match = energy_pattern.search(line)
        if energy_match:
            energies.append(float(energy_match.group(1)))

    return geometries, energies

def main():
    st.title('Gaussian IRC Visualization App')

    uploaded_file = st.file_uploader("Upload a Gaussian IRC .log file", type="log")

    if uploaded_file is not None:
        geometries, energies = parse_log_file(uploaded_file)
        num_geometries = len(geometries)

        if num_geometries % 2 == 0:
            st.error("The .log file does not seem to have an odd number of geometries (TS in the middle). Please check the file.")
            return

        # Calculate mid index for TS placement
        ts_index = num_geometries // 2

        # Plot energies
        fig, ax = plt.subplots()
        ax.plot(range(-ts_index, ts_index + 1), energies, marker='o')
        ax.set_xlabel('Step')
        ax.set_ylabel('Energy (Hartree)')
        ax.set_title('Energy Profile along IRC')
        st.pyplot(fig)

        # Slider for geometries
        step = st.slider('Select IRC Step', -ts_index, ts_index, 0)

        # Display selected geometry
        selected_geometry = geometries[ts_index + step]
        st.write(f'### Geometry at step {step} (Relative to TS)')
        st.text("\n".join([" ".join(atom) for atom in selected_geometry]))

if __name__ == "__main__":
    main()

