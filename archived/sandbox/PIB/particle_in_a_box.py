import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import hbar, pi

def main():
    # Sidebar for parameter inputs
    st.sidebar.title("Parameters")
    m = st.sidebar.number_input("Mass of Particle (m)", value=1.0, step=0.1)
    l = st.sidebar.number_input("Length of Box (l)", value=1.0, step=0.1)
    n_levels = st.sidebar.slider('Number of Energy Levels', 1, 10, 5)
    
    # Main content
    st.title('Quantum Particle in a Box Visualization')
    
    def wavefunction(n, x, l):
        """Calculate the wavefunction for a given level, position, and box length."""
        return np.sqrt(2 / l) * np.sin(n * pi * x / l)
    
    def energy_level(n, m, l):
        """Calculate the energy level for a given level, mass, and box length."""
        return n**2 * pi**2 * hbar**2 / (2 * m * l**2)
    
    # Plotting code
    fig, ax1 = plt.subplots()
    
    x = np.linspace(0, l, 1000)  # x-values within the box
    
    # Only one axis is needed since we're aligning everything according to energy levels
    ax1.set_xlabel('Position (x)')
    ax1.set_ylabel('Energy / Wave Amplitude')
    
    colors = plt.cm.viridis(np.linspace(0, 1, n_levels))
    
    for n, color in zip(range(1, n_levels + 1), colors):
        psi_n = wavefunction(n, x, l)
        energy = energy_level(n, m, l)
        
        # Normalize wavefunction for visualization purposes
        psi_n_normalized = psi_n / np.max(np.abs(psi_n)) * energy / 10  # Adjust scaling factor as needed
        
        # Offset wavefunction to start at its corresponding energy level
        psi_n_offset = psi_n_normalized + energy
        
        # Plot the energy level as a horizontal line
        ax1.hlines(energy, 0, l, colors='gray', linestyles='dashed', label=f'Energy n={n}' if n == 1 else "")
    
        # Plot wavefunction
        ax1.plot(x, psi_n_offset, label=f'Wavefunction n={n}', color=color)
        
        # Plot probability density (squared wavefunction), also offset
        probability_density_offset = psi_n**2 / np.max(np.abs(psi_n**2)) * energy / 10 + energy  # Adjust scaling
        ax1.fill_between(x, energy, probability_density_offset, color=color, alpha=0.3)
    
    # Adjust layout to make room for the legend
    #ax1.legend()
    ax1.grid(True)
    plt.tight_layout()
    
    st.pyplot(fig)
 

if __name__ == "__main__":
    main() 
