import streamlit as st
import py3Dmol

st.title('Simple Py3Dmol Visualization Test')

# A simple molecule in XYZ format for testing
xyz_content = """
3

O 0.0 0.0 0.0
H 0.0 0.0 1.0
H 1.0 0.0 0.0
"""

# Py3Dmol visualization
xyzview = py3Dmol.view(width=300, height=300)
xyzview.addModel(xyz_content, 'xyz')
xyzview.setStyle({'stick': {}})
xyzview.zoomTo()

# Attempt to display the visualization
html_str = xyzview.show()
st.components.v1.html(html_str, width=300, height=300, scrolling=False)

