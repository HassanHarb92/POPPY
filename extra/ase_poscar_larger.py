from ase.io import read, write
from ase.build import surface, make_supercell
from ase.visualize import view

# Step 1: Read the POSCAR file
bulk = read('POSCAR')

# Step 2: Create a supercell
# This example creates a 2x2x1 supercell of the original unit cell.
# Adjust the multipliers (2, 2, 1 here) as needed for your application.
supercell = make_supercell(bulk, [[2, 0, 0], [0, 2, 0], [0, 0, 1]])

# Step 3: Generate the surface
# Specify the Miller indices (hkl) of the surface and the number of layers to keep in the z direction.
# The 'vacuum' parameter adds vacuum padding above the surface (in Ångströms).
surface = surface(supercell, (1, 1, 1), layers=3, vacuum=10)

# Optional: Save the surface to a new POSCAR file or visualize it
write('surface_POSCAR', surface, format='vasp')

