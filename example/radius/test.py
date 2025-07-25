from qctools.radius import MinimumEnclosingBall
from ase.io import read
import numpy as np

image = read('POSCAR')
saved_elements = ['H', 'O']
elements_dict = image.symbols.indices()
newidx = np.concatenate([elements_dict[element] for element in saved_elements if element in elements_dict])
points = image[newidx].positions

meb = MinimumEnclosingBall(points)

welzl_sphere, welzl_history = meb.stable_welzl()
welzl_sphere_exact, welzl_history_exact = meb.stable_welzl(exact=True)
ritter_sphere = meb.ritter()
bb_sphere, bb_history = meb.stable_bouncing_bubble()

print("\n=== Ritter Algorithm (Approximate) ===")
print(welzl_sphere)
print(welzl_sphere_exact)
print("\n=== Ritter Algorithm (Approximate) ===")
print(ritter_sphere)
print(f"\n=== Bouncing Bubble Test with custom parameters ===")
print(bb_sphere)

