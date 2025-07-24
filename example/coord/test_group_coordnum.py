from ase.io import read
from qctools.coord import group_coordnum
a = read('test.dump', ':', format='lammps-dump-text')
print(len(a))

coord = group_coordnum(a, [i for i in range(1372)], [i for i in range(1372, len(a[0]))], r0=3.0, tolerance=0.001, cores=4)
print(coord)
