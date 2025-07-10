from qctools.adf import get_adf
from ase.io import read

a = read('test.xyz', index=':')
get_adf(a, rcut=3.0, bin_size=5, cores=4)