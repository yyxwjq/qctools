from qctools.element_tools import get_elements
from qctools.rdf import get_rdf
from ase.io import read

a = read('test.xyz', index=':')
get_rdf(a, cutoff=7.0, bin_size=0.05, first_neighbor=True, cores=4)
get_rdf(a, cutoff=7.0, bin_size=0.05, first_neighbor=False, cores=4)
