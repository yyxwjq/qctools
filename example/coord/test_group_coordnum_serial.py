# -*- coding: utf-8 -*-
from ase.io import read
from qctools.coord import _group_coordnum_serial

a = read('00.vasp')
coord_num = _group_coordnum_serial(a, 
                           group1=[i for i in range(1372)], 
                           group2=[i for i in range(1372, len(a))], 
                           r0=3.0,
                           tolerance=0.001)
print(coord_num)
