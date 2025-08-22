#!/usr/bin/env python3
from ase.io import read, write
from ase import Atoms
import numpy as np
import sys


def remove(f, elements):

    try:
        atoms = read(f)
    except Exception as e:
        print(f"Error: Can't read {f}")
        print(f"[Error]: {str(e)}")
        sys.exit(1)
    
    elements_set = set(atoms.symbols)
    elements_dict = atoms.symbols.indices()
    saved_elements = elements_set - set(elements)


    newidx = np.concatenate([elements_dict[element] for element in saved_elements if element in elements_dict])
    new_atoms = atoms[newidx]
    
    new_atoms.calc = atoms.calc

    out_name = f.rsplit('.', 1)[0]
    new_file = f"{out_name}_remove.vasp"
    
    is_direct = 'Direct' in atoms.info.get('labels', ['Direct'])[0]

    write(new_file, new_atoms, format='vasp', vasp5=True, direct=is_direct)
    
    print(f"Generate new file: {new_file}")
    print(f"Deleted elements: {', '.join(elements)}")
    print(f"Saved elements: {', '.join(set(new_atoms.symbols))}")
    print(f"Saved atomic number: {len(new_atoms)}")

if __name__ == "__main__":
    if sys.argv[1] == 'remove':
        remove(sys.argv[2], sys.argv[3:])