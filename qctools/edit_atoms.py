#!/usr/bin/env python3
from ase.io import read, write
from ase import Atoms
import numpy as np
import sys
import argparse


def remove(f, elements):
    """Remove specified elements from atomic structure.
    
    Parameters:
    -----------
    f : str
        Input file path
    elements : list
        List of element symbols to remove
    """
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

def main():
    """Main entry point for qctools-edit command line tool."""
    parser = argparse.ArgumentParser(
        description='QCTools structure editing utility',
        prog='qctools-edit'
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Remove command
    remove_parser = subparsers.add_parser('remove', help='Remove elements from structure')
    remove_parser.add_argument('input_file', help='Input structure file')
    remove_parser.add_argument('elements', nargs='+', help='Elements to remove')
    
    if len(sys.argv) == 1:
        parser.print_help()
        return
    
    args = parser.parse_args()
    
    if args.command == 'remove':
        remove(args.input_file, args.elements)
    else:
        parser.print_help()

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == 'remove':
        remove(sys.argv[2], sys.argv[3:])
    else:
        main()