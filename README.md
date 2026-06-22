# QCTools - Quantum Chemistry Analysis Toolkit

[简体中文](README.zh-CN.md)

QCTools is a small Python toolkit for routine quantum chemistry and atomistic
simulation analysis. It currently focuses on trajectory structure analysis,
machine-learning potential error inspection, simple structure editing, and a
shell helper for monitoring VASP job folders.

Current version: `0.1.0`

## Installation

Install the development version from GitHub:

```bash
git clone https://github.com/yyxwjq/qctools.git
cd qctools
pip install -e .
```

Or install directly from GitHub:

```bash
pip install git+https://github.com/yyxwjq/qctools.git
```

Optional NEP support requires `pynep`:

```bash
pip install -e .[nep]
```

## Python Tools

### Logging

```python
import qctools

qctools.qctools_logging(filename="qctools.log", overwrite=False)
```

`overwrite=True` removes an existing log file before configuring logging.

### RDF

```python
from ase.io import read
from qctools.rdf import get_rdf

images = read("trajectory.xyz", ":")
get_rdf(images, cutoff=5.0, bin_size=0.1, first_neighbor=False, cores=4)
```

RDF output is written to `RDF/` or `RDF_first/`. The current implementation
uses a frame-averaged, density-normalized `g(r)`.

### ADF

```python
from ase.io import read
from qctools.adf import get_adf

images = read("trajectory.xyz", ":")
get_adf(images, rcut=4.0, bin_size=5.0, cores=4)
```

ADF output is written to `ADF/`. The cutoff is applied to the actual pair
distance from the central atom to each neighbor.

### Coordination Number

```python
from qctools.coord import group_coordnum

coord_nums = group_coordnum(
    traj=images,
    group1=[0, 1, 2],
    group2=[3, 4, 5],
    r0=4.0,
    cores=4,
)
```

### ML Potential Error Analysis

```python
from qctools.ml.error_img import main

main(
    trajname="trajectory.xyz",
    apps="nep",          # "nep" or "n2p2"
    resource="software", # "software" or "images"
    fontsize=12,
    data={"energy": "energy.txt", "force": "force.txt"},
    er_bar=1.5,
    show_marginals=True,
)
```

Main generated files include:

- `energy.data`, `force.data`
- `energy_error.txt`, `force_error.txt`
- `Err-energy.xyz`
- `Err-force-ini.xyz`, `Err-force-replaced.xyz`
- `leave-E-img.xyz`, `leave-F-img.xyz`
- `energy_error_analysis.png`, `force_error_analysis.png`

### Structure Editing

```python
from qctools.edit_atoms import remove

remove("structure.vasp", ["H"])
```

Command-line entry point:

```bash
qctools-edit remove structure.vasp H
```

## VASP Job Monitor

`qctools/vaspjob` is a Bash helper for scanning VASP calculation directories.
It searches for `INCAR` files, checks related `POSCAR`, `OSZICAR`, and `OUTCAR`
files, and writes a summary table to `vaspjob.log`.

Run from the directory that contains your VASP tasks:

```bash
qctools/vaspjob
qctools/vaspjob -v
```

Options:

- `-v`: also print the report to the terminal.
- `-m`: generate `wavecar_list.txt` with all `WAVECAR` files and sizes.
- `-r`: remove `WAVECAR` files listed in `wavecar_list.txt`.

`-r` is destructive. Review `wavecar_list.txt` before running it.

The monitor recognizes common task types such as single-point, ion
optimization, cell optimization, frequency, DOS, dimer, COHP marked by
`lobsterin`/`LOBSTERIN`, and NEB directories.

## Dependencies

- Python >= 3.7
- ASE
- NumPy
- SciPy
- Matplotlib
- PyNEP, optional for NEP calculator workflows

## Tests

The regression tests use the standard library `unittest`:

```bash
python -m unittest discover -s tests -v
python -m compileall -q qctools tests
```

## License

MIT License.
