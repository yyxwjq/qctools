# ML Error Analysis - Data Parameter Usage

## Overview
The `data` parameter in the ML error analysis tool now supports both string and dictionary formats for maximum flexibility.

## Usage Examples

### 1. Dictionary Format (Recommended for software resource)
```python
from qctools.ml.error_img import main

# Analyze both energy and force simultaneously
main(trajname='train.xyz',
     apps='nep',
     resource='software',
     fontsize=12,
     data={'energy': 'energy_results.txt', 'force': 'force_results.txt'},
     er_bar=2.5,
     ra={'Pt':'Pd', 'H':'He', 'O':'F'},
     cut_img=True,
     comment=False,
     show_marginals=True
)

# Analyze only energy
main(trajname='train.xyz', apps='nep', resource='software',
     data={'energy': 'energy_results.txt'}, ...)

# Analyze only force  
main(trajname='train.xyz', apps='nep', resource='software',
     data={'force': 'force_results.txt'}, ...)
```

### 2. Legacy String Format (Still Supported)
```python
# Auto-detect energy or force based on file format
main(trajname='train.xyz', apps='nep', resource='software',
     data='results.txt', ...)
```

### 3. Images Resource (Unchanged)
```python
# For images resource, data parameter is not used
main(trajname='train.xyz', apps='nep', resource='images',
     pot='nep.txt', ...)
```

## Data Parameter Options

| Format | Usage | Description |
|--------|-------|-------------|
| `str` | `data='file.txt'` | Legacy format, auto-detect energy/force |
| `dict` | `data={'energy': 'e.txt', 'force': 'f.txt'}` | Specify both files |
| `dict` | `data={'energy': 'e.txt'}` | Energy only |
| `dict` | `data={'force': 'f.txt'}` | Force only |
| `None` | `data=None` | For images resource |

## File Format Requirements

### NEP Format
- **Energy files**: 2 columns (DFT_energy, NEP_energy)
- **Force files**: 6 columns (DFT_fx, DFT_fy, DFT_fz, NEP_fx, NEP_fy, NEP_fz)

### N2P2 Format  
- **Energy files**: 3 columns (structure_id, DFT_energy, NEP_energy)
- **Force files**: 4 columns (structure_id, atom_id, component, DFT_force, NEP_force)

## Benefits
- ✅ Simultaneous energy and force analysis for software resource
- ✅ Flexible data input options
- ✅ Backward compatibility with existing code
- ✅ Clear error messages for format issues