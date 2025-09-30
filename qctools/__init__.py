# -*- coding: utf-8 -*-
"""
QCTools - Quantum Chemistry Analysis Toolkit

A comprehensive Python toolkit for quantum chemistry calculations and analysis,
designed for efficient processing of atomic structures and molecular dynamics trajectories.

Author: Jiaqi Wang
Version: 1.0.0
"""

import logging
import os

# Version information
__version__ = "1.0.0"
__author__ = "Jiaqi Wang"
__email__ = "wangjiaqi@example.com"
__description__ = "A comprehensive Python toolkit for quantum chemistry calculations and analysis"

# Import main modules
try:
    from . import rdf
    from . import adf
    from . import coord
    from . import element_tools
    from . import edit_atoms
    from . import radius
    from . import ml
except ImportError as e:
    import warnings
    warnings.warn(f"Warning: Failed to import some modules: {e}", ImportWarning)

def qctools_logging(level=logging.INFO, filename='qctools.log'):
    """
    Initialize logging system for QCTools.
    
    Parameters:
    -----------
    level : int, optional
        Logging level (default: logging.INFO)
    filename : str, optional
        Log file name (default: 'qctools.log')
    """
    if os.path.exists('qctools.log'):
        os.remove('qctools.log')
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        filename=filename,
        filemode='a'
    )
    
_default_logger = logging.getLogger(__name__)
_default_logger.addHandler(logging.NullHandler())  # Forbid "No handlers"

# All available modules and functions
__all__ = [
    'qctools_logging',
    'rdf',
    'adf', 
    'coord',
    'element_tools',
    'edit_atoms',
    'radius',
    'ml',
    '__version__',
    '__author__',
    '__email__',
    '__description__'
]
