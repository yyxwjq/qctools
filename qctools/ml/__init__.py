# -*- coding: utf-8 -*-
"""
QCTools ML Module - Machine Learning Potential Analysis

This module provides tools for analyzing machine learning potentials,
including error analysis and validation for NEP and n2p2 potentials.
"""

try:
    from . import error_img
    __all__ = ['error_img']
except ImportError as e:
    import warnings
    warnings.warn(f"Warning: Failed to import ML modules: {e}", ImportWarning)
    __all__ = []