# -*- coding: utf-8 -*-

"""
Author: Jiaqi Wang
Date: 2025-07-15 21:03
This script calculates coordination number of atoms with resource-efficient parallel processing using starmap.
"""

from ase.geometry import get_distances
from qctools import qctools_logging
from multiprocessing import Pool
from ase.io import read
import numpy as np
import os
import logging

qctools_logging()

def _group_coordnum_serial(image, group1, group2, r0=4.0, en=6.0, ed=12.0, tolerance=0.0):
    """
    Calculate the coordination number between two groups of atoms using ASE's get_distances, vectorized.

    Parameters:
    image : ASE Atoms object
        The atomic configuration.
    group1 : list
        Atomic indices of atoms in group 1.
    group2 : list
        Atomic indices of atoms in group 2.
    r0 : float
        Cutoff distance (in Angstroms).
    en : int
        Numerator exponent (must be even).
    ed : int
        Denominator exponent (must be even).
    tolerance : float
        Pairlist tolerance (default 0.0).

    Returns:
    float
        The coordination number.
    """
    try:
        if en % 2 != 0 or ed % 2 != 0:
            raise ValueError("Exponents en and ed must be even integers.")
        if en / ed != 1/2:
            raise ValueError("Numerator / Denominator exponent must equal 1/2.")
        if en <= 0 or ed <= 0:
            raise ValueError("Exponents en and ed must be positive.")

        group1_pos = image[group1].positions
        group2_pos = image[group2].positions
        D, D_len = get_distances(group1_pos, group2_pos, cell=image.cell, pbc=image.pbc)

        en2 = en // 2
        ed2 = ed // 2
        l2 = (D_len / r0) ** 2
        xn = l2 ** en2
        xd = l2 ** ed2
        func = ((1.0 - xn) / (1.0 - xd) - tolerance) / (1.0 - tolerance)
        return np.sum(func[(func > 0) & (func < 1)])
    except Exception as e:
        logging.error(f"Error processing frame: {e}")
        raise


def _process_chunk(start_index, atoms_list, group1, group2, r0, en, ed, tolerance):
    """
    Process a chunk of trajectory frames, returning results with their original indices.

    Parameters:
    start_index : int
        Starting index of the chunk in the original trajectory.
    atoms_list : list
        List of ASE Atoms objects for the chunk.
    group1 : list
        Atomic indices of atoms in group 1.
    group2 : list
        Atomic indices of atoms in group 2.
    r0 : float
        Cutoff distance (in Angstroms).
    en : float
        Numerator exponent (must be even).
    ed : float
        Denominator exponent (must be even).
    tolerance : float
        Pairlist tolerance (default 0.0).

    Returns:
    list
        List of (index, coord_num) tuples for the chunk.
    """
    try:
        results = []
        for i, atoms in enumerate(atoms_list):
            coord_num = _group_coordnum_serial(atoms, group1, group2, r0, en, ed, tolerance)
            results.append((start_index + i, coord_num))
        return results
    except Exception as e:
        logging.error(f"Error processing chunk starting at index {start_index}: {e}")
        raise

def group_coordnum(traj, group1, group2, r0=4.0, en=6, ed=12, tolerance=0.0, cores=None):
    """
    Compute coordination numbers for an MD trajectory in parallel, preserving frame order.

    Parameters:
    traj : list[Atoms] or str
        ASE Atoms list or path to trajectory file (ASE-compatible format, e.g., XYZ, NetCDF).
    group1 : list
        Atomic indices of atoms in group 1.
    group2 : list
        Atomic indices of atoms in group 2.
    r0 : float
        Cutoff distance (in Angstroms).
    en : int
        Numerator exponent (must be even).
    ed : int
        Denominator exponent (must be even).
    tolerance : float
        Pairlist tolerance (default 0.0).
    cores : int or None
        Number of processes for parallel computation. If None, uses half CPU count.

    Returns:
    np.ndarray
        Array of coordination numbers in trajectory frame order.
    """
    try:
        # Handle trajectory input
        if isinstance(traj, str):
            logging.info(f"Loading trajectory from {traj}")
            traj = read(traj, index=':')
        elif not isinstance(traj, (list, tuple)) or not all(hasattr(atoms, 'positions') for atoms in traj):
            raise ValueError("traj must be a list of ASE Atoms objects or a valid trajectory file path.")

        n_frames = len(traj)
        logging.info(f"Processing {n_frames} frames with {len(group1)} and {len(group2)} atoms in groups")

        # Validate indices
        max_index = max(max(group1, default=-1), max(group2, default=-1))
        if max_index >= len(traj[0]):
            raise ValueError(f"Group indices exceed number of atoms ({len(traj[0])}) in trajectory")

        # Set number of processes
        if cores is None:
            cores = max(1, os.cpu_count() // 2)  # Use half the cores to avoid overloading
        cores = min(cores, n_frames)  # Don’t use more processes than frames
        logging.info(f"Using {cores} core(s) for parallel processing")

        # Evenly distribute frames across cores
        chunk_sizes = [n_frames // cores + (1 if i < n_frames % cores else 0) for i in range(cores)]
        chunks = []
        current_index = 0
        for i, size in enumerate(chunk_sizes):
            if size > 0:  # Only include non-empty chunks
                chunk = (current_index, traj[current_index:current_index + size], group1, group2, r0, en, ed, tolerance)
                chunks.append(chunk)
                logging.info(f"Core {i} will process {size} frames (indices {current_index} to {current_index + size - 1})")
                current_index += size

        # Parallel processing with starmap
        try:
            with Pool(processes=cores) as pool:
                chunk_results = pool.starmap(_process_chunk, chunks)
            logging.info("Parallel processing completed successfully")
        except Exception as e:
            logging.warning(f"Parallel processing failed: {e}. Falling back to serial processing")
            chunk_results = [_process_chunk(*chunk_data) for chunk_data in chunks]

        # Flatten and sort results by original index
        results = []
        for chunk_result in chunk_results:
            results.extend(chunk_result)
        results.sort(key=lambda x: x[0])  # Sort by original frame index
        coord_nums = np.array([result[1] for result in results])

        logging.info(f"Completed processing. Output array size: {len(coord_nums)}")
        return coord_nums

    except Exception as e:
        logging.error(f"Error in group_coordnum: {e}")
        raise