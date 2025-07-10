# -*- coding: utf-8 -*-
"""
Author: Jiaqi Wang
Date: 2025-07-09 23:45
This script calculates the angular distribution function (ADF) for a given set of atomic structures.
It uses ASE for atomic structure manipulation and analysis, and can handle multiple structures in parallel.
"""
from qctools.element_tools import get_elements, elements_iterators
from ase.geometry.analysis import Analysis
from ase.neighborlist import NeighborList
from qctools import qctools_logging
import matplotlib.pyplot as plt
import multiprocessing as mp
from ase.io import read
import numpy as np
import logging
import time
import os

qctools_logging()

def get_adf(images, rcut, bin_size, cores=None):
    """Process a list of images using multiprocessing
    :param images: List of ASE images to process
    :param rcut: Cutoff distance for ADF calculation
    :param bin_size: Size of the bins for ADF histogram
    :param cores: Number of CPU cores to use for parallel processing, defaults to None which uses all available cores
    :return: Present the ADF results of each element combination in the form of data file and graphs
    """  
    
    """
    Obtain radial distribution for images in parallel
    """
    
    combined = {}
    elements = get_elements(images)
    if cores is None:
        cores = 1
    elif cores > mp.cpu_count():
        logging.warning(f"Requested {cores} cores, but only {mp.cpu_count()} are available. Using all available cores.")
        cores = mp.cpu_count()
    # Allocate images to cores average evenly
    chunk_sizes = [len(images) // cores + (1 if i < len(images) % cores else 0) for i in range(cores)]
    chunks = []
    start = 0
    for size in chunk_sizes:
        end = min(start + size, len(images))
        chunks.append(images[start:end])
        start = end
    time_start = time.time()
    with mp.Pool(processes=cores) as pool:
        tasks = [(chunk, elements, rcut, bin_size) for chunk in chunks]

        # Use tqdm to display progress bar
        results = pool.starmap(get_angular_distribution, tasks)
        for result in results:
            for key, value in result.items():
                if key not in combined:
                    combined[key] = value
                else:
                    combined[key] += value
    """
    Obtain ADF  for each element combination
    """
    dirname = './ADF'
    if not os.path.exists(dirname):
        os.makedirs(dirname, exist_ok=True)
    for elements_comb, bins_values in combined.items():
        n_bins = len(bins_values)
        # Normalize the histogram
        bins = np.arange(bin_size, (n_bins + 1) * bin_size, bin_size)
        angular_func = bins_values
        a_afunc = np.vstack((bins, angular_func)).T
        e1, e2, e3 = elements_comb
        filename = f"adf_{e1}_{e2}_{e3}.dat"
        np.savetxt(os.path.join(dirname, filename), a_afunc, fmt='%.4f\t%.2f',header=f"ADF for {e1}-{e2}-{e3} with rcut {rcut} and bin size {bin_size}°.")
        plt.plot(bins, angular_func, label=f"{e1} - {e2} - {e3}")
        first_y = max(angular_func)
        first_x = bins[angular_func.tolist().index(first_y)]
        plt.xlabel('Angle (°)')
        plt.ylabel('Frequency')
        plt.title('Angular Distribution')
        plt.text(first_x, first_y, s="Max peak: {:2f}".format(first_x), fontsize=12)
        plt.legend()
        plt.grid()
        fig_name = f"adf_{e1}_{e2}_{e3}.png"
        plt.savefig(os.path.join(dirname, fig_name))
        plt.clf()
    time_end = time.time()
    logging.info(f"ADF running time: {time_end - time_start:.2f}s.")

def get_angular_distribution(images, elements, rcut, bin_size):
    """Calculate the radial distribution function (ADF) for a single image
    :param images: List of ASE images
    :param elements: List of elements to consider for ADF calculation
    :param rcut: Cutoff distance for neighbor search
    :param bin_size: Size of the bins for ADF histogram
    :return: Dictionary with element pairs as keys and ADF histograms as values
    """
    ang_dict = None
    for image in images:
        # Initial parameters
        min_b = 1  # Minimum bin value to avoid division by zero
        
        n_bins = int(180 / bin_size)
        counter_dict = {}
        natoms = len(image)
        
        # Obtain all element combinations
        combinations = elements_iterators(elements, 'adf')
        if not combinations:
            raise ValueError("No valid element combinations found for ADF calculation")

        # Obtain all neighbors
        for comb in combinations: 
            counter_dict[comb] = np.zeros(n_bins)
            # Initialize neighbor list
            cutoff = [rcut for _ in range(len(image))]
            nl = NeighborList(cutoff, skin=0)
            nl.update(image)
            ac = Analysis(image, nl)
            ea, ecenter, eb = comb
            # Just get first neighbors for target elements rdf analysis if needed
            if len(ac.get_angles(ea, ecenter, eb)[0]) == 0:
                logging.info(f"No statisfied angles. Skip {ea}-{ecenter}-{eb} angles.")
                continue
            angles = ac.get_values(ac.get_angles(ea, ecenter, eb, unique=True))[0]
            for a in angles:
                counter_dict[comb][int(a / bin_size)] += 1
            
        if ang_dict is None:
            ang_dict = {key: arr.copy() for key, arr in counter_dict.items()}
        else:
            # Combine results from multiple images
            for key, arr in counter_dict.items():
                ang_dict[key] += arr
    return ang_dict
