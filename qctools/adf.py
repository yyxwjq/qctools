# -*- coding: utf-8 -*-
"""
Author: Jiaqi Wang
Date: 2025-07-09 23:45
This script calculates the angular distribution function (ADF) for a given set of atomic structures.
It uses ASE for atomic structure manipulation and analysis, and can handle multiple structures in parallel.
"""
from qctools.element_tools import get_elements, elements_iterators
from ase.geometry import find_mic
import matplotlib.pyplot as plt
import multiprocessing as mp
import numpy as np
import logging
import time
import os

logger = logging.getLogger(__name__)

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
    tasks = [(chunk, elements, rcut, bin_size) for chunk in chunks if chunk]
    if cores <= 1:
        results = [get_angular_distribution(*task) for task in tasks]
    else:
        with mp.Pool(processes=cores) as pool:
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
    n_bins = int(np.ceil(180.0 / bin_size))
    for image in images:
        counter_dict = {}
        
        # Obtain all element combinations
        combinations = elements_iterators(elements, 'adf')
        if not combinations:
            raise ValueError("No valid element combinations found for ADF calculation")

        symbols = image.get_chemical_symbols()
        positions = image.get_positions()

        for comb in combinations: 
            counter_dict[comb] = np.zeros(n_bins)
            ea, ecenter, eb = comb

            for center_idx, center_symbol in enumerate(symbols):
                if center_symbol != ecenter:
                    continue

                neighbors_a = []
                neighbors_b = []
                center = positions[center_idx]
                for neighbor_idx, neighbor_symbol in enumerate(symbols):
                    if neighbor_idx == center_idx:
                        continue
                    vector = positions[neighbor_idx] - center
                    vector, distance = find_mic(vector, image.cell, image.pbc)
                    if distance > rcut:
                        continue
                    if neighbor_symbol == ea:
                        neighbors_a.append((neighbor_idx, vector))
                    if neighbor_symbol == eb:
                        neighbors_b.append((neighbor_idx, vector))

                for idx_a, vector_a in neighbors_a:
                    for idx_b, vector_b in neighbors_b:
                        if idx_a == idx_b:
                            continue
                        if ea == eb and idx_a >= idx_b:
                            continue
                        norm = np.linalg.norm(vector_a) * np.linalg.norm(vector_b)
                        if norm == 0:
                            continue
                        cos_angle = np.dot(vector_a, vector_b) / norm
                        angle = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
                        bin_index = min(int(angle / bin_size), n_bins - 1)
                        counter_dict[comb][bin_index] += 1

            if not counter_dict[comb].any():
                logging.info(f"No statisfied angles. Skip {ea}-{ecenter}-{eb} angles.")
            
        if ang_dict is None:
            ang_dict = {key: arr.copy() for key, arr in counter_dict.items()}
        else:
            # Combine results from multiple images
            for key, arr in counter_dict.items():
                ang_dict[key] += arr
    return ang_dict
