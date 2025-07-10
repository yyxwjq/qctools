# -*- coding: utf-8 -*-

"""
Author: Jiaqi Wang
Date: 2025-07-09 23:45
This script calculates the radial distribution function (RDF) for a given set of atomic structures.
It uses ASE for atomic structure manipulation and analysis, and can handle multiple structures in parallel.
"""
from qctools.element_tools import get_elements, elements_iterators
from ase.neighborlist import neighbor_list, first_neighbors
from qctools import qctools_logging
import matplotlib.pyplot as plt
import multiprocessing as mp
from ase.io import read
import numpy as np
import logging
import time
import os


qctools_logging()

def get_rdf(images, cutoff, bin_size, first_neighbor=False, cores=None):
    """Process a list of images using multiprocessing
    :param images: List of ASE images to process
    :param cutoff: Cutoff distance for RDF calculation
    :param bin_size: Size of the bins for RDF histogram
    :param first_neighbor: If True, only consider first neighbors for RDF calculation
    :param cores: Number of CPU cores to use for parallel processing, defaults to None which uses all available cores
    :return: Dictionary containing RDF results for each element combination
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
        tasks = [(chunk, elements, cutoff, bin_size, first_neighbor) for chunk in chunks]

        # Use tqdm to display progress bar
        results = pool.starmap(get_radial_distribution, tasks)
        for result in results:
            for key, value in result.items():
                if key not in combined:
                    combined[key] = value
                else:
                    combined[key] += value
    """
    Obtain g(r) functions for each element combination
    """
    if first_neighbor:
        dirname = "./RDF_first"
    else:
        dirname = "./RDF"
    if not os.path.exists(dirname):
        os.makedirs(dirname, exist_ok=True)
    for elements_comb, bins_values in combined.items():
        n_bins = len(bins_values)      
        # Normalize the histogram
        bins = np.arange(bin_size, (n_bins + 1) * bin_size, bin_size)
        grfunc = bins_values / (4 * np.pi * bins * bins * bin_size)
        r_gr = np.vstack((bins, grfunc)).T
        e1, e2 = elements_comb
        if first_neighbor:
            filename = f"rdf_first_{e1}_{e2}.dat"
        else:
            filename = f"rdf_{e1}_{e2}.dat"
        np.savetxt(os.path.join(dirname, filename), r_gr, fmt='%.4f\t%.2f',header=f"RDF for {e1} - {e2} with cutoff {cutoff} and bin size {bin_size}")
        plt.plot(bins, grfunc, label=f"{e1} - {e2}")
        first_y = max(grfunc)
        first_x = bins[grfunc.tolist().index(first_y)]
        plt.xlabel('Distance (Å)')
        plt.ylabel('g(r)')
        plt.title('Radial Distribution Function')
        plt.text(first_x, first_y, s="Max peak: {:2f}".format(first_x), fontsize=12)
        plt.legend()
        plt.grid()
        if first_neighbor:
            fig_name = f"rdf_first_{e1}_{e2}.png"
        else:
            fig_name = f"rdf_{e1}_{e2}.png"
        plt.savefig(os.path.join(dirname, fig_name))
        plt.clf()
    time_end = time.time()
    logging.info(f"RDF running time: {time_end - time_start:.2f}s.")

def get_radial_distribution(images, elements, cutoff, bin_size, first_neighbor=False):
    """Calculate the radial distribution function (RDF) for a single image
    :param images: List of ASE images
    :param elements: List of elements to consider for RDF calculation
    :param cutoff: Cutoff distance for neighbor search
    :param bin_size: Size of the bins for RDF histogram
    :param first_neighbor: If True, only consider the first neighbors for RDF calculation
    :return: Dictionary with element pairs as keys and RDF histograms as values
    """
    dist_dict = None
    for image in images:
        # Initial parameters
        min_b = 0.0000000001  # Minimum bin value to avoid division by zero
        
        n_bins = int(cutoff / bin_size)
        counter_dict = {}
        natoms = len(image)
        
        # Obtain all element combinations
        combinations = elements_iterators(elements, 'rdf')
        if not combinations:
            raise ValueError("No valid element combinations found for RDF calculation")
        
        # Obtain all neighbors
        for comb in combinations: 
            counter_dict[comb] = np.zeros(n_bins)                      
            # Initialize neighbor list
            rcut = {comb : cutoff}
            idxi, idxj, distances = neighbor_list('ijd', image, cutoff=rcut)
            # Just get first neighbors for target elements rdf analysis if needed
            if first_neighbor:
                # Check if we have the expected number of first neighbors
                e1, e2 = comb
                indices = image.symbols.indices()
                if e1 == e2:
                    target_natoms = len(indices[e1])
                else:
                    target_natoms = len(indices[e1]) + len(indices[e2])
                actual_natoms = len(np.unique(idxi))
                if actual_natoms != target_natoms:
                    logging.warning(f"Expected return {target_natoms} neighbors for {comb} system, "
                                    f"but only found {actual_natoms}, please check cutoff value.")
                # Obtain first neighbors for target element atoms.
                neighbors_matrix = first_neighbors(natoms=natoms, first_atom=idxi)
                for k in np.unique(idxi):
                    start, end = neighbors_matrix[k], neighbors_matrix[k + 1]
                    distances_k = distances[start:end]
                    first_neighbor_k = distances_k[np.argmin(distances_k)]
                    counter_dict[comb][int(first_neighbor_k / bin_size)] += 1
            else:
                for d in distances:
                    counter_dict[comb][int(d / bin_size)] += 1
            
        if dist_dict is None:
            dist_dict = {key: arr.copy() for key, arr in counter_dict.items()}
        else:
            # Combine results from multiple images
            for key, arr in counter_dict.items():
                dist_dict[key] += arr
    return dist_dict


