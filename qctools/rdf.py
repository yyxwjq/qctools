# -*- coding: utf-8 -*-

"""
Author: Jiaqi Wang
Date: 2025-07-09 23:45
This script calculates the radial distribution function (RDF) for a given set of atomic structures.
It uses ASE for atomic structure manipulation and analysis, and can handle multiple structures in parallel.
"""
from qctools.element_tools import get_elements, elements_iterators
from ase.neighborlist import neighbor_list, first_neighbors
import matplotlib.pyplot as plt
import multiprocessing as mp
import numpy as np
import logging
import time
import os


logger = logging.getLogger(__name__)

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
    tasks = [(chunk, elements, cutoff, bin_size, first_neighbor) for chunk in chunks if chunk]
    if cores <= 1:
        results = [get_radial_distribution(*task) for task in tasks]
    else:
        with mp.Pool(processes=cores) as pool:
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
    _, metadata = get_radial_distribution(
        images, elements, cutoff, bin_size, first_neighbor, collect_metadata=True
    )
    normalized = normalize_rdf(combined, metadata, cutoff, bin_size)
    for elements_comb, r_gr in normalized.items():
        bins = r_gr[:, 0]
        grfunc = r_gr[:, 1]
        e1, e2 = elements_comb
        if first_neighbor:
            filename = f"rdf_first_{e1}_{e2}.dat"
        else:
            filename = f"rdf_{e1}_{e2}.dat"
        np.savetxt(os.path.join(dirname, filename), r_gr, fmt='%.4f\t%.2f',header=f"RDF for {e1} - {e2} with cutoff {cutoff} and bin size {bin_size}")
        plt.plot(bins, grfunc, label=f"{e1} - {e2}")
        first_y = max(grfunc) if len(grfunc) else 0
        first_x = bins[grfunc.tolist().index(first_y)] if len(grfunc) else 0
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

def get_radial_distribution(images, elements, cutoff, bin_size, first_neighbor=False, collect_metadata=False):
    """Calculate the radial distribution function (RDF) for a single image
    :param images: List of ASE images
    :param elements: List of elements to consider for RDF calculation
    :param cutoff: Cutoff distance for neighbor search
    :param bin_size: Size of the bins for RDF histogram
    :param first_neighbor: If True, only consider the first neighbors for RDF calculation
    :return: Dictionary with element pairs as keys and RDF histograms as values
    """
    dist_dict = None
    metadata = {
        "n_frames": 0,
        "volumes": [],
        "element_counts": {},
    }
    for image in images:
        # Initial parameters
        n_bins = int(np.ceil(cutoff / bin_size))
        counter_dict = {}
        natoms = len(image)
        metadata["n_frames"] += 1
        metadata["volumes"].append(float(image.get_volume()) if image.cell.rank == 3 else np.nan)
        symbols = image.get_chemical_symbols()
        for element in elements:
            metadata["element_counts"][element] = metadata["element_counts"].get(element, 0) + symbols.count(element)
        
        # Obtain all element combinations
        combinations = elements_iterators(elements, 'rdf')
        if not combinations:
            raise ValueError("No valid element combinations found for RDF calculation")
        
        # Obtain all neighbors
        for comb in combinations: 
            counter_dict[comb] = np.zeros(n_bins)                      
            # Initialize neighbor list
            e1, e2 = comb
            rcut = {comb : cutoff}
            idxi, idxj, distances = neighbor_list('ijd', image, cutoff=rcut)
            pair_mask = np.array(
                [symbols[i] == e1 and symbols[j] == e2 for i, j in zip(idxi, idxj)],
                dtype=bool,
            )
            idxi = idxi[pair_mask]
            idxj = idxj[pair_mask]
            distances = distances[pair_mask]
            # Just get first neighbors for target elements rdf analysis if needed
            if first_neighbor:
                # Check if we have the expected number of first neighbors
                indices = image.symbols.indices()
                target_natoms = len(indices[e1])
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
                    if first_neighbor_k < cutoff:
                        bin_index = min(int(first_neighbor_k / bin_size), n_bins - 1)
                        counter_dict[comb][bin_index] += 1
            else:
                for d in distances:
                    if d < cutoff:
                        bin_index = min(int(d / bin_size), n_bins - 1)
                        counter_dict[comb][bin_index] += 1
            
        if dist_dict is None:
            dist_dict = {key: arr.copy() for key, arr in counter_dict.items()}
        else:
            # Combine results from multiple images
            for key, arr in counter_dict.items():
                dist_dict[key] += arr
    if collect_metadata:
        return dist_dict, metadata
    return dist_dict


def normalize_rdf(counts, metadata, cutoff, bin_size):
    """Normalize RDF histograms to standard dimensionless g(r)."""
    n_bins = int(np.ceil(cutoff / bin_size))
    radii = (np.arange(n_bins) + 0.5) * bin_size
    shell_volumes = 4.0 * np.pi * radii * radii * bin_size
    n_frames = metadata.get("n_frames", 0)
    volumes = np.array(metadata.get("volumes", []), dtype=float)
    finite_volumes = volumes[np.isfinite(volumes) & (volumes > 0)]
    avg_volume = float(np.mean(finite_volumes)) if len(finite_volumes) else np.nan
    element_counts = metadata.get("element_counts", {})
    normalized = {}

    for comb, values in counts.items():
        e1, e2 = comb
        avg_e1 = element_counts.get(e1, 0) / n_frames if n_frames else 0
        avg_e2 = element_counts.get(e2, 0) / n_frames if n_frames else 0
        if not np.isfinite(avg_volume) or avg_e1 == 0 or avg_e2 == 0:
            grfunc = np.zeros_like(values, dtype=float)
        else:
            number_density = avg_e2 / avg_volume
            denominator = n_frames * avg_e1 * number_density * shell_volumes
            grfunc = np.divide(
                values,
                denominator,
                out=np.zeros_like(values, dtype=float),
                where=denominator > 0,
            )
        normalized[comb] = np.vstack((radii, grfunc)).T
    return normalized
