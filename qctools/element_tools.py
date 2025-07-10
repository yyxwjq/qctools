# -coding: utf-8 -
from ase import Atoms
import logging
"""
Author: Jiaqi Wang
Date: 2025-07-09 23:58
This script provides utility functions to deal with elements.
"""
def get_elements(images):
    """Obtain unique elements from a list of images
    :param images: List of ASE images
    :return: Sorted list of unique elements found in the images
    """
    if not images:
        logging.warning("No images found")
        return []
    
    elements = set()
    for i in range(len(images)):
        try:
            # Detect if image is empty
            if len(images[i]) == 0:
                logging.warning(f"Found empty image at index {i}, skipping")
                continue

            # Use the correct attribute to get elements
            elements = elements | set(images[i].symbols.species())
        except AttributeError as e:
            logging.error(f"Error getting elements: {e}")
    return sorted(elements)

def elements_iterators(elements, func_type='rdf'):
    """Obtain combinations of elements for RDF or ADF calculations
    :param elements: List of elements to combine
    :param func_type: Type of function to generate combinations, either 'rdf' for radial distribution function or 'adf' for angular distribution function
    :return: Set of tuples containing combinations of elements
    """
    if not elements:
        logging.warning("No elements provided")
        return set()
    
    s_elements = sorted(elements)
    combinations = set()
    
    if func_type == 'rdf':
        for i, ele1 in enumerate(s_elements):
            for ele2 in s_elements[i:]:
                combinations.add((ele1, ele2))
    elif func_type == 'adf':
        for center_element in s_elements:
            for i, ele1 in enumerate(s_elements):
                for ele2 in s_elements[i:]:
                    # Include symmetric combinations with the central element
                    combinations.add((ele1, center_element, ele2))
    else:
        raise NotImplementedError(f"Function type '{func_type}' is not supported")
    return combinations