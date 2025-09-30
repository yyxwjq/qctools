#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
QCTools - Quantum Chemistry Analysis Toolkit
Setup configuration for pip installation
"""

from setuptools import setup, find_packages
import os

# Read the contents of README file
this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

# Read requirements from requirements.txt
def read_requirements():
    """Read requirements from requirements.txt file"""
    req_path = os.path.join(this_directory, 'requirements.txt')
    if os.path.exists(req_path):
        with open(req_path, 'r', encoding='utf-8') as f:
            return [line.strip() for line in f 
                   if line.strip() and not line.startswith('#')]
    return []

setup(
    name="qctools",
    version="1.0.0",
    author="Jiaqi Wang",
    author_email="wangjiaqi@example.com",
    description="A comprehensive Python toolkit for quantum chemistry calculations and analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/qctools",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
    install_requires=read_requirements(),
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov",
            "black",
            "flake8",
            "mypy",
        ],
        "nep": [
            "pynep",
        ],
        "all": [
            "pynep",
        ],
    },
    entry_points={
        "console_scripts": [
            "qctools-edit=qctools.edit_atoms:main",
        ],
    },
    include_package_data=True,
    package_data={
        "qctools": ["*.md", "*.txt"],
    },
    keywords=[
        "quantum chemistry",
        "molecular dynamics",
        "computational chemistry",
        "structure analysis",
        "machine learning potential",
        "RDF",
        "ADF",
        "coordination number",
        "ASE",
    ],
    project_urls={
        "Bug Reports": "https://github.com/yourusername/qctools/issues",
        "Source": "https://github.com/yourusername/qctools",
        "Documentation": "https://github.com/yourusername/qctools/wiki",
    },
)