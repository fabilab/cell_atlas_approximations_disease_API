"""
Setup script for atlasapprox_disease Python package
"""
import pathlib
from setuptools import setup, find_packages


# version
version_file = pathlib.Path(__file__).parent / "VERSION"
with open(version_file) as f:
    version = f.read().rstrip('\n')


long_description = """[![Documentation Status](https://readthedocs.org/projects/atlasapprox-disease/badge/?version=latest)](https://apidocs.atlasapprox-disease.org/en/latest/?badge=latest)

Python interface to disease-specific cell atlas approximations
==============================================================
This project provides access to disease-specific cell atlas approximations, which are lightweight and lossy compressions of cell atlas data focused on disease contexts. 

Built using the cellxgene census dataset as an initial source of data and the scquill package for data compression, the package enables biologists, doctors, and data scientists to quickly answer key questions such as:

- *What is the differential cell abundance of each cell type in lung tissue between normal and disease states?*
- *What are the top differentially expressed genes in diseases affecting the lung?*
- *How does cell type abundance change across organs for a given disease?*

In addition to this Python package, these questions can be addressed in R or through the REST API.

**Development**: https://github.com/fabilab/cell_atlas_approximations_disease_API
"""


setup(
    name="atlasapprox_disease",
    url="https://apidocs.atlasapprox-disease.org",
    description="Disease-specific cell atlas approximations, Python API",
    long_description=long_description,
    long_description_content_type='text/markdown',
    license="MIT",
    author="Ying Xu",
    author_email="ying.xu3@unsw.edu.au",
    version=version,
    packages=find_packages(),
    install_requires=[
        "requests",
        "numpy",
        "pandas",
    ],
    python_requires=">=3.8",
    platforms="ALL",
    keywords=[
        "single cell",
        "cell atlas",
        "omics",
        "biology",
        "disease",
        "bioinformatics",
    ],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3 :: Only",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
)
