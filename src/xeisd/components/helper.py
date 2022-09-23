"""
Contains functions to help X-EISD perform calculations with
third-party programs.

Initial creation made with integrated SPyCi-PDB calculators in mind.
"""
from functools import partial

from spycipdb.components.helpers import *

from idpconfgen.libs.libmulticore import pool_function

def selective_calculator(pdbfilepaths, modules, ncores):
    """
    Back-calculate experimental data with no pre-existing back-calculations.
    
    Output of SPyCi-PDB is matched to that of `Stack` objects so no additional
    "parsing" needed.
    
    Parameters
    ----------
    pdbfilepaths : list
        List of paths that point to PDBs to perform BC on
    
    modules : list
        List of experimental modules to perform BC on
    
    ncores : int
        Number of workers to initialize multiprocessing for    
    """