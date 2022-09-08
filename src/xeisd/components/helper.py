"""
Contains functions to help X-EISD perform calculations with
third-party programs.

Initial creation made with integrated SPyCi-PDB calculators in mind.
"""
from spycipdb.components.calculators import *

def selective_calculator(pdbfilepaths, modules):
    """
    Back-calculate experimental data with no pre-existing back-calculations.
    
    """
    