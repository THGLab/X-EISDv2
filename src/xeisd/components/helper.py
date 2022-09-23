"""
Contains functions to help X-EISD perform calculations with
third-party programs.

Initial creation made with integrated SPyCi-PDB calculators in mind.
"""
from functools import partial

from xeisd import log
from xeisd.logger import S, T, init_files, report_on_crash

from xeisd.components.parser import Stack
from xeisd.components import (
    default_bc_errors,
    saxs_name,
    cs_name,
    fret_name,
    jc_name,
    noe_name,
    pre_name,
    rdc_name,
    rh_name,
    )

from spycipdb.core.calculators import *
from spycipdb.components.helpers import *

from idpconfgen.libs.libmulticore import pool_function

def selective_calculator(
    pdbfilepaths,
    exp_fp,
    modules,
    bc_errors=default_bc_errors,
    ncores=1
    ):
    """
    Back-calculate experimental data with no pre-existing back-calculations.
    
    Output of SPyCi-PDB is matched to that of `Stack` objects so no additional
    "parsing" needed.
    
    Parameters
    ----------
    pdbfilepaths : list
        List of paths that point to PDBs to perform BC on
    
    exp_fp : list
        List of experimental files to use as reference
    
    modules : list
        List of experimental modules to perform BC on
    
    bc_errors : dict
        Dictionary of prescribed back-calculator errors
    
    ncores : int
        Number of workers to initialize multiprocessing for
    
    Returns
    -------
    new_bc : dict
        Dictionary of `Stack` objects associated with back-calculated
        data for each missing module
    """
    new_bc = {}
    
    init_files(log, ".xeisd_bc")
    log.info(T(f"Starting back-calculation of necessary modules using {ncores} workers"))
    
    for exp in modules:
        log.info(S(f"back-calculating {exp} datatypes..."))
        lists = []
        if exp == noe_name:
            execute = partial(
                report_on_crash,
                calc_noe,
                exp_fp[noe_name],
            )
        elif exp == pre_name:
            execute = partial(
                report_on_crash,
                calc_pre,
                exp_fp[pre_name],
            )
        # add more `elif` statements as we test more modules
        
        execute_pool = pool_function(execute, pdbfilepaths, ncores=ncores)
        for result in execute_pool:
            lists.append(result[1])
        data = pd.DataFrame(lists)
        new_bc[exp] = Stack(exp, data, bc_errors[exp], None)
        
        log.info(S("done"))
            
    return new_bc
