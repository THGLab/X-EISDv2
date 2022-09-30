"""
Help X-EISD perform calculations with third-party programs.

Initial creation made with integrated SPyCi-PDB calculators in mind.
"""
from functools import partial

import pandas as pd
from idpconfgen.libs.libmulticore import pool_function
from spycipdb.components.helpers import *  # noqa: F403
from spycipdb.core.calculators import *  # noqa: F403

from xeisd import log
from xeisd.components import (
    default_bc_errors,
    exp_idx,
    exp_val,
    fret_name,
    jc_bc_mu,
    jc_name,
    noe_name,
    pre_name,
    saxs_name,
    )
from xeisd.components.parser import Stack
from xeisd.logger import S, T, init_files, report_on_crash


def return_indices_of_bc_saxs(exp_idx, bc_idx):
    """
    Truncate indexes from back-calculated data to match exp data.
    
    Assumes back-calculated list of indicies are longer

    Parameters
    ----------
    exp_idx : list
        List of indicies from experiment.
    
    bc_idx : list
        List of indices from back-calculator.
    
    Returns
    -------
    to_keep : list of int
        List of positions of indices to keep
        to match experiment data.
    """
    exp_set = set(exp_idx)
    
    return [i for i, v in enumerate(bc_idx) if v in exp_set]


def selective_calculator(
        pdbfilepaths,
        exp_fp,
        modules,
        bc_errors=default_bc_errors,
        ncores=1
        ):
    """
    Back-calculate experimental data with no pre-existing back-calculations.
    
    Output of SPyCi-PDB is matched to that of `Stack` objects.
    
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
    log.info(T(f"Starting back-calculation of necessary modules using {ncores} workers"))  # noqa: E501
    
    for exp in modules:
        log.info(S(f"back-calculating {exp} datatypes..."))
        lists = []
        if exp == noe_name:
            execute = partial(
                report_on_crash,
                calc_noe,  # noqa: F405
                exp_fp[noe_name],
                )
        elif exp == pre_name:
            execute = partial(
                report_on_crash,
                calc_pre,  # noqa: F405
                exp_fp[pre_name],
                )
        elif exp == jc_name:
            execute = partial(
                report_on_crash,
                calc_jc,  # noqa: F405
                exp_fp[jc_name],
                )
        elif exp == fret_name:
            execute = partial(
                report_on_crash,
                calc_smfret,  # noqa: F405
                exp_fp[fret_name],
                )
        elif exp == saxs_name:
            execute = partial(
                report_on_crash,
                crysol_helper,  # noqa: F405
                lm=20,  # default number of harmonics for crysolv3
                )
            str_pdbfp = []
            for path in pdbfilepaths:
                str_pdbfp.append(str(path))
        # add more `elif` statements as we test more modules
        
        if exp == saxs_name:
            execute_pool = pool_function(execute, str_pdbfp, ncores=ncores)
        else:
            execute_pool = pool_function(execute, pdbfilepaths, ncores=ncores)
        
        for result in execute_pool:
            if exp == saxs_name:
                lists.append(result[1][exp_val])
                saxs_bc_idx = result[1][exp_idx]
            else:
                lists.append(result[1])
        data = pd.DataFrame(lists)
        
        if exp == jc_name:
            new_bc[exp] = Stack(exp, data, bc_errors[exp], jc_bc_mu)
        # saving saxs indexes to mu variable in stack
        elif exp == saxs_name:
            new_bc[exp] = Stack(exp, data, bc_errors[exp], saxs_bc_idx)
        else:
            new_bc[exp] = Stack(exp, data, bc_errors[exp], None)
        
        log.info(S("done"))
            
    return new_bc
