r"""
Main interface for optimization module of X-EISD.

Please note that if you provide a path to PDB structures AND
back-calculated files, X-EISD will attempt to parse the back-calculated
files first and if they're not in the right format, we will back-calculate
the information for you using SPyCi-PDB.

For `--data-files`, please name all your experimental files
ending with the extension(s) corresponding to the experiment
of interest and starting with `exp_`.
    
    For example:
    exp_*.saxs, exp_*.cs, exp_*.fret, exp_*.jc,
    exp_*.noe, exp_*.pre, exp_*.rdc, exp_*.rh
    
For back-calculated files (if they're already in the correct format),
please start them with the `bc_` prefix, ending with the same file
extensions as seen above.
The internal formatting can be .JSON however.
    
    For example:
    bc_*.saxs, bc_*.cs, bc_*.fret, etc.

USAGE:
    $ xeisd optimize \
        [--data-files] \
        [--pdb-structures] \
        [--nconfs] \
        [--nres] \
        [--custom-error] \
        [--epochs] \
        [--mode] \
        [--num-exchange] \
        [--mc-beta] \
        [--output] \
        [--ncores]

OUTPUT:

Please provide an output directory for the 3 .CSV files to be saved.
The three files will include final optimized results, indicies, and best J-coupling scores.
"""
import argparse
import json
import shutil
from functools import partial

import pandas as pd
from idpconfgen.libs.libio import extract_from_tar, read_path_bundle
from idpconfgen.libs.libmulticore import pool_function
from idpconfgen.libs.libstructure import Structure

from xeisd import Path, log
from xeisd.components import (
    add_optimization_mode,
    XEISD_TITLE,
    avg_bc_idx,
    cs_name,
    default_bc_errors,
    exp_atmID,
    exp_idx,
    exp_resnum,
    meta_data,
    parse_mode_back,
    parse_mode_exp,
    rmse_idx,
    saxs_name,
    score_idx,
    )
from xeisd.components.helper import (
    return_indices_of_bc_saxs,
    selective_calculator,
    )
from xeisd.components.optimizer import XEISD
from xeisd.components.parser import parse_bc_errors, parse_data
from xeisd.libs import libcli
from xeisd.logger import S, T, init_files, report_on_crash


LOGFILESNAME = '.xeisd_optimize'
TMPDIR = '__tmpoptimize__'

_name = 'optimize'
_help = 'Optimize conformational ensembles against experimental data.'

_prog, _des, _usage = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_usage,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_argument_data_files(ap)
libcli.add_argument_pdb_files(ap)
libcli.add_argument_number_conformers(ap)
libcli.add_argument_number_residues(ap)
libcli.add_argument_custom_bc_error(ap)
libcli.add_argument_epochs(ap)
add_optimization_mode(ap)

ap.add_argument(
    '-nex',
    '--num-exchange',
    help=(
        'Number of conformer exchange attempts. '
        'Defaults to 100.'
        ),
    default=100,
    type=int,
)

ap.add_argument(
    '-mct',
    '--mc-beta',
    help=(
        'Temperature parameter for Metropolis Monte Carlo optimization mode. '
        'Defaults to 0.1.'
        ),
    default=0.1,
    type=float,
)

libcli.add_argument_output(ap)
libcli.add_argument_ncores(ap)

ap.add_argument(
    '--tmpdir',
    help=(
        'Temporary directory to store data during calculation '
        'if needed.'
        ),
    type=Path,
    default=TMPDIR,
    )


def main(
        data_files,
        nconfs,
        nres,
        output,
        epochs,
        num_exchange=100,
        mc_beta=0.1,
        custom_error=None,
        pdb_files=None,
        ncores=1,
        tmpdir=TMPDIR,
        func=None,
        ):
    """
    Processes and optimizes ensembles.
    
    Parameters
    ----------
    data_files : str or Path, required
        Path to the folder containing experimental and back-calc
        data files.
    
    nconfs : int, required
        Number of conformations in ensemble.
    
    nres : int required
        Number of residues in protein.
    
    pdb_files : str or Path, optional
        Path to PDB files on the disk. Accepts TAR file.
    
    custom_error : str or Path, optional
        Path to the file containing custom back-calculator errors.
    
    output : str or Path, optional
        Path to the folder to store eisd outputs.
        
    ncores : int, optional
        Number of workers to use for multiprocessing.
        Defaults to 1.
    
    num_exchange : int, optional
        Number of conformer exchange attempts.
        Defaults to 100.
    
    epochs : int, required
        Number of times to run main optimization.
        Recommended at lesat the size of your ensemble.
    
    mc_beta : float, optional
        Temperature parameter for the Metropolis Monte Carlo
        optimization mode.
    """


if __name__ == '__main__':
    libcli.maincli(ap, main)
