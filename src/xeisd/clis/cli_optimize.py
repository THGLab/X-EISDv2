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
        [--mode] \
        [--num-exchange] \
        [--mc-temperature] \
        [--epochs] \
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