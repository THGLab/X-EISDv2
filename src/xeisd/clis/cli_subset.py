r"""
Interface for subsetting ensembles from optimization data.

Requires the `indices.csv` result file as well as the path to the folder
of PDBs used in obtaining the back-calculation data.

USAGE:
    $ xeisd subset \
        [--indices-file] \
        [--pdb-structures] \
        [--num-ensembles] \
        [--output-folder]
"""
import argparse
import shutil
import numpy as np

from idpconfgen.libs.libio import (
    extract_from_tar,
    make_folder_or_cwd,
    read_path_bundle,
)
from natsort import os_sorted

from xeisd import Path, log
from xeisd.libs import libcli
from xeisd.logger import S, T, init_files

LOGFILESNAME = '.xeisd_subset'
TMPDIR = '__tmpsubset__'

_name = 'subset'
_help = 'Subset initial pool ensembles with indices data from optimize.'

_prog, _des, _usage = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_usage,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )


