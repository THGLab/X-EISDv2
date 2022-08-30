"""
Main interface for scoring module of X-EISD.

Citation:
(1) James Lincoff, Mojtaba Haghighatlari, Mickael Krzeminski, Joao Teixeira, Gregory Neal-Gomes,
Claudu Gardinaru, Julie Forman-Kay, Teresa Head-Gordon, https://www.nature.com/articles/s42004-020-0323-0
(2) David H. Brookes, and Teresa Head-Gordon, Experimental Inferential Structure Determination of
Ensembles for Intrinsically Disordered Proteins, JACS 138, 2016, 4530-4538 DOI: 10.1021/jacs.6b00351

USAGE:
    $ xeisd score [--exp-files] [--back-files]
    $ xeisd score [--exp-files] [--pdb-structures] [--output] [--ncores]

OUTPUT:

"""
import argparse
import os
import shutil
from functools import partial

from xeisd import Path, log
from xeisd.libs import libcli
from xeisd.logger import S, T, init_files, report_on_crash

#TODO: idpconfgen imports and installation requirements

LOGFILESNAME = '.xeisd_score'
TMPDIR = '__tmpscore__'

_name = 'score'
_help = 'Score conformational ensembles against experimental data.'

_prog, _des, _usage = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_usage,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )
