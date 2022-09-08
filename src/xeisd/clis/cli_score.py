"""
Main interface for scoring module of X-EISD.

Please note that if you provide a path to PDB structures AND
back-calculated files, X-EISD will attempt to parse the back-calculated
files first and if they're not in the right format, we will back-calculate
the information for you using SPyCi-PDB.

For `--exp-files`, please name all your experimental files
ending with the extension(s) corresponding to the experiment
of interest and starting with `exp_`.
    
    For example:
    exp_*.saxs, exp_*.cs, exp_*.fret, exp_*.jc,
    exp_*.noe, exp_*.pre, exp_*.rdc, exp_*.rh
    
For back-calculated files (if they're already in the correct format),
please start them with the `bc_` prefix, ending with the same file extensions as seen above.
The internal formatting can be .JSON however.
    
    For example:
    bc_*.saxs, bc_*.cs, bc_*.fret, etc.

USAGE:
    $ xeisd score [--exp-files] [--back-files] [--epochs]
    $ xeisd score [--exp-files] [--pdb-structures] [--epochs] [--output] [--ncores]

OUTPUT:

"""
import argparse
import os
from select import select
import shutil
from functools import partial

from xeisd import Path, log
from xeisd.libs import libcli
from xeisd.logger import S, T, init_files, report_on_crash

from xeisd.components import (
    default_bc_errors,
    parse_mode_exp,
    parse_mode_back,
    meta_data,
    )
from xeisd.components.parser import parse_data, parse_bc_errors
from xeisd.components.optimizer import XEISD
from xeisd.components.helper import selective_calculator

from idpconfgen.libs.libstructure import Structure
from idpconfgen.libs.libio import extract_from_tar, read_path_bundle
from idpconfgen.libs.libmulticore import pool_function

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

libcli.add_argument_data_files(ap)
libcli.add_argument_pdb_files(ap)
libcli.add_argument_custom_bc_error(ap)
libcli.add_argument_epochs(ap)
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
        epochs,
        custom_error=None,
        pdb_files=None,
        output=None,
        ncores=1,
        tmpdir=TMPDIR,
        func=None,
    ):
    """
    Process and score ensembles.

    Parameters
    ----------
    exp_files : str or Path, required
        Path to the folder containing experimental
        data files.
        
    pdb_files : str or Path, optional
        Path to PDB files on the disk. Accepts TAR file.
        
    back_files : str or Path, optional
        Path to the folder containing experimental
        data files.
    
    custom_error : str or Path, optional
        Path to the file containing custom back-calculator errors.
    
    epochs : int, required
        Number of times to run main optimization.
    
    output : str or Path, optional
        Path to the folder to store eisd outputs.
        Defaults to working directory.
        
    ncores : int, optional
        Number of workers to use for multiprocessing.
        Defaults to 1.   
    """
    init_files(log, LOGFILESNAME)
    
    #if pdb_files == None and back_files == None:
    #    log.info(S('WARNING: you must provide either PDB files or back-calculated files.'))
    #    return
    
    log.info(T('Reading conformer ensemble paths'))
    _istarfile = False
    if pdb_files:
        try:
            pdbs2operate = extract_from_tar(pdb_files, output=tmpdir, ext='.pdb')
            _istarfile = True
        except (OSError, FileNotFoundError, TypeError):
            pdbs2operate = list(read_path_bundle(pdb_files, ext='pdb'))
            _istarfile = False
        log.info(S('done'))
        
        ens_size = len(pdbs2operate)
        
        s = Structure(pdbs2operate[0])
        s.build()
        nres = len(s.residues)
    
    log.info(T('Checking experimental and back-calculation data files'))
    filenames, tobc, _meta_errs = meta_data(data_files)
    if _meta_errs:
        for e in _meta_errs:
            log.info(S(e))
        if filenames == {}:
            log.info(S('done'))
            return
    log.info(S('done'))

    log.info(T('Parsing experimental and back-calculated files'))
    exp_data, _parse_errs = parse_data(filenames[parse_mode_exp], mode=parse_mode_exp)
    if _parse_errs:
        for e in _meta_errs:
            log.info(S(e))
    
    if custom_error:
        bc_errors = parse_bc_errors(custom_error)
    else:
        bc_errors = default_bc_errors
    
    #TODO: perform back-calculation via SPyCi-PDB
    # - make temporary files with the back-calculations to be removed later in tmpdir
    # - back_paths will point to files within tmpdir, or another directory
    # this way, we can save on the parsing steps
    if tobc:
        new_back_data = selective_calculator(pdbs2operate, tobc)
    if filenames[parse_mode_back] != {}:
        existing_back_data, _parse_errs = \
            parse_data(filenames[parse_mode_back], parse_mode_back, bc_errors)
    
    back_data = existing_back_data + new_back_data

    eisd_ens = XEISD(exp_data, back_data, ens_size, nres)
    log.info(T(f'Starting X-EISD Scoring'))
    execute = partial (
        report_on_crash,
        eisd_ens.calc_scores,
        epochs=epochs,
        ens_size=ens_size,        
        )
    execute_pool = pool_function(execute, ncores=ncores)
    #TODO: print output to stdout or save to disk
    
    log.info(S('done'))
    
    if _istarfile:
        shutil.rmtree(tmpdir)
        
    log.info(S('done'))
    
    return


if __name__ == '__main__':
    libcli.maincli(ap, main)
