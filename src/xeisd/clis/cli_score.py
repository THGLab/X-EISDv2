"""
Main interface for scoring module of X-EISD.

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
please start them with the `bc_` prefix, ending with the same file extensions as seen above.
The internal formatting can be .JSON however.
    
    For example:
    bc_*.saxs, bc_*.cs, bc_*.fret, etc.

USAGE:
    $ xeisd score [--data-files] [--epochs]
    $ xeisd score [--data-files] [--pdb-structures] [--epochs] [--output] [--ncores]

OUTPUT:

"""
import json
import argparse
import shutil
from functools import partial

from xeisd import Path, log
from xeisd.libs import libcli
from xeisd.logger import S, T, init_files, report_on_crash

from xeisd.components import (
    default_bc_errors,
    parse_mode_exp,
    parse_mode_back,
    rmse_idx,
    score_idx,
    avg_bc_idx,
    XEISD_TITLE,
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
libcli.add_argument_number_conformers(ap)
libcli.add_argument_number_residues(ap)
libcli.add_argument_custom_bc_error(ap)
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
        custom_error=None,
        pdb_files=None,
        ncores=1,
        tmpdir=TMPDIR,
        func=None,
    ):
    """
    Process and score ensembles.

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
    """
    init_files(log, LOGFILESNAME)
    
    # Sets back-calculator errors
    if custom_error:
        bc_errors = parse_bc_errors(custom_error)
    else:
        bc_errors = default_bc_errors
        
    
    # Processes PDB files if given
    _istarfile = False
    if pdb_files:
        log.info(T('Reading conformer ensemble paths'))
        try:
            pdbs2operate = extract_from_tar(pdb_files, output=tmpdir, ext='.pdb')
            _istarfile = True
        except (OSError, TypeError):
            pdbs2operate = list(read_path_bundle(pdb_files, ext='pdb'))
            _istarfile = False
        log.info(S('done'))
        
        ens_size = len(pdbs2operate)
        
        s = Structure(pdbs2operate[0])
        s.build()
        nres = len(s.residues)
    
    
    # Check to see experimental and back-calc files match
    log.info(T('Checking experimental and back-calculation data files'))
    filenames, tobc, _meta_errs = meta_data(data_files)
    if _meta_errs:
        for e in _meta_errs:
            log.info(S(str(e)))
        if filenames == {}:
            log.info(S('done'))
            return
    log.info(S('done'))


    # Parse experimental and given back-calculated files
    log.info(T('Parsing experimental and back-calculated files'))
    exp_data, _parse_errs = parse_data(
        filenames[parse_mode_exp],
        mode=parse_mode_exp
    )
    if _parse_errs:
        for e in _parse_errs:
            log.info(S(str(e)))
    if filenames[parse_mode_back]:
        back_data, _parse_errs = parse_data(
            filenames[parse_mode_back],
            parse_mode_back,
            bc_errors,
        )
        if _parse_errs:
            for e in _parse_errs:
                log.info(S(str(e)))
    log.info(S('done'))    
    
    #TODO: perform back-calculation via SPyCi-PDB
    if tobc:
        try:
            new_back_data = selective_calculator(pdbs2operate, tobc, ncores)
            back_data = back_data + new_back_data
        except UnboundLocalError:
            log.info(S('You must provide an ensemble of PDBs to perform back-calculations.'))
            return
    
    dtypes = [module for module in exp_data]
    eisd_ens = XEISD(exp_data, back_data, nres, nconfs)
    log.info(T(f'Starting Scoring Function'))
    execute = partial (
        report_on_crash,
        eisd_ens.calc_scores,
        ens_size=nconfs,
        )
    execute_pool = pool_function(execute, dtypes, ncores=ncores)
    _output = {}
    for result in execute_pool:
        _output[result[0]] = {
            "rmse": result[1][rmse_idx],
            "score": result[1][score_idx],
            "avg_bc_value": result[1][avg_bc_idx].tolist(),
        }
    log.info(S('done'))
    
    log.info(T('Writing output onto disk'))
    with open(output, mode='w') as fout:
        fout.write(json.dumps(_output, indent=4))
    log.info(S('done'))
    
    std_out = XEISD_TITLE + "\n\n"
    for module in _output:
        std_out += f"{module.upper()} Score = {_output[module]['score']}\n"
        std_out += f"{module.upper()} RMSE  = {_output[module]['rmse']}\n"
        std_out += "---\n"
    
    log.info(S(std_out))
    
    if _istarfile:
        shutil.rmtree(tmpdir)
    
    return


if __name__ == '__main__':
    libcli.maincli(ap, main)
