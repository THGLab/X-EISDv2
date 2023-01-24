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
The three files will include final optimized results, indicies, and
best J-coupling scores.
"""
import argparse
import json
import shutil

import pandas as pd
from idpconfgen.libs.libio import extract_from_tar, read_path_bundle
from idpconfgen.libs.libstructure import Structure

from xeisd import Path, log
from xeisd.components import (
    XEISD_TITLE,
    add_optimization_mode,
    cs_name,
    default_bc_errors,
    exp_atmID,
    exp_idx,
    exp_resnum,
    meta_data,
    opt_max,
    parse_mode_back,
    parse_mode_exp,
    saxs_name,
    )
from xeisd.components.helper import (
    return_indices_of_bc_saxs,
    selective_calculator,
    )
from xeisd.components.optimizer import XEISD
from xeisd.components.parser import parse_bc_errors, parse_data
from xeisd.libs import libcli
from xeisd.logger import S, T, init_files


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
        mode=opt_max,
        num_exchange=100,
        mc_beta=0.1,
        custom_error=None,
        pdb_files=None,
        ncores=1,
        tmpdir=TMPDIR,
        func=None,
        ):
    """
    Process and optimizes ensembles.
    
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
    init_files(log, LOGFILESNAME)
    back_data = None
    
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
            pdbs2operate = extract_from_tar(pdb_files, output=tmpdir, ext='.pdb')  # noqa: E501
            _istarfile = True
        except (OSError, TypeError):
            pdbs2operate = list(read_path_bundle(pdb_files, ext='pdb'))
            _istarfile = False
        log.info(S('done'))
        
        nconfs = len(pdbs2operate)
        
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
    
    # START----CHECKS SPECIFIC MODULES------

    # Checks if back-calculated SAXS contains more data points than experiment
    if saxs_name in filenames[parse_mode_back]:
        with open(filenames[parse_mode_back][saxs_name], 'r') as f:
            raw = json.loads(f.read())
            bc_idx = raw['format']
            exp_indices = exp_data[saxs_name].data[exp_idx].values.tolist()
            
            if len(bc_idx) != len(exp_indices):
                aligned = return_indices_of_bc_saxs(exp_indices, bc_idx)
                back_data[saxs_name].data = \
                    back_data[saxs_name].data.iloc[:, aligned]

    # Checks if back-calculated CS data contains more atom names than exp.
    # Also aligns back-calc data to residues and atoms from exp
    if cs_name in filenames[parse_mode_back]:
        with open(filenames[parse_mode_back][cs_name], 'r') as f:
            raw = json.loads(f.read())
            res_list = raw['format']['res']
            raw.pop('format', None)
            
            exp_cs = exp_data[cs_name].data
            realigned_bc = []
            
            for conf in raw:
                per_conf = []
                cs_info = raw[conf]
                
                for _idx, row in exp_cs.iterrows():
                    resnum = row[exp_resnum]
                    atmid = row[exp_atmID]
                    res_idx = res_list.index(resnum)
                    per_conf.append(cs_info[atmid][res_idx])
                
                realigned_bc.append(per_conf)
            
            back_data[cs_name].data = pd.DataFrame(realigned_bc)
    
    # END-----CHECKS SPECIFIC MODULES-------
    
    if tobc:
        try:
            new_back_data = selective_calculator(
                pdbs2operate,
                filenames[parse_mode_exp],
                tobc,
                bc_errors,
                ncores,
                )
            
            if saxs_name in new_back_data:
                bc_idx = new_back_data[saxs_name].mu
                exp_indices = exp_data[saxs_name].data[exp_idx].values.tolist()
                if len(bc_idx) != len(exp_indices):
                    aligned = return_indices_of_bc_saxs(exp_indices, bc_idx)
                    new_back_data[saxs_name].data = \
                        new_back_data[saxs_name].data.iloc[:, aligned]
            
            if back_data:
                back_data = {**back_data, **new_back_data}
            else:
                back_data = new_back_data
        except UnboundLocalError:
            log.info(S('You must provide an ensemble of PDBs to perform back-calculations.'))  # noqa: E501
            return
    
    # TODO: make optimization multiprocessing
    # also why is it taking so long when scoring is really fast?
    eisd_ens = XEISD(exp_data, back_data, nres, nconfs)
    log.info(T('Starting Optimization Function'))
    results, indices, best_jcoups = \
        eisd_ens.optimize(
            epochs=epochs,
            ens_size=nconfs,
            opt_type=mode,
            beta=mc_beta,
            iters=num_exchange,
            )
    std_out = XEISD_TITLE + "\n\n"  # noqa: F841
    
    if _istarfile:
        shutil.rmtree(tmpdir)
    
    return


if __name__ == '__main__':
    libcli.maincli(ap, main)
