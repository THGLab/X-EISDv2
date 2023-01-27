r"""
Interface for subsetting ensembles from optimization data.

Requires the `indices.csv` result file as well as the path to the folder
of PDBs used in obtaining the back-calculation data.

USAGE:
    $ xeisd subset \
        [--indices-file] \
        [--pdb-structures] \
        [--num-ensembles] \
        [--random-seed] \
        [--output-folder]
"""
import argparse
import os
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

ap.add_argument(
    '-i',
    '--indices',
    help='Path to indices.csv file from xeisd optimize subclient.',
    required=True,
    type=str,
    )

libcli.add_argument_pdb_files(ap)

ap.add_argument(
    '-ne',
    '--num-ensembles',
    help='Number of ensembles to obtain from subsetting.',
    required=True,
    type=int,
    )

libcli.add_argument_random_seed(ap)
libcli.add_argument_output_folder(ap)


def main(
        indices,
        pdb_files,
        num_ensembles,
        random_seed=0,
        output_folder=None,
        tmpdir=TMPDIR,
        func=None,
        ):
    """
    Subsets ensembles given indices.csv output and initial pool of interest.
    
    Parameters
    ----------
    indices : str or Path, required
        Path to the indices.csv file generated from optimize subclient.
    
    pdb_files : str or Path, required
        Path to PDB files on the disk. Accepts TAR file.
    
    num_ensembles : int, required
        Integer number of ensembles to subset from indices.csv.
    
    random_seed : int, optional
        Random seed for reproducibility of ensemble selection from indices.csv.
        Defaults to 0.
        
    output_folder : str or Path, optional
        Path to folder to store ensembles.
        Defaults to current working directory.
    """
    np.random.seed(random_seed)
    output_folder = make_folder_or_cwd(output_folder)
    init_files(log, Path(output_folder, LOGFILESNAME))
    
    log.info(T('Reading conformer ensemble paths'))
    try:
        pdbs2operate = extract_from_tar(pdb_files, output=tmpdir, ext='.pdb')  # noqa: E501
        _istarfile = True
    except (OSError, TypeError):
        pdbs2operate = list(read_path_bundle(pdb_files, ext='pdb'))
        _istarfile = False
    pdbs2operate = os_sorted(pdbs2operate)
    log.info(S('done'))
    
    indices_list = []
    with open(indices) as ifile:
        for line in ifile:
            line = line.strip().split(',')
            int_line = [eval(i) for i in line]
            indices_list.append(int_line)
    
    if num_ensembles > len(indices_list):
        log.info(
            "Warning, your number of desired ensembles must be "
            "less than or equal to the number of lines in indices.csv"
            )
        log.info(f"Setting `-ne` to {len(indices_list)}.")
    
    i_arr = np.array(indices_list)
    selected_ens = i_arr[np.random.choice(len(i_arr), size=num_ensembles, replace=False)]  # noqa: E501
    
    for i in range(num_ensembles):
        path = os.path.join(output_folder, f"ensemble_{i + 1}")
        try:
            os.mkdir(path)
        except FileExistsError:
            pass
        
        for j in selected_ens[i]:
            shutil.copy2(pdbs2operate[j], path)
        
    if _istarfile:
        shutil.rmtree(tmpdir)
    
    return


if __name__ == '__main__':
    libcli.maincli(ap, main)
