"""
Module contains variable definitions and misc utility functions.

Functions and logic inspired/imported from https://github.com/THGLab/X-EISD/
"""
import os
import numpy as np

exp_idx = 'index'
exp_resnum = 'resnum'
exp_val = 'value'
exp_dist_val = 'dist_value'
exp_max = 'upper'
exp_atmID = 'atomname'
exp_min = 'lower'
exp_err = 'error'

eisd_run_all = 'all'
eisd_run_single = 'single'
eisd_run_pairs = 'pairs'
eisd_modes = (
    eisd_run_all,
    eisd_run_single,
    eisd_run_pairs
    )
default_mode = eisd_run_all

opt_max = 'max'
opt_mc = 'mc'
eisd_optimization_types = (
    opt_max,
    opt_mc,
    None,
    )
default_type = opt_max

parse_mode_exp = 'exp'
parse_mode_back = 'bc'
parse_modes = (
    parse_mode_exp,
    parse_mode_back,
    )

saxs_name = 'saxs'
cs_name = 'cs'
fret_name = 'fret'
jc_name = 'jc'
noe_name = 'noe'
pre_name = 'pre'
rdc_name = 'rdc'
rh_name = 'rh'
eisd_modules = (
    saxs_name,
    cs_name,
    fret_name,
    jc_name,
    noe_name,
    pre_name,
    rdc_name,
    rh_name,
    )

# define back calculation uncertainties
# refer to Lincoff et al. 2020 for details
default_bc_errors = {
    pre_name: 0.0001,
    noe_name: 0.0001,
    saxs_name: 0.006,
    fret_name: 0.0074,
    rh_name: 0.812,
    rdc_name: 0.88,
    # cs error reported from UCBShift, units are in ppm
    cs_name: {'C': 0.84, 'CA': 0.81, 'CB': 1.00, 'H': 0.31, 'HA': 0.19, 'N': 1.81},  # noqa: E501
    jc_name: {'A': np.sqrt(0.14), 'B': np.sqrt(0.03), 'C': np.sqrt(0.08)},
    }
default_weights = {
    pre_name: 1,
    noe_name: 1,
    saxs_name: 1,
    fret_name: 1,
    rh_name: 1,
    rdc_name: 1,
    cs_name: 1,
    jc_name: 1,
    }
jc_bc_mu = {'A': 6.51, 'B': -1.76, 'C': 1.6}

rmse_idx = 0
score_idx = 1
avg_bc_idx = 2

XEISD_TITLE = (
    "\n _  _     ____  ____  ___  ____  \n"  # noqa: W605
    "( \\/ )___( ___)(_  _)/ __)(  _ \\ \n"  # noqa: W605
    " )  ((___))__)  _)(_ \\__ \\ )(_) )\n"  # noqa: W605
    "(_/\\_)   (____)(____)(___/(____/ \n"  # noqa: W605
    "================================="
    )


def add_optimization_mode(parser):
    """Add option to choose from optimization methods."""
    parser.add_argument(
        '-m',
        '--mode',
        help=(
            'Optimization mode. Could be Metropolis Monte Carlo (mc) '
            'or score maximization method (max).'
            'Defaults to max.'
            ),
        default=opt_max,
        choices=tuple((opt_max, opt_mc)),
        )


# The following two functions have been imported from:
# https://github.com/THGLab/X-EISD/blob/master/eisd/utils/miscell.py
# https://github.com/Oufan75/X-EISD/blob/master/eisd/utils.py
def make_pairs(all):
    """Make pairs of items in a list."""
    pairs = []
    for i in range(len(all)):
        for j in range(i + 1, len(all)):
            pairs.append([all[i], all[j]])
    return pairs


def modes(mode, all):
    """Initialize a mode from a list."""
    flags = {}
    for prop in all:
        flags[prop] = False

    if mode is eisd_run_all:
        return {flag: True for flag in flags}
    elif type(mode) is list:
        for flag in mode:
            flags[flag] = True
        return flags
    elif type(mode) is str:
        flags[mode] = True
        return flags

    else:
        raise ValueError("The mode in the main function is not recognized.")


def meta_data(fpath):
    """
    Filter through exp and bc ata files.
    
    Return only the paths that exist in both cases.

    Automatically removes paths to datafiles if e.g. there are more
    back-calculated data than experimental.
    
    Parameters
    ----------
    fpath : str or Path
        Path to the directory containing both experimental
        and back-calculated data files.

    Returns
    -------
    meta : dict
        First layer keys of `exp` and `bc` having values of dictionaries
        with keys being the module and the value being the path to the data.
    
    errlog : list
        List of errors to relay to the user if there are any.
    """
    exp_paths = []
    back_paths = []
    valid_exp_modules = []
    valid_back_modules = []
    
    meta = {}
    tobc = None
    errlog = []
    
    all_files = [f for f in os.listdir(fpath) if os.path.isfile(os.path.join(fpath, f))]  # noqa: E501
    for f in all_files:
        if f.startswith(parse_mode_exp):
            if f.endswith(eisd_modules):
                exp_paths.append(os.path.join(fpath, f))
                _ext = f[f.rindex('.') + 1:]
                valid_exp_modules.append(f'.{_ext}')
        elif f.startswith(parse_mode_back):
            if f.endswith(eisd_modules):
                back_paths.append(os.path.join(fpath, f))
                _ext = f[f.rindex('.') + 1:]
                valid_back_modules.append(f'.{_ext}')
    
    valid_exp_modules.sort()
    valid_back_modules.sort()
    
    if valid_exp_modules == []:
        errlog.append(
            'WARNING: no valid experimental files found.'
            ' Please refer to the help documentation for'
            ' this module.'
            )
        return meta, errlog
    else:
        differences = tuple(set(valid_exp_modules) ^ (set(valid_back_modules)))
        
        # Do not clear `exp_paths` because we do back-calculations internally
        if differences:
            #exp_paths = [exp for exp in exp_paths if not exp.endswith(differences)]  # noqa: E265, E501
            back_paths = [bck for bck in back_paths if not bck.endswith(differences)]  # noqa: E501
            errlog.append(
                'Note: found inconsistent experimental and back-calculation'
                ' data pairs.'
                )
            errlog.append(f'Back-calculating for: {differences}...')
            tobc = [module[1:] for module in differences]
    
    EXP_DATA_FILENAMES = {}
    BACK_DATA_FILENAMES = {}
    
    for module in eisd_modules:
        for exp in exp_paths:
            if exp.endswith(f'.{module}'):
                EXP_DATA_FILENAMES[module] = exp
        for bc in back_paths:
            if bc.endswith(f'.{module}'):
                BACK_DATA_FILENAMES[module] = bc

    meta = {
        parse_mode_exp: EXP_DATA_FILENAMES,
        parse_mode_back: BACK_DATA_FILENAMES,
        }
    
    return meta, tobc, errlog
