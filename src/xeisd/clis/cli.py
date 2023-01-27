"""
X-EISDv2.

Command-line interface for scoring and reweighting functions
for protein ensembles. Core logic based on X-EISD: extended
experimental inferential structure determination.

Please see: https://github.com/THGLab/X-EISD
for reference and citing.

USAGE:
    For help:
    $ xeisd -h
"""
import argparse
import sys

from xeisd import __version__, log
from xeisd.clis import cli_optimize, cli_score, cli_subset
from xeisd.libs import libcli
from xeisd.logger import S


_prog, _description, _usageage = libcli.parse_doc_params(__doc__)

description = f"""
{_description}

Core functions:
    * {cli_score._name}
    * {cli_optimize._name}
    * {cli_subset._name}
"""

ap = libcli.CustomParser(
    prog='xeisd',  # _prog,
    description=libcli.detailed.format(description),
    usage=_usageage,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_version(ap)

subparsers = ap.add_subparsers(
    title='X-EISD routines',
    help='Short description:',
    )

libcli.add_subparser(subparsers, cli_score)
libcli.add_subparser(subparsers, cli_optimize)
libcli.add_subparser(subparsers, cli_subset)


def load_args():
    """Load user input arguments."""
    return ap.parse_args()


def maincli():
    """
    Execute subroutine.
    
    Arguments are read from user command line input.
    """
    # prints help if not arguments are passed
    # if >2 prints help of subroutines.
    if len(sys.argv) < 2:
        ap.print_help()
        ap.exit()

    cmd = load_args()

    with open('xeisd.version', 'w') as fout:
        fout.write(f'version: {__version__}')

    cmd.func(**vars(cmd))
    log.info(S('finished properly'))


if __name__ == '__main__':
    maincli()
