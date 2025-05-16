from tpatools.plot import tpaplot
import numpy as np
import argparse
import sys
from tpatools.parse import tpa_table
from pathlib import Path
from tpatools.tools import filepath_searcher

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'filepath', 
        default=None,
        nargs="?",
        help='Filepath containing a turbomole 2PA calculation output file. If no filepath is specified, will search for files named escf.out, bse.out, ricc2.out or tpa.out within the current directory',
    )
    parser.add_argument(
        '-b',
        '--broadening',
        default=0.1,
        type=float,
        help='Lifetime broadening (in eV) for simulated spectrum. Default 0.1',
    )
    parser.add_argument(
        '-f',
        '--xmin',
        default=None,
        type=float,
        help='Minimum value on x axis (in eV)',
    )
    parser.add_argument(
        '-t',
        '--xmax',
        default=None,
        type=float,
        help='Maximum value on x axis (in eV)',
    )
    parser.add_argument(
        '-s',
        '--savefile',
        default=None,
        help='An output filename to save the generated plot (no file saved if not specified)'
    )
    parser.add_argument(
        '-n',
        '--nm',
        default=False,
        action='store_true',
        help='Whether to convert units from eV to nm (false by default)',
    )
    args = parser.parse_args()
 
    filepath = filepath_searcher(args.filepath)
    tab = tpa_table(filepath)

    if args.xmin is None:
        xmin = np.min(tab['Excitation Energy /eV']) - 1
    else:
        xmin = args.xmin

    if args.xmax is None:
        xmax = np.max(tab['Excitation Energy /eV']) + 1
    else:
        xmax = args.xmax

    tpaplot(
        tab,
        width=args.broadening,
        xmin = xmin,
        xmax= xmax,
        nm=args.nm,
        save=args.savefile,
    )

