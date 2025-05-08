import tpatools
import argparse
from pathlib import Path
import sys
from tpatools.parse import gather_state_data

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-d',
        '--basedir',
        default='./',
        help='directory to search in',
    )
    parser.add_argument(
        '-o',
        '--escfname',
        default='escf.out',
        help='name of the escf output file(s)',
    )
    parser.add_argument(
        '-e',
        '--egradname',
        default='egrad.out',
        help='name of the grad output file(s)',
    )
    parser.add_argument(
        '-t',
        '--tex',
        action='store_true',
        help='Print output in LaTeX format',
    )
    parser.add_argument(
        '-f',
        '--fullnames',
        action='store_true',
        help='include full paths to files in column',
    )
    parser.add_argument(
        '-l',
        '--latexnames',
        action='store_true',
        help='name parameters using latex-ready notation',
    )
    parser.add_argument(
        '-r',
        '--round',
        type=int,
        default=6,
        help="number of digits to round decimals to - default is 3",
    )
    parser.add_argument(
        '-n',
        '--stateno',
        type=int,
        default=1,
        help='excited state number'
    )
    parser.add_argument(
        '-sort',
        '--sortlist', 
        default=None,
        nargs="?",
        help='filepath pointing to a text file containing the desired keys/filenames sorted in the order for the output table',
    )
    args = parser.parse_args()

    if args.sortlist is not None:
        with Path(args.sortlist).open() as sortfile:
            orderedkeys = [x.strip() for x in sortfile if x.strip != '']
    else:
        orderedkeys = None

    df = gather_state_data(
        args.basedir,
        outfilename = args.escfname,
        egradoutname = args.egradname,
        state=args.stateno,
        egradavail=False,
        fulldirnames = args.fullnames,
        suppress_egrad_warning = True,
        tabulate = True,
        latexnames=args.latexnames,
        orderedkeys=orderedkeys,
    )
    if args.tex == True:
        print(
            df.to_latex(
                float_format=f'%.{args.round}f',
            )
        )
    else:
        print(
            df.to_string(
                float_format=f'%.{args.round}f',
            )
        )
