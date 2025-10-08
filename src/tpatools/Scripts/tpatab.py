import tpatools
import argparse
from pathlib import Path
import sys
from tpatools.parse import tpa_table
from tpatools.tools import filepath_searcher

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'filepath', 
        default=None,
        nargs="?",
        help='Filepath containing a turbomole 2PA calculation output file. If no filepath is specified, will search for files named escf.out, bse.out, ricc2.out, or tpa.out within the current directory',
    )
    parser.add_argument(
        '-v', 
        '--verbose',
        action='store_true',
        help='Print excitation energies and transition strengths in a.u.',
    )
    parser.add_argument(
        '-c', 
        '--compact',
        action='store_true',
        help='Print output names in a more compact format',
    )
    parser.add_argument(
        '-t',
        '--tex',
        action='store_true',
        help='Print output in LaTeX format',
    )
    parser.add_argument(
        '-r',
        '--round',
        type=int,
        default=3,
        help="number of digits to round decimals to - default is 3",
    )
    parser.add_argument(
        '-n',
        '--noirrep',
        action='store_false',
        help="whether or not to include the irreducible representation in state labels (enabled by default)",
    )
    parser.add_argument(
        '-m',
        '--mult',
        action='store_true',
        help="whether or not to include the multiplicity in labels (ex S1) - disabled by default",
    )

    args = parser.parse_args()

    filepath = filepath_searcher(args.filepath)
    df = tpa_table(filepath, irrep=args.noirrep, mult=args.mult, compact=args.compact)
    if args.verbose == False:
        df.drop(
            columns=['Excitation Energy /eV', '2PA Strength /a.u.'],
            inplace=True,
        )

    if args.tex == True:
        print(
            df.to_latex(
                index=False,
                float_format=f'%.{args.round}f',
            )
        )
    else:
        print(
            df.to_string(
                index=False,
                float_format=f'%.{args.round}f',
            )
        )




    

