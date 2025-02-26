import tpatools
import argparse
from pathlib import Path
import sys
from tpatools.parse import escf_table

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'filepath', 
        default=None,
        nargs="?",
        help='Filepath containing a turbomole escf 2PA calculation output file. If no filepath is specified, will search for files named escf.out, bse.out, or tpa.out within the current directory',
    )
    parser.add_argument(
        '-v', 
        '--verbose',
        action='store_true',
        help='Print excitation energies and transition strengths in a.u.',
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
    args = parser.parse_args()


    filepath = None
    if args.filepath is None:
        potential_names = ['escf.out', 'bse.out', 'tpa.out']
        for name in potential_names:
            if Path(name).is_file():
                filepath = Path(name)
                break
        
        if filepath is None: 
            print(f'No output filename has been provided. As a fallback, the program has checked and found no files from the following default namelist in the directory:\n\n{"\n".join(potential_names)}\n\nPlease retry and specify the name of your output file' )
            sys.exit()
    else:
        if Path(args.filepath).is_file():
            filepath = Path(args.filepath)
        else:
            print('The provided filepath does not exist, please check your input and try again')
            sys.exit()

    df = escf_table(filepath)
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




    

