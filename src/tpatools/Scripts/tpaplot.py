from tpatools.interactive import tpaplot
import numpy as np
import argparse
import sys
from tpatools.parse import escf_table
from pathlib import Path

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'filepath', 
        default=None,
        nargs="?",
        help='Filepath containing a turbomole escf 2PA calculation output file. If no filepath is specified, will search for files named escf.out, bse.out, or tpa.out within the current directory',
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
 
    filepath = None
    fallback_names = ['escf.out', 'bse.out', 'tpa.out']
    if args.filepath is None:
        for name in fallback_names:
            if Path(name).is_file():
                filepath = Path(name)
                break
        
        if filepath is None: 
            print(f'No output filename has been provided. As a fallback, the program has checked and found no files from the following default namelist in the directory:\n\n{"\n".join(fallback_names)}\n\nPlease retry and specify the name of your output file' )
            sys.exit()
    else:
        if Path(args.filepath).is_file():
            filepath = Path(args.filepath)
        else:
            for name in fallback_names:
                if (Path(args.filepath) / name).is_file():
                    filepath = Path(args.filepath) / name
                    break

            if filepath is None:
                print('The provided filepath does not exist, please check your input and try again')
                sys.exit() 

    tabdict = {'table' : escf_table(filepath)}

    if args.xmin is None:
        xmin = np.min(tabdict['table']['Excitation Energy /eV']) - 1
    else:
        xmin = args.xmin

    if args.xmax is None:
        xmax = np.max(tabdict['table']['Excitation Energy /eV']) + 1
    else:
        xmax = args.xmax

    tpaplot_multi(
        tabdict,
        width=args.broadening,
        x_offset=0,
        y_offset=0,
        xmin = xmin,
        xmax= xmax,
        fromentry=1,
        toentry=1,
        justone=True,
        show_y=True,
        nm=args.nm,
        save=args.savefile,
        showlegend=False,
        show_labels=False,
    )

