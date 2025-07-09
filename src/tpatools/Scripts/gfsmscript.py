from pathlib import Path
import pandas as pd
from tpatools.gfsm import gfsm, extract_gfsm_dipole_data
import sys
import argparse
import numpy as np


def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-s',
        '--states',
        nargs='+',
        default=['0', '1'],
        help=(
            'States to be used in the GFSM. i.e, "0 1" for a two state model to the'
            + ' first excited state, "0 1 4" for a three state model to the fourth '
            + 'excited state, etc. By default will perform a 2SM to the first excited state'
        )
    )
    parser.add_argument(
        '-d',
        '--basedir',
        default='./',
        help='directory containing the egrad output file. Defaults to the current directory',
    )
    parser.add_argument(
        '-f',
        '--egradname',
        default='egrad.out',
        help='Name of the egrad logfile. Defaults to egrad.out',
    )
    parser.add_argument(
        '-o',
        '--excelname',
        default='GFSM_data.xlsx',
        help='Name of the output excel sheet. Defaults to GFSM_data.xlsx',
    )
    parser.add_argument(
        '-n',
        '--noexcel',
        action='store_false',
        help='Suppress creation of the excel sheet',
    )
    args = parser.parse_args()
    
    basedir = Path(args.basedir)
    filepath = basedir / args.egradname
    
    parsed_data = extract_gfsm_dipole_data(
        filepath,
    )
    
    
    if len(args.states) < 2:
        print(
            'the --states option requires at least two states. Note that'
            + 'the ground state must be stated explicitly, i.e. --states 0 1'
        )
        sys.exit()
    
    delta_gfsm, sigma_gfsm, term_table = gfsm(
        parsed_data,
        states = [int(x) for x in args.states]
    )
    
    dipole_tab = pd.DataFrame.from_dict(
        parsed_data['dipoles'],
        columns = [
            'x /a.u.',
            'y /a.u.',
            'z /a.u.',
        ],
        orient = 'index',
    )
    dipole_tab['norm /a.u.'] = np.sqrt(
        (dipole_tab['x /a.u.'])**2
        + (dipole_tab['y /a.u.'])**2
        + (dipole_tab['z /a.u.'])**2
    )
    
    excitation_tab = pd.DataFrame(
        parsed_data['excitations'],
        columns = ['Excitation energy /a.u.']
    )
    
    summary_tab = pd.DataFrame(
        [
            len(args.states),
            ' '.join(args.states),
            args.states[-1],
            parsed_data['excitations'][int(args.states[-1])],
            delta_gfsm,
            sigma_gfsm,
        ],
        index = [
            'Number of states in model',
            'States used',
            'State of interest',
            'Excitation energy /a.u.',
            'delta_TPA /a.u.',
            'sigma_TPA /a.u.',
        ]
    )
    
    with pd.ExcelWriter(args.excelname, engine='openpyxl') as writer:
        summary_tab.to_excel(writer, sheet_name='Summary')
        term_table.to_excel(writer, sheet_name='GFSM terms')
        dipole_tab.to_excel(writer, sheet_name='Dipole moments')
        excitation_tab.to_excel(writer, sheet_name='Excitation energies')
    
    print(f'Excitation Energy /a.u.: {parsed_data['excitations'][int(args.states[-1])]:.6f}')
    print(f'delta_TPA /a.u.: {delta_gfsm:.2f}')
    print(f'sigma_TPA /a.u.: {sigma_gfsm:.2f}')
