## LEARN REGULAR EXPRESSIONS FOR REWRITE
import re
import numpy as np
import pandas as pd
from pathlib import Path
from tpatools.state import State

def _multno_to_text(nom):
    try:
        nom = int(nom)
    except TypeError:
        print('give me a number next time idiot')
    match nom:
        case 1:
            return 'singlet'
        case 2:
            return 'doublet'
        case 3:
            return 'triplet'
        case 4:
            return 'quartet'
        case 5:
            return 'quintet'
        case 6:
            return 'sextet'
        case 7:
            return 'septet'
        case 8:
            return 'octet'
        case 9:
            return 'nonet'
        case _:
            return 'jeepers how many unpaired spins are in this thing. figure it out yourself'

def _state_key(statedat):
    key=f'{_multno_to_text(statedat['multiplicity'])[0].capitalize()}{statedat['stateno']}/{statedat['irrep']}'
    return key

def _parse_ricc2_name(statedat):
    """
        Convert data into the '1 singlet a excitation' type format that the State class currently expects (based on names in escf output)
        Eventually should rework State class to remove this
    """
    name = f'{statedat['stateno']} {_multno_to_text(statedat['multiplicity'])} {statedat['irrep']} excitation'
    return name


def _osc_to_mu_01(osc_strength, e_hartree):
    return np.sqrt((3*osc_strength) / (2*e_hartree)) * 2.541746472

def parse_ricc2(filepath):
    filepath = Path(filepath)
    with filepath.open() as logfile:
        log = logfile.read()


    au_debye_convfactor = 2.541746472

    pattern_gs = re.compile(
        r""" 
            GROUND\sSTATE\sFIRST-ORDER\sPROPERTIES.*
            Analysis\sof\srelaxed\sproperties:.*
            dipole\smoment:.*
            x\s+(?P<mu_00_x>-?\d+\.\d+)\s+
            y\s+(?P<mu_00_y>-?\d+\.\d+)\s+
            z\s+(?P<mu_00_z>-?\d+\.\d+)\s+
            \|\sdipole\smoment\s\|.*=\s+
            (?P<mu_00_norm>-?\d+\.\d+)
            \s+debye
            .*Analysis\sof\sunrelaxed\sproperties.*
            ONE-PHOTON\sABSORPTION\sSTRENGTHS
        """,
        re.DOTALL | re.X
    )

    pattern_1PA = re.compile(
        r"""Transition\s+model:\sCC2.*?
            number,\ssymmetry,\smultiplicity:\s
            +(?P<stateno>\d+)\s+
            (?P<irrep>\S+)\s+
            (?P<multiplicity>\d).*?
            frequency\s+:\s+
            (?P<excitation_energy>\d+\.\d+)
            \s+a\.u.*?
            oscillator\sstrength\s\(length\sgauge\)
            \s+:\s+(?P<oscillator_strength>\d+\.\d+)
        """,
        re.DOTALL | re.X
    )

    pattern_2PA = re.compile(
        r"""
            STATE\sNO\.:\s+.*?
            (?P<stateno>\d)\s*
            SYMMETRY:\s
            (?P<irrep>\S+)\s+
            MULTIPLICITY:\s+
            (?P<multiplicity>\d+).*?
            EXCI\.\sENERGY:\s+
            (?P<excitation_energy>\d+\.\d+)
            \sa\.u.*
            1ST\sPHOTON:\s+
            (?P<photon_1>\d+\.\d+)\sa\.u.*?
            2ND\sPHOTON:\s+
            (?P<photon_2>\d+\.\d+)\sa\.u.*?
            ROTATIONALLY\sAVERAGED\sVALUES:.*?
            Linear:\s+
            (?P<tpa_strength>-?\d+\.\d+)
        """,
        re.DOTALL | re.X
    )

    pattern_exprop = re.compile(
        r"""
            Excited\sstate\sreached\sby\stransition:\s+
            model:\sCC2.*number,\ssymmetry,\smultiplicity:\s+
            (?P<stateno>\d+)\s+
            (?P<irrep>\S+)\s+
            (?P<multiplicity>\d+).*
            frequency\s+:\s+
            (?P<excitation_energy>\d+\.\d+).*
            Analysis\sof\srelaxed\sproperties:.*
            dipole\smoment:.*
            x\s+(?P<mu_11_x>-?\d+\.\d+)\s+
            y\s+(?P<mu_11_y>-?\d+\.\d+)\s+
            z\s+(?P<mu_11_z>-?\d+\.\d+)\s+
            \|\sdipole\smoment\s\|.*=\s+
            (?P<mu_11_norm>-?\d+\.\d+)
            \s+debye
            .*Analysis\sof\sunrelaxed\sproperties
        """,
        re.DOTALL | re.X
    )

    data = {
        'mu_00' : {},
    }
    mu_00dat = pattern_gs.search(log).groupdict()
    # Norm is parsed in debye: add to data
    data['mu_00']['norm'] = float(mu_00dat['mu_00_norm'])

    # convert x,y,z components of mu_00 to debye
    for vect in ['x', 'y', 'z']:
        data['mu_00'][vect] = float(mu_00dat[f'mu_00_{vect}']) * au_debye_convfactor

    states = {}

    for match in pattern_1PA.finditer(log):
        statedat = match.groupdict()
        key = _state_key(statedat)

        if key not in states:
            states[key] = State(_parse_ricc2_name(statedat))
            states[key].set_excitation_energy(
                float(statedat['excitation_energy']),
                degenerateTPA = False,
            )
            states[key].set_osc(
                float(statedat['oscillator_strength'])
            )
            states[key].set_transition_dipole(
                _osc_to_mu_01(
                    states[key].oscillator_strength,
                    states[key].excitation_energy,
                ),
                *['NA', 'NA', 'NA']
            )

    for match in pattern_2PA.finditer(log):
        statedat = match.groupdict()
        key = _state_key(statedat)
        if key not in states:
            ## not sure if this si really needed, shoudl always be present
            ## as far as i know
            states[key] = State(_parse_ricc2_name(statedat))
            states[key].set_excitation_energy(
                float(statedat['excitation_energy']),
                degenerateTPA = False,
            )
        
        states[key].set_photon_energies([
            float(statedat['photon_1']),
            float(statedat['photon_2']),
        ])

        states[key].set_strength(
            float(statedat['tpa_strength'])
        )


    for match in pattern_exprop.finditer(log):
        statedat = match.groupdict()
        key = _state_key(statedat)
        if key not in states:
            states[key] = State(_parse_ricc2_name(statedat))
            states[key].set_excitation_energy(
                float(statedat['excitation_energy']),
                degenerateTPA = False,
            )

        permdip_axes = []
        for cart in ['x', 'y', 'z']:
            permdip_axes.append(
                float(statedat[f'mu_11_{cart}']) * au_debye_convfactor
            )
        states[key].set_permanent_dipole(
            float(statedat['mu_11_norm']),
            *permdip_axes
        )

    data['states'] = states
    return data

        
        
