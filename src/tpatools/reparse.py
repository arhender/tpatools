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

def _state_key(statedat, overallno=None):
    #key=f'{_multno_to_text(statedat['multiplicity'])[0].capitalize()}{statedat['stateno']}/{statedat['irrep']}'
    if overallno is None:
        key = f'{statedat["stateno"]}{statedat["irrep"]}'
    else:
        key = f'{_multno_to_text(statedat["multiplicity"])[0].capitalize()}' 
        key = key + f'{overallno} {statedat["stateno"]}{statedat["irrep"]}'
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
            \sa\.u.*?
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
            model:\sCC2.*?number,\ssymmetry,\smultiplicity:\s+
            (?P<stateno>\d+)\s+
            (?P<irrep>\S+)\s+
            (?P<multiplicity>\d+).*?
            frequency\s+:\s+
            (?P<excitation_energy>\d+\.\d+).*?
            Analysis\sof\srelaxed\sproperties:.*?
            dipole\smoment:.*?
            x\s+(?P<mu_11_x>-?\d+\.\d+)\s+
            y\s+(?P<mu_11_y>-?\d+\.\d+)\s+
            z\s+(?P<mu_11_z>-?\d+\.\d+)\s+
            \|\sdipole\smoment\s\|.*?=\s+
            (?P<mu_11_norm>-?\d+\.\d+)
            \s+debye
            .*?Analysis\sof\sunrelaxed\sproperties
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

    states = []
    # following logic needed in case of for ex singlets and triplets in same calc

    for i, match in enumerate(pattern_1PA.finditer(log)):

        statedat = match.groupdict()

        states.append(
            State(_parse_ricc2_name(statedat))
        )
        states[i].set_excitation_energy(
            float(statedat['excitation_energy']),
            degenerateTPA=False,
        )
        states[i].set_osc(
            float(statedat['oscillator_strength'])
        )

        states[i].set_transition_dipole(
            _osc_to_mu_01(
                states[i].oscillator_strength,
                states[i].excitation_energy,
            ),
            *['NA', 'NA', 'NA']
        )


    for i, match in enumerate(pattern_2PA.finditer(log)):
        statedat = match.groupdict()

        states[i].set_photon_energies([
            float(statedat['photon_1']),
            float(statedat['photon_2']),
        ])
        states[i].set_strength(
            float(statedat['tpa_strength'])
        )

    for i, match in enumerate(pattern_exprop.finditer(log)):
        statedat = match.groupdict()

        print(statedat)
        permdip_axes = []
        for cart in ['x', 'y', 'z']:
            permdip_axes.append(
                float(statedat[f'mu_11_{cart}']) * au_debye_convfactor
            )
        states[i].set_permanent_dipole(
            float(statedat['mu_11_norm']),
            *permdip_axes
        )

    data['states'] = states
    return data

        
def parse_escf(
        filepath,
        search_for_egrad = True,
        states_only = False, # for backwards compatibility to my old version, which just returned a list of State objects
        egradname = 'egrad.out',
    ):
    filepath = Path(filepath)
    with filepath.open() as logfile:
        log = logfile.read()

    if search_for_egrad and (filepath.parent / egradname).is_file():
        egradavail = True
    else:
        egradavail = False

    data = {}
    au_debye_convfactor = 2.541746472

    pattern_gs = re.compile(
        r""" 
            number\sof\soccupied\s
            orbitals\s+:\s+
            (?P<homo>\d+).*?
            Ground\sstate.*?
            Electric\sdipole\smoment:.*?
            x.*?(?P<mu_00_x>-?\d+\.\d+)
            \s+Norm.*?
            y\s+[\d.-]+\s+[\d.-]+\s+
            (?P<mu_00_y>-?\d+\.\d+).*?
            (?P<mu_00_z>-?\d+\.\d+)
            \s+Norm\s/\sdebye:\s+
            (?P<mu_00_norm>-?\d+\.\d+)
        """,
        re.DOTALL | re.X
    )

    pattern_1PA = re.compile(
        r""" 
            (?P<escfname>
            (?P<stateno>\d+)\s
            (?P<multiplicity_string>\S+)\s+
            (?P<irrep>\S+)\s+
            excitation
            ).*?
            Excitation\senergy:\s+?
            (?P<excitation_energy>-?\d+.\d+).\s+?
            Excitation\senergy\s/\seV.*?
            Oscillator\sstrength:.*?
            length\srepresentation:\s+
            (?P<oscillator_strength>\d+\.\d+).*?
            Dominant\scontributions:.*?
            (?P<mo_contributions>(?:\d+\s+\S+\s+-?[\d.]+\s+\d+\s+\S+\s+-?[\d.]+\s+[\d*.]*\s+)+)
            .*?Electric\stransition\sdipole\smoment\s\(length\srep\.\):\s+
            x\s+(?P<mu_01_x>-?\d+\.\d+).*?
            y\s+(?P<mu_01_y>-?\d+\.\d+).*?
            z\s+(?P<mu_01_z>-?\d+\.\d+)\s+
            Norm\s/\sdebye:\s+
            (?P<mu_01_norm>-?\d+\.\d+)
        """,
        re.DOTALL | re.X
    )
    pattern_2PA = re.compile(
        r""" 
            Two-photon\sabsorption\samplitudes\sfor
            \stransition\sto\sthe\s+
            (?P<stateno>\d+)[a-z]{2}\s+?
            .*?in\ssymmetry\s+?
            (?P<irrep>\S+)\s+
            Exc\.\senergy:\s+
            (?P<excitation_energy>-?\d\.\d+)
            \s+Hartree.*?
            omega_1\s+
            (?P<photon_1>[+-]?[\d.]+(?:[eE][+-]?\d+)?)
            .*?omega_2\s+
            (?P<photon_2>[+-]?[\d.]+(?:[eE][+-]?\d+)?)
            .*?transition\sstrength\s\[a\.u\.\]:\s+
            (?P<tpa_strength>-?[\d.]+)
        """,
        re.DOTALL | re.X
    )

    gsdat = pattern_gs.search(log).groupdict()

    data['mu_00'] = {}
    data['mu_00']['norm'] = float(gsdat['mu_00_norm'])

    # convert x,y,z components of mu_00 to debye
    for vect in ['x', 'y', 'z']:
        data['mu_00'][vect] = float(gsdat[f'mu_00_{vect}']) * au_debye_convfactor

    data['homo'] = int(gsdat['homo'])
    data['lumo'] = int(gsdat['homo']) + 1



    states = []

    for i, match in enumerate(pattern_1PA.finditer(log)):
        statedat = match.groupdict()

        states.append(
            State(
                statedat['escfname'],
                homo = data['homo']
            )
        )
        states[i].set_excitation_energy(
            float(statedat['excitation_energy']),
            degenerateTPA=False,
        )
        states[i].set_osc(
            float(statedat['oscillator_strength'])
        )

        trandip_axes = []
        for cart in ['x', 'y', 'z']:
            trandip_axes.append(
                float(statedat[f'mu_01_{cart}']) * au_debye_convfactor
            )
        states[i].set_transition_dipole(
            float(statedat['mu_01_norm']),
            *trandip_axes
        )
        ## Add TD-DFT transition amplitudes
        for cont in [x.strip() for x in statedat['mo_contributions'].split('\n') if x.strip() != '']:
            states[i].add_contribution(cont)


    for i, match in enumerate(pattern_2PA.finditer(log)):
        statedat = match.groupdict()

        states[i].set_photon_energies([
            float(statedat['photon_1']),
            float(statedat['photon_2']),
        ])
        states[i].set_strength(
            float(statedat['tpa_strength'])
        )

    if states_only:
        return states
    else:
        data['states'] = states
        return data


        
def parse_results(filepath):
    """
        General parser for escf, egrad, and ricc2 OPA and TPA calculations
    """
    filepath = Path(filepath)
    with filepath.open() as logfile:
        log = logfile.read()

    if re.search('e s c f', log):
        mode = 'escf'
        return parse_escf(filepath)

    elif re.search('r i c c 2', log):
        mode = 'ricc2'
        return parse_ricc2(filepath)

    elif re.search('e g r a d', log):
        mode = 'egrad'
    else:
        print(f'ERROR: Output file {log} unrecognized, does not seem to correspond to any of the accepted calculations (escf, ricc2, or egrad).')
        return None




        
        
