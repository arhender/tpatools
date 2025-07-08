## LEARN REGULAR EXPRESSIONS FOR REWRITE
import re
import numpy as np
import pandas as pd
from pathlib import Path
from tpatools.state import State
import itertools

def _is_int(value):
    try:
        x = int(value)
        return True
    except:
        return False

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
            x\s+(?P<mu_00_x>-?\d+\.\d+(?:[eE][+-]?\d+)?)\s+
            y\s+(?P<mu_00_y>-?\d+\.\d+(?:[eE][+-]?\d+)?)\s+
            z\s+(?P<mu_00_z>-?\d+\.\d+(?:[eE][+-]?\d+)?)\s+
            \|\sdipole\smoment\s\|.*=\s+
            (?P<mu_00_norm>-?\d+\.\d+(?:[eE][+-]?\d+)?)
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
            (?P<excitation_energy>\d+\.\d+(?:[eE][+-]?\d+)?)
            \s+a\.u.*?
            oscillator\sstrength\s\(length\sgauge\)
            \s+:\s+(?P<oscillator_strength>\d+\.\d+(?:[eE][+-]?\d+)?)
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
            (?P<excitation_energy>\d+\.\d+(?:[eE][+-]?\d+)?)
            \sa\.u.*?
            1ST\sPHOTON:\s+
            (?P<photon_1>\d+\.\d+(?:[eE][+-]?\d+)?)\sa\.u.*?
            2ND\sPHOTON:\s+
            (?P<photon_2>\d+\.\d+(?:[eE][+-]?\d+)?)\sa\.u.*?
            ROTATIONALLY\sAVERAGED\sVALUES:.*?
            Linear:\s+
            (?P<tpa_strength>-?\d+\.\d+(?:[eE][+-]?\d+)?)
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
            (?P<excitation_energy>\d+\.\d+(?:[eE][+-]?\d+)?).*?
            Analysis\sof\srelaxed\sproperties:.*?
            dipole\smoment:.*?
            x\s+(?P<mu_11_x>-?\d+\.\d+(?:[eE][+-]?\d+)?)\s+
            y\s+(?P<mu_11_y>-?\d+\.\d+(?:[eE][+-]?\d+)?)\s+
            z\s+(?P<mu_11_z>-?\d+\.\d+(?:[eE][+-]?\d+)?)\s+
            \|\sdipole\smoment\s\|.*?=\s+
            (?P<mu_11_norm>-?\d+\.\d+(?:[eE][+-]?\d+)?)
            \s+debye
            .*?Analysis\sof\sunrelaxed\sproperties
        """,
        re.DOTALL | re.X
    )

    data = {
        'mu_00' : {},
    }
    try:
        mu_00dat = pattern_gs.search(log).groupdict()
        # Norm is parsed in debye: add to data
        data['mu_00']['norm'] = float(mu_00dat['mu_00_norm'])
    
        # convert x,y,z components of mu_00 to debye
        for vect in ['x', 'y', 'z']:
            data['mu_00'][vect] = float(mu_00dat[f'mu_00_{vect}']) * au_debye_convfactor
    except AttributeError:
        print('ground state dipole not parsed')
    
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

def _convert_dipole_axes(
        statedat,
        stringstart,
        stringend='',
        as_dict = True,
    ):
    au_debye_convfactor = 2.541746472
    if as_dict:
        converted = {}
        for cart in ['x', 'y', 'z']:
            converted[cart] = float(statedat[f'{stringstart}{cart}{stringend}']) * au_debye_convfactor
    else:
        converted = []
        for cart in ['x', 'y', 'z']:
            converted.append(
                float(statedat[f'{stringstart}{cart}{stringend}']) * au_debye_convfactor
            )
    return converted


def parse_egrad(
        filepath,
    ):
    filepath = Path(filepath)
    with filepath.open() as logfile:
        log = logfile.read()

    pattern_egrad = re.compile(
        r""" 
        Ground\sstate.*?
        Electric\sdipole\smoment:.*?
        x.*?(?P<mu_00_x>-?\d+\.\d+(?:[eE][+-]?\d+)?)
        \s+Norm.*?
        y\s+[\d.-]+\s+[\d.-]+\s+
        (?P<mu_00_y>-?\d+\.\d+(?:[eE][+-]?\d+)?).*?
        (?P<mu_00_z>-?\d+\.\d+(?:[eE][+-]?\d+)?)
        \s+Norm\s/\sdebye:\s+
        (?P<mu_00_norm>-?\d+\.\d+(?:[eE][+-]?\d+)?).*?
        Excited\sstate\sno\.\s+
        (?P<exno>\d+)\s+chosen\s
        for\soptimization.*?
        electrostatic\smoments.*?
        dipole\smoment.*?
        x.[\s\d.-]+\s+(?P<mu_11_x>-?\d+\.\d+(?:[eE][+-]?\d+)?)\s+
        y.[\s\d.-]+\s+(?P<mu_11_y>-?\d+\.\d+(?:[eE][+-]?\d+)?)\s+
        z.[\s\d.-]+\s+(?P<mu_11_z>-?\d+\.\d+(?:[eE][+-]?\d+)?)\s+
        \|\sdipole\smoment\s\|\s=.*?
        (?P<mu_11_norm>-?\d+\.\d+)\sdebye
        """,
        re.DOTALL | re. X
    )
    scan_excitations = re.compile(
        r""" 
        (?P<stateno>\d+)\s
        (?P<multiplicity>[a-z]+)\s
        (?P<irrep>\S+)\s
        excitation.*?
        Electric\stransition\sdipole\s
        moment\s\(length\srep\.\):\s+
        x\s+(?P<mu_01_x>-?\d+\.\d+(?:[eE][+-]?\d+)?).*?
        y\s+(?P<mu_01_y>-?\d+\.\d+(?:[eE][+-]?\d+)?)\s+
        z\s+(?P<mu_01_z>-?\d+\.\d+(?:[eE][+-]?\d+)?)\s+
        Norm\s/\sdebye:\s+
        (?P<mu_01_norm>-?\d+\.\d+(?:[eE][+-]?\d+)?)
        """,
        re.DOTALL | re.X
    )
    egradata = pattern_egrad.search(log).groupdict()
    excitation_list = [x.groupdict() for x in scan_excitations.finditer(log)]
    
    chosen_state = int(egradata['exno']) - 1

    data = {
        'chosen_state' : int(egradata['exno']),
        'stateno' : int(excitation_list[chosen_state]['stateno']),
        'irrep' : excitation_list[chosen_state]['irrep'],
        'mu_00' : {
            'norm' : float(egradata['mu_00_norm']),
            **_convert_dipole_axes(
                egradata,
                'mu_00_',
            ),
        },
        'mu_01' : {
            'norm' : float(excitation_list[chosen_state]['mu_01_norm']),
            **_convert_dipole_axes(
                excitation_list[chosen_state],
                'mu_01_',
            ),
        },
        'mu_11' : {
            'norm' : float(egradata['mu_11_norm']),
            **_convert_dipole_axes(
                egradata,
                'mu_11_',
            ),
        },
    }
    return data

    



        
def parse_escf(
        filepath,
        search_for_egrad = True,
        states_only = False, # for backwards compatibility to my old version, which just returned a list of State objects
        egradname = 'egrad.out',
        suppress_egrad_notification = False,
    ):
    filepath = Path(filepath)
    with filepath.open() as logfile:
        log = logfile.read()

    if search_for_egrad and (filepath.parent / egradname).is_file():
        egradavail = True
        egradpath = filepath.parent / egradname
        if suppress_egrad_notification == False:
            print('egrad output found in escf directory! Will parse!')
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
            x.*?(?P<mu_00_x>-?\d+\.\d+(?:[eE][+-]?\d+)?)
            \s+Norm.*?
            y\s+[\d.-]+\s+[\d.-]+\s+
            (?P<mu_00_y>-?\d+\.\d+(?:[eE][+-]?\d+)?).*?
            (?P<mu_00_z>-?\d+\.\d+(?:[eE][+-]?\d+)?)
            \s+Norm\s/\sdebye:\s+
            (?P<mu_00_norm>-?\d+\.\d+(?:[eE][+-]?\d+)?)
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
            (?P<excitation_energy>-?\d+.\d+(?:[eE][+-]?\d+)?).\s+?
            Excitation\senergy\s/\seV.*?
            Oscillator\sstrength:.*?
            length\srepresentation:\s+
            (?P<oscillator_strength>\d+\.\d+(?:[eE][+-]?\d+)?).*?
            Dominant\scontributions:.*?
            (?P<mo_contributions>(?:\d+\s+\S+\s+-?[\d.]+\s+\d+\s+\S+\s+-?[\d.]+\s+[\d*.]*\s+)+)
            .*?Electric\stransition\sdipole\smoment\s\(length\srep\.\):\s+
            x\s+(?P<mu_01_x>-?\d+\.\d+(?:[eE][+-]?\d+)?).*?
            y\s+(?P<mu_01_y>-?\d+\.\d+(?:[eE][+-]?\d+)?).*?
            z\s+(?P<mu_01_z>-?\d+\.\d+(?:[eE][+-]?\d+)?)\s+
            Norm\s/\sdebye:\s+
            (?P<mu_01_norm>-?\d+\.\d+(?:[eE][+-]?\d+)?)
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
            (?P<excitation_energy>-?\d\.\d+(?:[eE][+-]?\d+)?)
            \s+Hartree.*?
            omega_1\s+
            (?P<photon_1>[+-]?[\d.]+(?:[eE][+-]?\d+)?)
            .*?omega_2\s+
            (?P<photon_2>[+-]?[\d.]+(?:[eE][+-]?\d+)?)
            .*?transition\sstrength\s\[a\.u\.\]:\s+
        (?P<tpa_strength>-?\d+\.\d+(?:[eE][-+]?\d+)?)
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

    if egradavail:
        egrad_data = parse_egrad(egradpath)
        egrad_state = egrad_data['stateno']
        egrad_irrep = egrad_data['irrep']
        matching_states = [x for x in states if x.number == egrad_state and x.irrep == egrad_irrep]
        if len(matching_states) > 1:
            print('ERROR: Unexpected double state match for egrad, something is funky (I coded it wrong, clearly)')
        match_index = states.index(matching_states[0])
        states[match_index].set_permanent_dipole(
            **egrad_data['mu_11']
        )

    if states_only:
        return states
    else:
        data['states'] = states
        return data


        
def parse_results(
        filepath,
        egradoutname='egrad.out', ## for escf parsing, this allows nonstandard checks for the egrad name
        search_for_egrad = True, # passed to escf parser is escf file supplied
        suppress_egrad_notification = False,
    ):
    """
        General parser for escf, egrad, and ricc2 OPA and TPA calculations
    """
    filepath = Path(filepath)
    with filepath.open() as logfile:
        log = logfile.read()

    if re.search('e s c f', log):
        return parse_escf(
            filepath, 
            egradname=egradoutname, 
            search_for_egrad=search_for_egrad,
            suppress_egrad_notification=suppress_egrad_notification,
        )

    elif re.search('R I C C 2', log):
        return parse_ricc2(filepath)
        
    elif re.search('e g r a d', log):
        return parse_egrad(filepath)

    else:
        print(f'ERROR: Output file {filepath.name} unrecognized, does not seem to correspond to any of the accepted calculations (escf, ricc2, or egrad).')
        return None

def _get_state_label(state, irrep=False, mult=False):
    if irrep and mult:
        mult_letter = state.mult[0].capitalize()
        state_id = f'{mult_letter}{state.number}{state.irrep}'
    elif irrep:
        state_id = f'{state.number}{state.irrep}'
    elif mult:
        mult_letter = state.mult[0].capitalize()
        state_id = f'{mult_letter}{state.number}'
    else:
        state_id = state.number

    return state_id


def tpa_table(
        filepath,
        irrep = True,
        mult = False,
    ):
    states = parse_results(filepath, search_for_egrad=False)['states']

    values = {
        'State' : [_get_state_label(x, irrep=irrep, mult=mult) for x in states],
        'Excitation Energy /eV' : [x.excitation_energy * 27.2114 for x in states],
        '2PA Strength /a.u.' : [x.transition_strength for x in states],
        'Cross Section /GM' : [x.get_cross_section() for x in states],
    }
    df = pd.DataFrame(values)
    return df

        
def parse_exspec(filepath):
    """
        For the expspectrum file which provides vertical excitations (UV-Vis)
    """
    filepath = Path(filepath)
    excno = []
    hartree = []
    eV = []
    ecm = []
    enm = []
    oscvel = []
    osclen = []
    with filepath.open() as log:
        for line in log.readlines():
            if _is_int(line.split()[0]):
                excno.append(int(line.split()[0]))
                hartree.append(float(line.split()[2]))
                eV.append(float(line.split()[3]))
                ecm.append(
                    float(
                        line.split()[4].replace('D','E')
                    )
                )
                enm.append(float(line.split()[5]))
                oscvel.append(float(line.split()[6]))
                osclen.append(float(line.split()[7]))

    return {
        'Excitation numbers' : excno,
        'E /hartree' : hartree,
        'E /eV' : eV,
        'E /cm$^{-1}$' : ecm,
        'E /nm' : enm,
        'Oscillator (vel)' : oscvel,
        'Oscillator (len)' : osclen,
    }


def gather_state_data(
        basedir, 
        outfilename = None, 
        egradoutname = 'egrad.out',
        state=1, 
        egradavail=False,
        orderedkeys=None,
        fulldirnames = False,
        suppress_egrad_notification = False,
        tabulate = False,
        latexnames = False,
        compactnames = False,
        verbose_output = False,
        osc = False,
        sortfunc = None,
    ):
    """
        Recursively gather excitation energies, transition dipoles, cross sections, and dipole moments for all output files in a given directory and compile them into a dictionary, with keys provided by the directory names
    """

    basedir = Path(basedir)
    if outfilename is None:
        # Search for any files in the fallback list
        fallbacknames = ['escf.out', 'bse.out', 'tpa.out', 'td-dft.out', 'ricc2.out']
        filelist = []
        for name in fallbacknames:
            filelist.extend(basedir.rglob(name))
    else:
        filelist = list(basedir.rglob(outfilename))

    if orderedkeys is not None:
        filelist = sorted(
            filelist,
            key=lambda x: orderedkeys.index(x.parent.name)
        )
    elif sortfunc is not None:
        filelist = sorted(
            filelist,
            key=sortfunc,
        )
    elif fulldirnames:
        filelist = sorted(filelist, key=lambda x: str(x.relative_to(basedir).parent))
    else:
        filelist = sorted(filelist, key=lambda x: str(x.parent.name))


    collectdata = {}
    for logfile in filelist:
        if verbose_output:
            print(f'Analyzing log file {logfile.resolve()}')
        try:
            parsed_data = parse_results(logfile, egradoutname=egradoutname, suppress_egrad_notification=suppress_egrad_notification)
            statedata = parsed_data['states'][state - 1]
        except IndexError:
            print(f'Excited state {state} not available for file {logfile.resolve()}')
            print(parsed_data)
            continue
        except:
            print(f'Error parsing output file {logfile.resolve()}, check for issues')
            continue


        #try:
        #    dipdata = parse_egrad(logfile.parent / egradoutname)
        #    mu_00 = dipdata['mu_00']['norm']
        #    mu_01 = dipdata['mu_01']['norm']
        #    mu_11 = dipdata['mu_11']['norm']
        #except FileNotFoundError:
        #    if suppress_egrad_warning == False:
        #        print(f'No egrad outfile in directory {logfile.parent}')
        #    mu_00 = get_ground_state_dipole(logfile)
        #    mu_01 = statedata.transition_dipole['norm']
        #    mu_11 = 'NA'

        mu_00 = parsed_data['mu_00']['norm']
        mu_01 = statedata.transition_dipole['norm']
        if statedata.permanent_dipole == {}:
            mu_11 = 'NA'
        else:
            mu_11 = statedata.permanent_dipole['norm']

        if fulldirnames:
            key = str(logfile.relative_to(basedir).parent)
        else:
            key = logfile.parent.name

        try:
            if osc:
                if latexnames:
                    collectdata[key] = {
                        '$\\Delta E$ /eV' : statedata.excitation_energy * 27.2114,
                        '$\delta^{\\textrm{2PA}}$ /a.u.' : statedata.transition_strength,
                        '$\\sigma^{\\textrm{2PA}}$ /GM' : statedata.get_cross_section(),
                        '$|\\mu_{00}|$ /D' : mu_00,
                        '$|\\mu_{01}|$ /D' : mu_01,
                        '$|\\mu_{11}|$ /D' : mu_11,
                        '$f$' : statedata.oscillator_strength,
                    }
                elif compactnames:
                    collectdata[key] = {
                        'ex. E' : statedata.excitation_energy * 27.2114,
                        'delta /a.u.' : statedata.transition_strength,
                        'sigma /GM' : statedata.get_cross_section(),
                        'mu_00' : mu_00,
                        'mu_01' : mu_01,
                        'mu_11' : mu_11,
                        'osc' : statedata.oscillator_strength,
                    }

                else:
                    collectdata[key] = {
                        'Excitation Energy /eV' : statedata.excitation_energy * 27.2114,
                        '2PA Strength /a.u.' : statedata.transition_strength,
                        'Cross Section /GM' : statedata.get_cross_section(),
                        'Ground State Dipole Moment /D' : mu_00,
                        'Transition Dipole Moment /D' : mu_01,
                        'Excited State Dipole Moment /D' : mu_11,
                        'Oscillator Strength' : statedata.oscillator_strength,
                    }
            else:
                if latexnames:
                    collectdata[key] = {
                        '$\\Delta E$ /eV' : statedata.excitation_energy * 27.2114,
                        '$\delta^{\\textrm{2PA}}$ /a.u.' : statedata.transition_strength,
                        '$\\sigma^{\\textrm{2PA}}$ /GM' : statedata.get_cross_section(),
                        '$|\\mu_{00}|$ /D' : mu_00,
                        '$|\\mu_{01}|$ /D' : mu_01,
                        '$|\\mu_{11}|$ /D' : mu_11,
                    }
                elif compactnames:
                    collectdata[key] = {
                        'ex. E' : statedata.excitation_energy * 27.2114,
                        'delta /a.u.' : statedata.transition_strength,
                        'sigma /GM' : statedata.get_cross_section(),
                        'mu_00' : mu_00,
                        'mu_01' : mu_01,
                        'mu_11' : mu_11,
                    }

                else:
                    collectdata[key] = {
                        'Excitation Energy /eV' : statedata.excitation_energy * 27.2114,
                        '2PA Strength /a.u.' : statedata.transition_strength,
                        'Cross Section /GM' : statedata.get_cross_section(),
                        'Ground State Dipole Moment /D' : mu_00,
                        'Transition Dipole Moment /D' : mu_01,
                        'Excited State Dipole Moment /D' : mu_11,
                    }

        except:
            print(f'Errors in extracting data in subdir {logfile.parent}, skipping for now')
            continue

    if tabulate:
        return pd.DataFrame.from_dict(
            collectdata,
            orient='index',
        )
    else:
        return collectdata


