from tpatools.state import State
from pathlib import Path
import numpy as np
import pandas as pd
import re


def is_int(value):
    try:
        x = int(value)
        return True
    except:
        return False

def parse_escf(filepath, nontpa=False):
    filepath = Path(filepath)
    states = []
    homo = None
    with filepath.open() as log:
        reading_state = None
        state_index = 0
        reading_dipole = False
        reading_strength = None
        reading_contributions = False
        logfile = log.readlines()
        for i, line in enumerate(logfile):
            if 'number of occupied orbitals' in line:
                homo = int(line.split()[-1])
            if 'excitation' in line and is_int(line.split()[0]):
                if reading_state is not None:
                    states.append(reading_state)
                reading_state = State(line.strip(), homo=homo)
            if 'Excitation energy:' in line and reading_state is not None:
                excitation_energy = float(line.split()[2])
                reading_state.set_excitation_energy(excitation_energy)

            if 'length representation' in line and reading_state is not None:
                if reading_state.oscillator_strength is None:
                    reading_state.set_osc(float(line.split()[2]))

            if reading_contributions:
                #print(line)
                try:
                    if is_int(line.split()[0]):
                        reading_state.add_contribution(line)
                except IndexError:
                    blanklinecount += 1
                    if blanklinecount >= 2:
                        reading_contributions = False


            if 'Dominant contributions:' in line:
                blanklinecount = 0
                reading_contributions = True

            if reading_dipole: 
                if 'x ' in line:
                    current_dipole['x'] = float(line.split()[1])
                if 'y ' in line:
                    current_dipole['y'] = float(line.split()[1])
                if 'Norm / debye:' in line and 'z' in line:
                    current_dipole['z'] = float(line.split()[1])
                    current_dipole['norm'] = float(line.split()[5])
                    reading_state.set_transition_dipole(**current_dipole)
                    reading_dipole = False


            if 'Electric transition dipole moment (length rep.)' in line:
                reading_dipole = True
                current_dipole = {}

            if 'Two-photon absorption amplitudes for transition to the' in line:
                #print(reading_state.name)
                #print(states)
                if reading_state is not None and reading_state not in states:
                    states.append(reading_state)
                    #print('TEST')
                #index = int(re.sub('\\D', '', line.split()[7])) - 1
                #print(index)
                reading_state = states[state_index]
                state_index += 1

            if 'transition strength [a.u.]:' in line:
                #print(line)
                reading_state.set_strength(float(line.split()[3]))

            if 'all done' in line and reading_state is not None and reading_state not in states:
                states.append(reading_state)

    return states

def escf_table(filepath):
    states = parse_escf(filepath)
    values = {
        'State' : [x.number for x in states],
        'Excitation Energy /eV' : [x.excitation_energy * 27.2114 for x in states],
        '2PA Strength /a.u.' : [x.transition_strength for x in states],
        'Cross Section /GM' : [x.get_cross_section() for x in states],
    }
    df = pd.DataFrame(values)
    return df

#   def escf_table(filepath):
#       filepath = Path(filepath)
#       with filepath.open() as log:
#           reading_state = None
#           reading_strength = None
#           states = []
#           for line in log.readlines():
#               if 'number of occupied orbitals' in line:
#                   homo = int(line.split()[-1])
#               if 'excitation' in line and is_int(line.split()[0]):
#                   stateno = line.strip()
#                   reading_state = State(stateno, homo=homo)
#               if 'Excitation energy:' in line and reading_state is not None:
#                   excitation_energy = float(line.split()[2])
#                   reading_state.set_excitation_energy(excitation_energy)
#                   states.append(reading_state)
#               if 'Two-photon absorption amplitudes for transition to the ' in line:
#                   index = int(re.sub('\\D', '', line.split()[7])) - 1
#                   reading_state = states[index]
#   
#               if 'transition strength [a.u.]:' in line:
#                   reading_state.set_strength(float(line.split()[3]))
#      
#           values = {
#               'State' : [x.number for x in states],
#               'Excitation Energy /eV' : [x.excitation_energy * 27.2114 for x in states],
#               '2PA Strength /a.u.' : [x.transition_strength for x in states],
#               'Cross Section /GM' : [x.get_cross_section() for x in states],
#           }
#           df = pd.DataFrame(values)
#           return df



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
            if is_int(line.split()[0]):
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



def parse_egrad(filepath):
    """
        Obtain permanent and transition dipole moments from a TPA egrad calculation
    """
    data = {}
    au_debye_convfactor = 2.54158
    filepath = Path(filepath)
    with filepath.open() as logfile:
        log = [x.strip() for x in logfile.readlines() if x.strip() != '']

        grndex = log.index('Electric dipole moment:')
        data['mu_00'] = {
            'x' : round(float(log[grndex + 2].split()[3]) * au_debye_convfactor, 6),
            'y' : round(float(log[grndex + 3].split()[3]) * au_debye_convfactor, 6),
            'z' : round(float(log[grndex + 4].split()[3]) * au_debye_convfactor, 6),
            'norm' : float(log[grndex + 4].split()[-1]),
        }

        # Determine which excited state has been selected by $exopt, then grab the relevant transition dipole moment
        transition_indices = [x for x, line in enumerate(log) if line == 'Electric transition dipole moment (length rep.):']
        exstateno = [int(x.split()[3]) for x in log if re.search('Excited state no\.    .+ chosen for optimization', x)][0]
        tndex = transition_indices[exstateno - 1]
        data['mu_01'] = {
            'x' : round(float(log[tndex + 1].split()[1]) * au_debye_convfactor, 6),
            'y' : round(float(log[tndex + 2].split()[1]) * au_debye_convfactor, 6),
            'z' : round(float(log[tndex + 3].split()[1]) * au_debye_convfactor, 6),
            'norm' : float(log[tndex + 3].split()[-1]),
        }

        pmdex = log.index('dipole moment')
        data['mu_11'] = {
            'x' : round(float(log[pmdex + 2].split()[-1]) * au_debye_convfactor, 6),
            'y' : round(float(log[pmdex + 3].split()[-1]) * au_debye_convfactor, 6),
            'z' : round(float(log[pmdex + 4].split()[-1]) * au_debye_convfactor, 6),
            'norm' : float(log[pmdex + 5].split()[8]),
        }

    return data
        


def get_ground_state_dipole(logfile):
    """
        Grab the ground state dipole moment from an escf output file
    """
    with logfile.open() as log:
        outlines = [x.strip() for x in log.readlines()]
        grindex = outlines.index('Ground state')
        for line in outlines[grindex:]:
            if 'Norm / debye' in line:
                try:
                    return float(line.split()[7])
                except IndexError:
                    print(f"""
                    hmm...
                    there should be a dipole moment here...
                    but something is WRONG

                    Check your output file ({logfile.resolve()})
                    or maybe write a better parser next time, bozo
                    """)

    print(f'something wrong, check outfile ({logfile.resolve()})')
    return 'NA'

def gather_state_data(
        basedir, 
        outfilename = 'escf.out', 
        egradoutname = 'egrad.out',
        state=1, 
        egradavail=False,
        orderedkeys=None,
        fulldirnames = False,
        suppress_egrad_warning = False,
        tabulate = False,
        latexnames = False,
    ):
    """
        Recursively gather excitation energies, transition dipoles, cross sections, and dipole moments for all output files in a given directory and compile them into a dictionary, with keys provided by the directory names
    """

    basedir = Path(basedir)
    filelist = list(basedir.rglob(outfilename))

    if orderedkeys is not None:
        filelist = sorted(
            filelist,
            key=lambda x: orderedkeys.index(x.parent.name)
        )
    elif fulldirnames:
        filelist = sorted(filelist, key=lambda x: str(x.relative_to(basedir).parent))
    else:
        filelist = sorted(filelist, key=lambda x: str(x.parent.name))


    collectdata = {}
    for logfile in filelist:
        try:
            statedata = parse_escf(logfile)[state - 1]
        except IndexError:
            print(f'Excited state {state} not available for file {logfile.resolve()}')
            continue

        try:
            dipdata = parse_egrad(logfile.parent / egradoutname)
            mu_00 = dipdata['mu_00']['norm']
            mu_01 = dipdata['mu_01']['norm']
            mu_11 = dipdata['mu_11']['norm']
        except FileNotFoundError:
            if suppress_egrad_warning == False:
                print(f'No egrad outfile in directory {logfile.parent}')
            mu_00 = get_ground_state_dipole(logfile)
            mu_01 = statedata.transition_dipole['norm']
            mu_11 = 'NA'

        if fulldirnames:
            key = str(logfile.relative_to(basedir).parent)
        else:
            key = logfile.parent.name

        try:
            if latexnames:
                collectdata[key] = {
                    '$\\Delta E$ /eV' : statedata.excitation_energy * 27.2114,
                    '$\delta^{\\textrm{2PA}}$ /a.u.' : statedata.transition_strength,
                    '$\\sigma^{\\textrm{2PA}}$ /GM' : statedata.get_cross_section(),
                    '$|\\mu_{00}|$ /D' : mu_00,
                    '$|\\mu_{01}|$ /D' : mu_01,
                    '$|\\mu_{11}|$ /D' : mu_11,
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


