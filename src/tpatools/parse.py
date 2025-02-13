from state import State
import numpy as np
import pandas as pd

def is_int(value):
    try:
        x = int(value)
        return True
    except:
        return False

def parse_escf(filepath):
    filepath = Path(filepath)
    states = []
    homo = None
    with filepath.open() as log:
        reading_state = None
        state_index = 0
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


            if 'Two-photon absorption amplitudes for transition to the ' in line:
                if reading_state is not None  and reading_state not in states:
                    states.append(reading_state)
                #index = int(re.sub('\\D', '', line.split()[7])) - 1
                #print(index)
                reading_state = states[state_index]
                state_index += 1

            if 'transition strength [a.u.]:' in line:
                #print(line)
                reading_state.set_strength(float(line.split()[3]))


    return states

def escf_table(filepath):
    filepath = Path(filepath)
    with filepath.open() as log:
        reading_state = None
        reading_strength = None
        states = []
        for line in log.readlines():
            if 'number of occupied orbitals' in line:
                homo = int(line.split()[-1])
            if 'excitation' in line and is_int(line.split()[0]):
                stateno = line.strip()
                reading_state = State(stateno, homo=homo)
            if 'Excitation energy:' in line and reading_state is not None:
                excitation_energy = float(line.split()[2])
                reading_state.set_excitation_energy(excitation_energy)
                states.append(reading_state)
            if 'Two-photon absorption amplitudes for transition to the ' in line:
                index = int(re.sub('\\D', '', line.split()[7])) - 1
                reading_state = states[index]

            if 'transition strength [a.u.]:' in line:
                reading_state.set_strength(float(line.split()[3]))
   
        values = {
            'State' : [x.number for x in states],
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


