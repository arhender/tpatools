import numpy as np
import pandas as pd
from pathlib import Path
import re

def get_cos(
    a,
    b,
):
    """
    Returns the cosine of an angle between two dipole moments
    
    Args:
        a (list): the first dipole moment vector, as a list [x,y,z]
        b (list): the second dipole moment vector, as a list [x,y,z]
    Returns:
        cosine of the angle between a and b
    """
    norm_a = np.linalg.norm(a)
    norm_b = np.linalg.norm(b)
    dot_product = np.dot(a, b)
    cos_theta = dot_product / (norm_a * norm_b)
    return cos_theta

def _dip(k, l):
    """
    Get dipole dictionary key for dipole moment (ie mu_01 for a transition
        dipole moment to the first excited state)

    Args:
        k (int): state 1 in the dipole moment
        l (int): state 2 in the dipole moment
    Returns: 
        string with key to access the dipole as parsed in the 
            extract_gfsm_dipole_data function
    """
    return f'mu_{k}{l}'

def au_to_GM(
    au_transition_strength,
    delta_E,
    N=4,
):
    """
    Convert 2PA transition strength in a.u. to the cross section in GM

    Args:
        au_transition_strength (float): delta_2PA value in au
        delta_E (float): excitation energy to excited state of 
            interest
    Returns:
        2PA cross section in GM
    """
    fine_structure = 0.007297
    bohr_radius = 5.29E-9 # in cm
    light_speed = 2.9979E10 # in cm/s

    cross_section = (
        4 * (np.pi**2) 
        * bohr_radius**5
        * fine_structure 
        * (delta_E / 2)**2 
        / (light_speed * 0.1 / 27.2) 
    ) * au_transition_strength * 10**50
    return cross_section



def _parse_exenergy(
    logstring,
):
    """
    Parse the single excitation portions of the egrad output file

    Args:
        logstring (string): text of the output file
    Returns:
        list of excitation energies, dict of transition dipole moments
    """
    excitations = []
    dipoles = {}
    irrep_re_group = r"\w'\""
    float_number_all_re = r"[+-]?[0-9]*[.]?[0-9]+(?:[ED][+-]\d{2})?"
    r_excited_state = (
        r'(?P<stateno>[0-9]+)\s+[\w]+\s(?P<irrep>['
        + irrep_re_group + r']+)\s+excitation.*?'
        + r'Excitation energy:\s+(?P<ex_energy>'
        + float_number_all_re
        + r').*?Electric transition dipole moment'
        + r' \(length rep.\):\s+x\s+(?P<mu_x>'
        + float_number_all_re
        + r').*?y\s+(?P<mu_y>'
        + float_number_all_re
        + r')\s+z\s+(?P<mu_z>'
        + float_number_all_re + r')'
    )
    regex_excited_state = re.compile(r_excited_state, re.DOTALL)
    match = regex_excited_state.finditer(logstring)
    if match:
        for block in match:
            params = block.groupdict()
            excitations.append(
                float(params['ex_energy'])
            )
            dipname = f'mu_0{params['stateno']}'
            dipoles[dipname] = [
                float(params[f'mu_{ax}']) for ax in 'xyz'
            ]

    return excitations, dipoles


def extract_gfsm_dipole_data(
    egrad_output, 
):
    """
    extract dipoles needed for a GFSM calculation.

    Args:
        egrad_output (str or Path): path to file containing egrad results
            with $nacme block to get state-to-state transition dipoles.
    Returns:
        a dictionary containing excitation energies and dipole moments
    """
    dipoles = {}
    excitations = [0]


    # Setup regex for ground state (may be used by either condition,
    # better to set here
    float_number_all_re = r"[+-]?[0-9]*[.]?[0-9]+(?:[ED][+-]\d{2})?"
    irrep_re_group = r"\w'\""

    regex_gs_dip = re.compile(
        (
            r'Ground state.*?'
            + r'Electric dipole moment:.*?('
            + float_number_all_re
            + r')\s+Norm:.*?('
            + float_number_all_re
            + r')\s+z.*?('
            + float_number_all_re
            + r')\s+Norm / debye'
        ),
        re.DOTALL
    )

    if egrad_output is not None:
        egrad_output = Path(egrad_output)
        with egrad_output.open() as file:
            egrad_output = file.read()

        match = regex_gs_dip.search(egrad_output)
        if match:
            # add gs dipole data
            dipoles['mu_00'] = [float(x) for x in match.groups()]
        
        # regex to find state-to-state transition dipoles
        r_betweenstate_transition_moment = (
            r'<\s+([0-9]+)\|\s+W\s+\|\s+([0-9]+)>\s+transition moments'
            + r'.*?Relaxed electric transition dipole moment \(length rep.\):\s+x\s+('
            + float_number_all_re + r').*?y\s+(' + float_number_all_re + r')\s+z\s+('
            + float_number_all_re + r')'
        )
        regex_betweenstate_transition_moment = re.compile(r_betweenstate_transition_moment, re.DOTALL)
        match = regex_betweenstate_transition_moment.finditer(egrad_output)
        for block in match:
            dipname = f'mu_{block.group(1)}{block.group(2)}'
            dipoles[dipname] = [float(x) for x in block.groups()[2:]]

        parsed_ex, parsed_mu0x = _parse_exenergy(egrad_output)
        

        r_difference_moments = (
            r'<\s+([0-9]+)\|\s+W\s+\|\s+([0-9]+)>\s+'
            + r'-\s+<\s*[0-9]+\s*\|\s*W\s*\|[0-9]+>'
            + r'\s+difference moments'
            + r'.*?Relaxed electric dipole moment \(length rep.\):\s+x\s+('
            + float_number_all_re + r').*?y\s+(' + float_number_all_re + r')\s+z\s+('
            + float_number_all_re + r')'
        )
        regex_difference_moments = re.compile(r_difference_moments, re.DOTALL)
        match = regex_difference_moments.finditer(egrad_output)
        for block in match:
            dipname = f'mu_{block.group(1)}{block.group(2)}'
            mu_to_gs_diff = [float(x) for x in block.groups()[2:]]
            perm_dip = np.array(mu_to_gs_diff) + np.array(dipoles['mu_00'])
            dipoles[dipname] = [round(float(x), 6) for x in perm_dip]

        excitations.extend(parsed_ex)
        dipoles.update(parsed_mu0x)


    dipoles = dict(sorted(dipoles.items()))
    
    return {
        'dipoles' : dipoles,
        'excitations' : excitations
    }
        


def gfsm(
    data,
    states = [0, 1]
): 
    """
    Perform a calculation of the TPA cross section using the GFSM (Generalized Few State
        Model) equation

    Args:
        data: dictionary dipole moments and excitation energies, in the form output by the
            extract_gfsm_dipole_data function
        states: the states (including the ground state!!) to be considered within the GFSM
            calculation. ex [0,1] for a two-state model to the first excited state,
            [0,1,4] for a three-state model to the 4th excited state, etc.

    Returns:
        delta_tpa (calculated TPA transition moment), sigma_tpa (calculated TPA cross section)
            and term_table (pandas dataframe containing the relevant values used in calculation
            of each term within the GFSM calculation)
    """
    delta_gfsm = 0
    term_table = pd.DataFrame(
        index = [
            'K',
            'L',
            'Delta E_K',
            'Delta E_L',
            'mu_0K',
            'mu_0L',
            'mu_KJ',
            'mu_LJ',
            'dipole product',
            'cos_theta_KJ_0K',
            'cos_theta_LJ_0L',
            'cos_theta_0L_0K',
            'cos_theta_LJ_KJ',
            'cos_theta_LJ_0K',
            'cos_theta_KJ_0L',
            'angle term',
            'delta_TPA contribution /a.u.',
            'sigma_TPA contribution /GM'
        ]
    )
    j = states[-1]
    max_state = data['excitations'][j]
    dipoles = data['dipoles']
    firstpass = True
    for kno, k in enumerate(states):
        for lno, l in enumerate(states):
            gfsm_term_name = f'delta_0{j}{k}{l}'
            delta_ek = max_state * 0.5 - data['excitations'][k]
            delta_el = max_state * 0.5 - data['excitations'][l]
            #print(delta_ek)

            # Calculate the relevant cosines
            dipole_norms = [
                np.linalg.norm(dipoles[_dip(0, k)]),
                np.linalg.norm(dipoles[_dip(0, l)]),
                np.linalg.norm(dipoles[_dip(k, j)]),
                np.linalg.norm(dipoles[_dip(l, j)]),
            ]
            dipole_prod_term = np.prod(dipole_norms)

            cos_angle_vals = [
                get_cos(
                    dipoles[_dip(k, j)],
                    dipoles[_dip(0, k)],
                ),
                get_cos(
                    dipoles[_dip(l, j)],
                    dipoles[_dip(0, l)],
                ),
                get_cos(
                    dipoles[_dip(0, l)],
                    dipoles[_dip(0, k)],
                ),
                get_cos(
                    dipoles[_dip(l, j)],
                    dipoles[_dip(k, j)],
                ),
                get_cos(
                    dipoles[_dip(l, j)],
                    dipoles[_dip(0, k)],
                ),
                get_cos(
                    dipoles[_dip(k, j)],
                    dipoles[_dip(0, l)],
                ),       
            ]
            angle_term = (
                cos_angle_vals[0] * cos_angle_vals[1]
                + cos_angle_vals[2] * cos_angle_vals[3]
                + cos_angle_vals[4] * cos_angle_vals[5]
            )


            addition = (
                (4 / (15 * delta_ek * delta_el))
                * dipole_prod_term
                * angle_term
            )
            delta_gfsm += addition 

            term_data = [
                k,
                l,
                delta_ek,
                delta_el,
            ]
            term_data.extend(dipole_norms)
            term_data.append(dipole_prod_term)
            term_data.extend(cos_angle_vals)
            term_data.append(angle_term)
            term_data.append(addition)
            term_data.append(au_to_GM(addition, max_state))
            term_table[gfsm_term_name] = term_data
            

    sigma_gfsm = au_to_GM(
        delta_gfsm,
        max_state
    )

    return delta_gfsm, sigma_gfsm, term_table
