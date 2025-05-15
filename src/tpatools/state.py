import numpy as np

class State():
    """ a Class to record the data for a given excitation state"""

    instances = []
    fine_structure = 0.007297
    bohr_radius = 5.29E-9 # in cm
    light_speed = 2.9979E10 # in cm/s

    def __init__(self, name, homo=None):
        self.name = name
        State.instances.append(self)

        pieces = name.split()
        self.number = int(pieces[0])
        self.mult = pieces[1]
        self.irrep = pieces[2]

        if homo is not None:
            self.homo = int(homo)
            self.lumo = homo + 1

        self.excitation_energy = None
        self.photon_energies = [None, None]
        self.transition_strength = None

        self.oscillator_strength = None

        self.dominant_contributions = []
        self.fmo_contributions = []
        self.transition_dipole = {}

        self.permanent_dipole = {}

    def set_strength(self, value):
        self.transition_strength = value

    def set_osc(self, value):
        self.oscillator_strength = value

    def set_excitation_energy(self, excitation_energy, degenerateTPA=True):
        self.excitation_energy = excitation_energy
        if degenerateTPA:
            self.photon_energies = [excitation_energy/2, excitation_energy/2]

    def set_photon_energies(self, photon_energies):
        self.photon_energies = photon_energies

    def set_transition_dipole(self, norm, x, y, z):
        self.transition_dipole = {
            'norm' : norm,
            'x' : x,
            'y' : y,
            'z' : z,
        }
    def set_permanent_dipole(self, norm, x, y, z):
        self.permanent_dipole = {
            'norm' : norm,
            'x' : x,
            'y' : y,
            'z' : z,
        }


    def add_contribution(self, line):
        occ = int(line.split()[0])
        virt = int(line.split()[3])
        symm = {
            'occ' : line.split()[1],
            'virt' : line.split()[4],
        }
        energies = {
            'occ' : line.split()[2],
            'virt' : line.split()[5],
        }
        coeff = float(line.split()[6])

        self.dominant_contributions.append({
            'occ' : occ,
            'virt' : virt,
            'coeff' : coeff,
            'energies' : energies,
            'symm' : symm,
        })
        
        if occ == self.homo:
            occ_id = 'HOMO'
        elif occ < self.homo:
            occ_id = f'HOMO - {self.homo - occ}'
        else:
            occ_id = f'HOMO + {occ - self.homo}'
        
        if virt == self.lumo:
            virt_id = 'LUMO'
        elif virt > self.lumo:
            virt_id = f'LUMO + {virt - self.lumo}'
        else:
            virt_id = f'LUMO - {self.lumo - virt}'

        self.fmo_contributions.append({
            'occ' : occ_id,
            'virt' : virt_id,
            'coeff' : coeff,
            'energies' : energies,
            'symm' : symm,
        })


    def get_contributions(self, fmo_relative = True):
        if fmo_relative:
            return self.fmo_contributions
        else: 
            return self.dominant_contributions

    def get_coeff(self, occ, virt, fmo_relative = True):
        if fmo_relative:
            targdat = self.fmo_contributions
        else:
            targdat = self.dominant_contributions
        for cont in targdat:
            if cont['occ'] == occ and cont['virt'] == virt:
                return cont['coeff']
        
        return 0


    def get_cross_section(self):
        cross_section = (
            4 * (np.pi**2) * (State.bohr_radius**5) * State.fine_structure 
            * self.photon_energies[0] * self.photon_energies[1] 
            * self.transition_strength
            / (State.light_speed * 0.1 / 27.2)
        )
        cross_section = cross_section * 10**50
        return cross_section

    def as_dict(self):
        if self.permanent_dipole != {}:
            return {
                'Excitation Energy /eV' : self.excitation_energy * 27.2114,
                '2PA Strength /a.u.' : self.transition_strength,
                '2PA Cross Section /GM' : self.get_cross_section(),
                'Oscillator Strength /a.u.' : self.oscillator_strength,
                'Transition Dipole Moment /D' : self.transition_dipole,
                'Permanent Dipole Moment /D' : self.permanent_dipole,
            }
        else:
            return {
                'Excitation Energy /eV' : self.excitation_energy * 27.2114,
                '2PA Strength /a.u.' : self.transition_strength,
                '2PA Cross Section /GM' : self.get_cross_section(),
                'Oscillator Strength /a.u.' : self.oscillator_strength,
                'Transition Dipole Moment /D' : self.transition_dipole,
            }

