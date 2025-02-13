
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

        self.homo = int(homo)
        self.lumo = homo + 1

        self.excitation_energy = None
        self.photon_energies = [None, None]
        self.transition_strength = None

        self.dominant_contributions = []
        self.fmo_contributions = []

    def set_strength(self, value):
        self.transition_strength = value

    def set_excitation_energy(self, excitation_energy):
        self.excitation_energy = excitation_energy
        self.photon_energies = [excitation_energy/2, excitation_energy/2]

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
