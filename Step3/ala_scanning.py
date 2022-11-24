'''
The idea of this step is to do an alanine scanning, for that, our code does the same that we did in the previous steps but chainging un residue at the time for alanine
so that we can see how changing a residue changes the energy of the interaction.
'''

import argparse
import sys
import os
import math
import matplotlib.pyplot as plt

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NACCESS import NACCESS_atomic
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBIO import PDBIO, Select

class ResiduesDataLib():
    def __init__(self, fname):
        self.residue_data = {}
        try:
            fh = open(fname, "r")
        except OSError:
            print("#ERROR while loading library file (", fname, ")")
            sys.exit(2)
        for line in fh:
            if line[0] == '#':
                continue
            data = line.split()
            r = Residue(data)
            self.residue_data[r.id] = r
        self.nres = len(self.residue_data)

    def get_params(self, resid, atid):
        atom_id = resid + ':' + atid
        if atom_id in self.residue_data:
            return self.residue_data[atom_id]
        else:
            print("WARNING: atom not found in library (", atom_id, ')')
            return None

class Residue():
    def __init__(self,data):
        self.id = data[0]+':'+data[1]
        self.at_type = data[2]
        self.charge  = float(data[3])
        
class AtomType():
    def __init__(self, data):
        self.id   = data[0]
        self.eps  = float(data[1])
        self.sig  = float(data[2])
        self.mass = float(data[3])
        self.fsrf = float(data[4])
        self.rvdw = self.sig * 0.5612
        
class VdwParamset(): # Extracted from GELPI's github
    # Parameters for the VdW
    def __init__ (self, file_name):
        self.at_types = {}
        try:
            fh = open(file_name, "r")
        except OSError:
            print ("#ERROR while loading parameter file (", file_name, ")")
            sys.exit(2)
        for line in fh:
            if line[0] == '#':
                continue
            data = line.split()
            self.at_types[data[0]] = AtomType(data)
        self.ntypes = len(self.at_types)
        fh.close()

def atom_id(at):
    # Returns readable atom id
    return '{}.{}'.format(residue_id(at.get_parent()), at.id)

def residue_id(res):
    # Returns readable residue id
    return '{} {}{}'.format(res.get_resname(), res.get_parent().id, res.id[1])

def add_atom_parameters(st, res_lib, ff_params):
    for at in st.get_atoms():
        resname = at.get_parent().get_resname()
        params = res_lib.get_params(resname, at.id)
        if not params:
            print("WARNING: residue/atom pair not in library ("+atom_id(at) + ')')
            at.xtra['atom_type'] = at.element
            at.xtra['charge'] = 0
        else:
            at.xtra['atom_type'] = params.at_type
            at.xtra['charge'] = params.charge
        at.xtra['vdw'] = ff_params.at_types[at.xtra['atom_type']]

def MH_diel(r):
    # Mehler-Solmajer dielectric
    return 86.9525 / (1 - 7.7839 * math.exp(-0.3153 * r)) - 8.5525

def elec_int(at1, at2, r):
    # Electrostatic interaction energy between two atoms at r distance
    return 332.16 * at1.xtra['charge'] * at2.xtra['charge'] / MH_diel(r) / r

def vdw_int(at1, at2, r):
    # Vdw interaction energy between two atoms
    eps12 = math.sqrt(at1.xtra['vdw'].eps * at2.xtra['vdw'].eps)
    sig12_2 = at1.xtra['vdw'].sig * at2.xtra['vdw'].sig
    return 4 * eps12 * (sig12_2**6/r**12 - sig12_2**3/r**6)

def calc_solvation(res): # Computes the solvation energies for alanine
    # Returns the solvation energies (residue against other chains) for Ala atoms
    solv_ala = 0.
    for at in res.get_atoms():
        if 'EXP_NACCESS' not in at.xtra:
            continue
        s = float(at.xtra['EXP_NACCESS'])* at.xtra['vdw'].fsrf
        if at.id in ala_atoms:
            solv_ala += s
    return solv_ala

def calc_int_energies(st, res): # Computes the interaction energies for alanine
    # Returns interaction energies (residue against other chains) for Ala atoms
    elec_ala = 0.
    vdw_ala = 0.
    for at1 in res.get_atoms():
        for at2 in st.get_atoms():
        # Skip same chain atom pairs
            if at2.get_parent().get_parent() != res.get_parent():
                r = at1 - at2
                e = elec_int(at1, at2, r)
                if at1.id in ala_atoms: # GLY are included implicitly
                    elec_ala += e
                e = vdw_int(at1, at2, r)
                if at1.id in ala_atoms: # GLY are included implicitly
                    vdw_ala += e
    return elec_ala, vdw_ala

def get_interface(st, dist):
    # Detects interface residues within a distance(dist). Assumes two chains, i.e. a unique interface set per chain.
    select_ats = []
    for at in st.get_atoms():
        if at.element != 'H': # Skip Hydrogens to reduce time
            select_ats.append(at)
    nbsearch = NeighborSearch(select_ats)
    interface = {} # Sets are more efficient than lists. Use sets when order is not relevant
    for ch in st[0]:
        interface[ch.id] = set()

    for at1, at2 in nbsearch.search_all(dist):
        # Only different chains
        res1 = at1.get_parent() # Pointer to the parent node, here the parent node of an atom is the residue
        ch1 = res1.get_parent() # Here the parent node of a residue is a chain
        res2 = at2.get_parent()
        ch2 = res2.get_parent()
        if ch1 != ch2:
            interface[ch1.id].add(res1)
            interface[ch2.id].add(res2)
    return interface

# Load all the paramenters that we will need to execute all the previous functions
residue_library = ResiduesDataLib('/home/nuria/Downloads/biophysics.project/step2/aaLib.lib')
ff_params = VdwParamset('/home/nuria/Downloads/biophysics.project/step2/vdwprm.txt')
pdb_path = "/home/nuria/Downloads/biophysics.project/6m0j_fixed.pdb"
parser = PDBParser(PERMISSIVE=1)
st = parser.get_structure('st', pdb_path)
add_atom_parameters(st, residue_library, ff_params)

# Here we create the variable alanine atoms, with all the characteristics of alanine.
ala_atoms = {'N', 'H', 'CA', 'HA', 'C', 'O', 'CB', 'HB', 'HB1', 'HB2', 'HB3', 'HA1', 'HA2', 'HA3'}

# To be able to the the solvation related calculations:
NACCESS_BINARY = '/home/nuria/Downloads/biophysics.project/step2/soft/NACCESS/naccess'
srf = NACCESS_atomic(st[0], naccess_binary=NACCESS_BINARY)
io = PDBIO()
st_chains = {}
# Using BioIO trick to select chains
class SelectChain(Select):
    def __init__(self, chid):
        self.id = chid

    def accept_chain(self, chain):
        if chain.id == self.id:
            return 1
        else:
            return 0

for ch in st[0]:
    io.set_structure(st)
    io.save('tmp.pdb', SelectChain(ch.id))
    st_chains[ch.id] = parser.get_structure('stA', 'tmp.pdb')
    add_atom_parameters(st_chains[ch.id], residue_library, ff_params)
    srfA = NACCESS_atomic(st_chains[ch.id][0], naccess_binary=NACCESS_BINARY)
os.remove('tmp.pdb')


Telec = 0
Tvdw = 0
Tsolvation = 0
TsubunitA = 0
TsubunitE = 0
diccionari = {}

'''
In order to better understand energy related to the present amino acids, we use of dictionaries to enable 
us to take the energies with higher difference within eachother to compare the corresponding amino acids
'''

for model in st.get_models(): # We get the models of our structure
    chainA = model["A"] # we select the chain A of our structure
    chainE = model["E"] # we select the chain E of our structure

for residues in chainA.get_residues(): # For every residue in the chain A of our protein
    Telec += calc_int_energies(st[0], residues)[0] # Add electrostatic energy to the value of Ielec
    Tvdw += calc_int_energies(st[0], residues)[1] # Add Van der Waals energies to the value of Ielec
    Tsolvation += calc_solvation(residues) # Add the solvation for all the residues
    TsubunitA += calc_solvation(st[0][chainA.id][residues.id[1]]) # Add the solvarion at chain A
    diccionari[residues.id[1]] = calc_int_energies(st[0], residues)[0]+calc_int_energies(st[0], residues)[1]+calc_solvation(residues)-calc_solvation(st[0][chainA.id][residues.id[1]])

for residues in chainE.get_residues(): # For every residue in the chain E of our protein
    Telec += calc_int_energies(st[0], residues)[0] # For every residue in the interface in a distance of max 3.7
    Tvdw += calc_int_energies(st[0], residues)[1] # Add electrostatic energy to the value of Ielec
    Tsolvation += calc_solvation(residues) # Add the solvation for all the residues
    TsubunitE += calc_solvation(st[0][chainE.id][residues.id[1]]) # Add the solvarion at chain E
    diccionari[residues.id[1]] = calc_int_energies(st[0], residues)[0]+calc_int_energies(st[0], residues)[1]+calc_solvation(residues)-calc_solvation(st[0][chainA.id][residues.id[1]])

print("Total energy of all residues: ", Telec + Tvdw + Tsolvation - TsubunitA - TsubunitE) # Global energy


Ielec = 0
Ivdw = 0
Isolvation = 0
IsubunitA = 0
IsubunitE = 0

for residues in get_interface(st, 3.7)["A"]: # For every residue in the interface in a distance of max 3.7
    Ielec += calc_int_energies(st[0], residues)[0] # Add electrostatic energy to the value of Ielec
    Ivdw += calc_int_energies(st[0], residues)[1] # Add Van der Waals energies to the value of Ielec
    Isolvation += calc_solvation(residues) # Add the solvation for all the residues
    IsubunitA += calc_solvation(st[0][chainA.id][residues.id[1]]) # Add the solvarion at chain A

for residues in get_interface(st, 3.7)["E"]: # For every residue in the interface in a distance of max 3.7
    Ielec += calc_int_energies(st[0], residues)[0] # Add electrostatic energy to the value of Ielec
    Ivdw += calc_int_energies(st[0], residues)[1] # Add Van der Waals energies to the value of Ielec
    Isolvation += calc_solvation(residues) # Add the solvation for all the residues
    IsubunitE += calc_solvation(st[0][chainE.id][residues.id[1]]) # Add the solvarion at chain E

print("Total energy of interface residues: ", Ielec + Ivdw + Isolvation - IsubunitA - IsubunitE) # Interface energy

# Create two lists and select the amino acids that are more and less relevant in terms of energy
pos_important = []
neg_important = []
for key, values in diccionari.items():
    if values > 1:
        pos_important.append(key)
    elif values < -3.25:
        neg_important.append(key)

# Print to the terminal the residues with more positive energy
print('Positive energy residues: ')
for x in pos_important:
    for residues in st.get_residues():
        if residues.id[1] == x:
            print(residues.get_resname())

print('')
# Print to the terminal the residues with more negative energy
print('Negative energy residues: ')
for x in neg_important:
    for residues in st.get_residues():
        if residues.id[1] == x:
            print(residues.get_resname())

# Observe all the energies of the residues in a plot
names = list(diccionari.keys())
values = list(diccionari.values())
plt.bar(range(len(diccionari)), values, tick_label=names)
plt.show()
