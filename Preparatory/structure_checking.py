# A complete checking analysis of a single structure. This script is modified from an original script:
# https://github.com/jlgelpi/BioPhysics/blob/master/Notebooks/6m0j_check.ipynb

# 1. We import the biobb_structure_checking module and the constants module from this module.
import biobb_structure_checking
import biobb_structure_checking.constants as cts
from biobb_structure_checking.structure_checking import StructureChecking

# 2. We get the base path of the biobb_structure_checking module and set the arguments for the StructureChecking class.
base_dir_path=biobb_structure_checking.__path__[0]
args = cts.set_defaults(base_dir_path,{'notebook':True})
with open(args['commands_help_path']) as help_file:
    print(help_file.read())

# 3. We define all the input and output paths.
base_path = '/home/nuria/Downloads/'
args['output_format'] = 'pdb'
args['keep_canonical'] = False
args['input_structure_path'] = base_path + '6m0j.pdb'
args['output_structure_path'] = base_path + '6m0j_fixed.pdb'
args['output_structure_path_charges'] = base_path + '6m0j_fixed.pdbqt'
args['debug'] = False
args['verbose'] = False

# 4. Initializing checking engine, loading structure and showing statistics
st_c = StructureChecking(base_dir_path, args)

# 5. Checks for the presence of models, chains and residues in the structure
st_c.models()
st_c.chains()
st_c.altloc()
st_c.altloc('occupancy')
st_c.altloc()

# 6. Checks for the presence of heteroatoms in the structure
st_c.metals()
st_c.ligands()
st_c.ligands('All')
st_c.ligands()

# 7. Detects and remove hydrogen atoms and water. 
st_c.rem_hydrogen()
st_c.water()
st_c.water("yes")

# 8. Amide terminal atoms can be labelled incorrectly, so observe possible fixes.
st_c.amide()
st_c.amide('all')
st_c.amide('A42,A103')
st_c.amide('E394')

# 9. Side chains of Thr and Ile are chiral, incorrect atom labelling lead to the wrong chirality.  
st_c.chiral()

# 10. Detects and fixes several problems with the backbone use
st_c.backbone()
st_c.backbone('--fix_atoms All --fix_chain none --add_caps none')

# 11. Detects and re-built missing protein side chains.   
st_c.fixside()

# 12. Detects possible -S-S- bonds based on distance criteria.
st_c.getss()
st_c.getss('all')

# 13. Add Hydrogen Atoms
st_c.add_hydrogen()
st_c.add_hydrogen('auto')

# 14. Detects steric clashes based on distance criteria.  
st_c.clashes()
st_c.checkall()

# 15. Save the fixed structure
st_c._save_structure(args['output_structure_path'])
st_c.rem_hydrogen('yes')
import os
os.system('check_structure -i ' + args['output_structure_path'] + ' -o ' + args['output_structure_path_charges'] + ' add_hydrogen --add_charges --add_mode auto')
Footer
