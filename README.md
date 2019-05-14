# atom-types-bc-thesis

About

ATTYC - ATom TYpe Classification Documentation

Author: Radka Sedláková, 05-15-2019
	460676@mail.muni.cz
	www.github.com/radenka/atom-types-bc-thesis
	
This program was developed as implementation part of "Atom types in methods for calculation of partial 
atomic charges" thesis at Faculty of Science, Masaryk University in Brno, the Czech Republic. 


ATTYC library provides atom type classification of SDF and PDB files based on hybridization, highest bond 
order, structural motif membership and bonding partners of atom; special classification 'peptide' for PDB
files is also added. These atom type classifications can be used for parameterization of empirical methods 
for calculation of atomic partial charges.


Directory 'attyc' contains program implementation. Atom type classification can be run through command line 
using module main.py, located in this directory. In case of calling ATTYC from different program, proper 
integration of ATTYC library into program is required. After that, classification can be run using command 
'attyc.classify_atoms(input_file, classifier, file_output, screen_output)'. In case of running ATTYC through
command line, command line arguments are passed on this function.


Running the program

To show help type './main.py --help' or './main.py -h' to command line in this directory.

usage: ATTYC: ATom TYpe Classification [-h] [--file_output] [--screen_output]
                                       input_file classifier

positional arguments:
  input_file       Path to file to process. 
  classifier       Atomic property that atom types assigned to atoms are
                   derived from. Select one of these:'substruct' 'hbo'
                   'peptide' 'hybrid' 'partners'.

optional arguments:
  -h, --help       show this help message and exit
  --file_output    Atom types will be written into text file.
  --screen_output  Atom types will be printed on screen.

When running classification from different program, always assign boolean to output arguments and use molecule
files with appropriate file extensions (molset.sdf, molecule.pdb). ATTYC raises exception otherwise. When 
setting --file_output to True, attyc_outputs directory is created automatically containing text file with 
assigned atom types. All further outputs will be saved into this directory.


Classifiers 'peptide' and 'substruct'

Classifiers 'substruct' and 'peptide' classify atoms using data from external text files where appropriate atom 
types of specific atoms are defined (substruct: SMARTS_atom_types.txt, peptide: PDB_atom_types.txt). 

Classifications based on data in SMARTS_atom_types.txt file and PDB_atom_types.txt are set as default, however, 
these classifications resulted computationally demanding. Simplified atom type classification text files can be 
found in 'simplified_classifiers' directory. Replace original classification files located in 'attyc' directory 
and rename them to 'SMARTS_atom_types.txt' and 'PDB_atom_types.txt' to run simplified classification. 
To run simplified 'substruct' classification correctly:

file substruct.py:
Follow instructions in comments in methods "classifify_remaining_H(self, hydrogen)" and 
					   "analyze_aromatic_rings(self, molecule, atom_types)".


Wish you good luck with ATTYC!

