import os
import numpy as np
import pandas as pd

import tempfile

from rdkit import Chem
from rdkit.Chem import AllChem
import deepchem as dc

from deepchem.utils import download_url, load_from_disk 

from openmm.app import PDBFile
from pdbfixer import PDBFixer
import mdtraj as md
import nglview
from deepchem.utils.vina_utils import prepare_inputs
from IPython.display import display, Image
from pkg_resources import resource_filename
import os
folder_path = './upload'
file_list = os.listdir(folder_path)
for filename in file_list:
    if "protein.pdb" in filename:
        fixer = PDBFixer(filename="./upload/"+filename)
        filename = filename.replace("_protein.pdb", "")

        PDBFile.writeFile(fixer.topology, fixer.positions, open('templates/display/%s.pdb' % (filename), 'w'))
        p, m = None, None
        p, m = prepare_inputs('templates/display/%s.pdb' % (filename), "C")
        if p and m:  # protein and molecule are readable by RDKit
            Chem.rdmolfiles.MolToPDBFile(p, 'templates/display/protein_%s.pdb' % (filename))
            
                
        protein_mdtraj = md.load_pdb( 'templates/display/protein_%s.pdb' % (filename))
        p = nglview.show_mdtraj(protein_mdtraj)
        nglview.write_html("templates/display/protein_%s.html" % (filename),[p])
    