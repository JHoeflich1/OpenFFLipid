from openff.toolkit import Molecule, Quantity, RDKitToolkitWrapper, Topology, unit #type:ignore
from openff.units import Quantity, unit#type:ignore
import pandas as pd #type:ignore
import numpy as np
import mdtraj#type:ignore
import os
import shutil 
import time
import subprocess
from typing import List
import json
from datetime import datetime
from openff.interchange import Interchange#type:ignore
from openff.toolkit import ForceField, Molecule, Topology#type:ignore


path = 'testSM102.pdb'
smiles = 'OCCN(CCCCCCCC(OC(CCCCCCCC)CCCCCCCC)=O)CCCCCC(OCCCCCCCCCCC)=O'

molecule = Topology.from_pdb(path)  # Load the molecule from the PDB file
print(f"Lipid: , Atoms: {len(molecule.atoms)}")
