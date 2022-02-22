#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import numpy as np
from rdkit import Chem
import random
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
from rdkit.Chem import DataStructs
from rdkit.Chem import AllChem
IPythonConsole.ipython_useSVG=True
import sys




name = sys.argv[1]
mols = [ mol for mol in Chem.SDMolSupplier( name + ".sdf", removeHs=False ) ]
random_mol = random.sample(mols, 15000)
w = Chem.SDWriter(name + '15KRandom'+'.sdf')
for m in random_mol:
    w.write(m)
w.close()




