#!/usr/bin/env python
# coding: utf-8

# In[11]:

import sys
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
from rdkit.Chem import DataStructs
from rdkit.Chem import AllChem
IPythonConsole.ipython_useSVG=True


def CrossTC(classi, classj):
    TC= np.zeros( shape=(len(classi), len(classj)), dtype = np.float16)
    for i in range( len(classi) ):
        for j in range( len(classj) ):
            SimTC = DataStructs.TanimotoSimilarity( classi[i], classj[j] )
            TC[i][j] = SimTC
    return TC



#filename = sys.argv[1]
#Name = sys.argv[1]
#filename=['[2]', '[3]', '[5]'. '[10]'. '[11]', '[13]', '[14]']

mols2 = [ mol for mol in Chem.SDMolSupplier( "[2].sdf", removeHs=False ) ]
mols3 = [ mol for mol in Chem.SDMolSupplier( "[3].sdf", removeHs=False ) ]
mols5 = [ mol for mol in Chem.SDMolSupplier( "[5].sdf", removeHs=False ) ]
mols10 = [ mol for mol in Chem.SDMolSupplier( "[10].sdf", removeHs=False ) ]
mols11 = [ mol for mol in Chem.SDMolSupplier( "[11].sdf", removeHs=False ) ]
mols13 = [ mol for mol in Chem.SDMolSupplier( "[13].sdf", removeHs=False ) ]
mols14 = [ mol for mol in Chem.SDMolSupplier( "[14].sdf", removeHs=False ) ]
mols5 = [ mol for mol in Chem.SDMolSupplier( "[5].sdf", removeHs=False ) ]
mols11 = [ mol for mol in Chem.SDMolSupplier( "[11].sdf", removeHs=False ) ]


fps2 = [ AllChem.GetMorganFingerprintAsBitVect(mol,2) for mol in mols2 ]
fps3 = [ AllChem.GetMorganFingerprintAsBitVect(mol,2) for mol in mols3 ]
fps5 = [ AllChem.GetMorganFingerprintAsBitVect(mol,2) for mol in mols5 ]
fps10 = [ AllChem.GetMorganFingerprintAsBitVect(mol,2) for mol in mols10 ]
fps11 = [ AllChem.GetMorganFingerprintAsBitVect(mol,2) for mol in mols11 ]
fps13 = [ AllChem.GetMorganFingerprintAsBitVect(mol,2) for mol in mols13 ]
fps14 = [ AllChem.GetMorganFingerprintAsBitVect(mol,2) for mol in mols14 ]

fps5 = [ AllChem.GetMorganFingerprintAsBitVect(mol,2) for mol in mols5 ]
fps11 = [ AllChem.GetMorganFingerprintAsBitVect(mol,2) for mol in mols11 ]



names2 = [x.GetProp("_Name") for x in mols2]  
names3 = [x.GetProp("_Name") for x in mols3]  
names5 = [x.GetProp("_Name") for x in mols5]  
names10 = [x.GetProp("_Name") for x in mols10]  
names11 = [x.GetProp("_Name") for x in mols11]  
names13 = [x.GetProp("_Name") for x in mols13]  
names14 = [x.GetProp("_Name") for x in mols14]  

names5 = [x.GetProp("_Name") for x in mols5]  
names11 = [x.GetProp("_Name") for x in mols11] 


pd.DataFrame(CrossTC(fps3, fps3), index = names3, columns =  names3).to_csv('cross3-3'+'.csv.zip')
pd.DataFrame(CrossTC(fps13, fps13), index = names13, columns =  names13).to_csv('cross13-13'+'.csv.zip')
pd.DataFrame(CrossTC(fps3, fps13), index = names3, columns =  names13).to_csv('cross3-13'+'.csv.zip')
pd.DataFrame(CrossTC(fps13, fps3), index = names13, columns =  names3).to_csv('cross13-3'+'.csv.zip')
pd.DataFrame(CrossTC(fps2, fps2), index = names2, columns =  names2).to_csv('cross2-2'+'.csv.zip')
pd.DataFrame(CrossTC(fps10, fps10), index = names10, columns =  names10).to_csv('cross10-10'+'.csv.zip')
pd.DataFrame(CrossTC(fps14, fps14), index = names14, columns =  names14).to_csv('cross14-14'+'.csv.zip')
pd.DataFrame(CrossTC(fps2, fps10), index = names2, columns =  names10).to_csv('cross2-10'+'.csv.zip')
pd.DataFrame(CrossTC(fps2, fps14), index = names2, columns =  names14).to_csv('cross2-14'+'.csv.zip')
pd.DataFrame(CrossTC(fps10, fps2), index = names10, columns =  names2).to_csv('cross10-2'+'.csv.zip')
pd.DataFrame(CrossTC(fps10, fps14), index = names10, columns =  names14).to_csv('cross10-14'+'.csv.zip')
pd.DataFrame(CrossTC(fps14, fps2), index = names14, columns =  names2).to_csv('cross14-2'+'.csv.zip')
pd.DataFrame(CrossTC(fps14, fps10), index = names14, columns =  names10).to_csv('cross14-10'+'.csv.zip')

pd.DataFrame(CrossTC(fps5, fps5), index = names5, columns =  names5).to_csv('cross5-5'+'.csv.zip')
pd.DataFrame(CrossTC(fps5, fps11), index = names5, columns =  names11).to_csv('cross5-11'+'.csv.zip')
pd.DataFrame(CrossTC(fps11, fps5), index = names11, columns =  names5).to_csv('cross11-5'+'.csv.zip')
pd.DataFrame(CrossTC(fps11, fps11), index = names11, columns =  names11).to_csv('cross11-11'+'.csv.zip')

