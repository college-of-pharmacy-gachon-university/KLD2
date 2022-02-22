#!/usr/bin/env python
# coding: utf-8



import pandas as pd
import numpy as np
from rdkit import Chem
from e3fp.fingerprint.generate import fp, fprints_dict_from_mol, fprints_dict_from_sdf
from e3fp.conformer.generate import generate_conformers
from e3fp.pipeline import fprints_from_sdf, fprints_from_mol
from e3fp.fingerprint.fprint import save, savez
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
from rdkit.Chem import DataStructs
from rdkit.Chem import AllChem
IPythonConsole.ipython_useSVG=True
from e3fp.fingerprint.db import FingerprintDatabase
from e3fp.fingerprint.fprint import Fingerprint
import numpy as np



def CrossTC(classi, classj):
    TC= np.zeros( shape=(len(classi), len(classj)), dtype = np.float16)
    for i in range( len(classi) ):
        for j in range( len(classj) ):
            e3fpTC = DataStructs.TanimotoSimilarity( classi[i], classj[j] )
            TC[i][j] = e3fpTC
    return TC




db1 = FingerprintDatabase.load("Class1fingerprint.fps.bz2")
db2 =  FingerprintDatabase.load("Class2fingerprint.fps.bz2")
db3 = FingerprintDatabase.load("Class3fingerprint.fps.bz2")
db4 =  FingerprintDatabase.load("Class4fingerprint.fps.bz2")
db5 =  FingerprintDatabase.load("Class5fingerprint.fps.bz2")
db6 = FingerprintDatabase.load("Class6fingerprint.fps.bz2")
db7 =  FingerprintDatabase.load("Class7fingerprint.fps.bz2")
db8 =  FingerprintDatabase.load("Class8fingerprint.fps.bz2")
db9 =  FingerprintDatabase.load("Class9fingerprint.fps.bz2")
db10 =  FingerprintDatabase.load("Class10fingerprint.fps.bz2")
db11 =  FingerprintDatabase.load("Class11fingerprint.fps.bz2")
db12 =  FingerprintDatabase.load("Class12fingerprint.fps.bz2")
db13 =  FingerprintDatabase.load("Class13fingerprint.fps.bz2")
db14 =  FingerprintDatabase.load("Class14fingerprint.fps.bz2")
db15 =  FingerprintDatabase.load("Class15fingerprint.fps.bz2")
db16 =  FingerprintDatabase.load("Class16fingerprint.fps.bz2")
db17 =  FingerprintDatabase.load("Class17fingerprint.fps.bz2")
#db18 =  FingerprintDatabase.load("Class18fingerprint.fps.bz2")



binfp1 = [fp.fold().to_rdkit() for fp in db1]
binfp2 = [fp.fold().to_rdkit() for fp in db2]
binfp3 = [fp.fold().to_rdkit() for fp in db3]
binfp4 = [ fp.fold().to_rdkit() for fp in db4 ]
binfp5 = [fp.fold().to_rdkit() for fp in db5]
binfp6 = [fp.fold().to_rdkit() for fp in db6]
binfp7 = [fp.fold().to_rdkit() for fp in db7]
binfp8 = [ fp.fold().to_rdkit() for fp in db8 ]
binfp9 = [fp.fold().to_rdkit() for fp in db9]
binfp10 = [ fp.fold().to_rdkit() for fp in db10 ]
binfp11 = [ fp.fold().to_rdkit() for fp in db11 ]
binfp12 = [ fp.fold().to_rdkit() for fp in db12 ]
binfp13 = [ fp.fold().to_rdkit() for fp in db13 ]
binfp14 = [ fp.fold().to_rdkit() for fp in db14 ]
binfp15 = [ fp.fold().to_rdkit() for fp in db15 ]
binfp16 = [ fp.fold().to_rdkit() for fp in db16 ]
binfp17 = [ fp.fold().to_rdkit() for fp in db17 ]
#binfp18 = [ fp.fold().to_rdkit() for fp in db18 ]



pd.DataFrame(CrossTC(binfp1, binfp1),index = db1.fp_names, columns = db1.fp_names).to_csv('cross1-1.csv.zip')
pd.DataFrame(CrossTC(binfp1, binfp2),index = db1.fp_names, columns = db2.fp_names).to_csv('cross1-2.csv.zip')
pd.DataFrame(CrossTC(binfp1, binfp3),index = db1.fp_names, columns = db3.fp_names).to_csv('cross1-3.csv.zip')
pd.DataFrame(CrossTC(binfp1, binfp4),index = db1.fp_names, columns = db4.fp_names).to_csv('cross1-4.csv.zip')
pd.DataFrame(CrossTC(binfp1, binfp5),index = db1.fp_names, columns = db5.fp_names).to_csv('cross1-5.csv.zip')
pd.DataFrame(CrossTC(binfp1, binfp6),index = db1.fp_names, columns = db6.fp_names).to_csv('cross1-6.csv.zip')
pd.DataFrame(CrossTC(binfp1, binfp7),index = db1.fp_names, columns = db7.fp_names).to_csv('cross1-7.csv.zip')
pd.DataFrame(CrossTC(binfp1, binfp8),index = db1.fp_names, columns = db8.fp_names).to_csv('cross1-8.csv.zip')
pd.DataFrame(CrossTC(binfp1, binfp9),index = db1.fp_names, columns = db9.fp_names).to_csv('cross1-9.csv.zip')
pd.DataFrame(CrossTC(binfp1, binfp10),index = db1.fp_names, columns = db10.fp_names).to_csv('cross1-10.csv.zip')
pd.DataFrame(CrossTC(binfp1, binfp11),index = db1.fp_names, columns = db11.fp_names).to_csv('cross1-11.csv.zip')
pd.DataFrame(CrossTC(binfp1, binfp12),index = db1.fp_names, columns = db12.fp_names).to_csv('cross1-12.csv.zip')
pd.DataFrame(CrossTC(binfp1, binfp13),index = db1.fp_names, columns = db13.fp_names).to_csv('cross1-13.csv.zip')
pd.DataFrame(CrossTC(binfp1, binfp14),index = db1.fp_names, columns = db14.fp_names).to_csv('cross1-14.csv.zip')
pd.DataFrame(CrossTC(binfp1, binfp15),index = db1.fp_names, columns = db15.fp_names).to_csv('cross1-15.csv.zip')
pd.DataFrame(CrossTC(binfp1, binfp16),index = db1.fp_names, columns = db16.fp_names).to_csv('cross1-16.csv.zip')
pd.DataFrame(CrossTC(binfp1, binfp17),index = db1.fp_names, columns = db17.fp_names).to_csv('cross1-17.csv.zip')
#pd.DataFrame(CrossTC(binfp1, binfp18),index = db1.fp_names, columns = db18.fp_names).to_csv('cross1-18.csv.zip')





