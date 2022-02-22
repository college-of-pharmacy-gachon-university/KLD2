# A Random Forest Model for DTI(Drug-Target Interaction) Prediction via KL-Divergence

This repository contains the useful script, used in our manuscript. 
Collection of scripts for Kullback-Leibler divergence aided Random forest model


This repository is divided into four parts;

- The Molecular Fingerprints from E3FP descriptor using the E3FP python interface (https://github.com/keiserlab/e3fp).
- Calculate Similarity score(e3fp Rdkit)
- Density Estimation scheme(KDE) and KL-Divergence Calculation 
- The Drug-Target classifier based on Random-Forest Classifier 

(Requirements)

a) This script is designed to run under Python 3.6 
e3fp, NumPy, SciPy, Pandas, Matplotlib, RDKit Scikit Learn, Seaboarn 

(usage)

1. Generate 3d-Fingerprint(by Molecular descriptor) using GenerateFingerprint.ipynb
2. Calculate similarity score using Calc_Sim_score.py
3. Generate Densities using KDE_Rep.py
4. Calculate Divergence using CalcKLD.py
5. Train and Test the Dataset using Random-forest Classifier (KLD-RF.ipynb)

(Result)

The KL-divergence from estimated density(by KDE) works as an 
feature for RF-Classifier that characterizes the structures of pharmachological target.


(Acknowledgement)

This Study was Conducted as a Research Project of M.H.Kim Lab in the School of Pharmacy, Gachon University.
