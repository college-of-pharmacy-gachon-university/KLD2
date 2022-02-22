# A Random Forest Model for DTI(Drug-Target Interaction) Prediction via KL-Divergence

Collection of scripts for distinguishing two Pharmacological classes using Kullback-Leibler divergence from Tanimoto Coefficient structure of each classes.


(Requirements)

This script is designed to run under Python 3.6 and Anaconda version 5.2.0


(usage)

1. Generate pdf (probability density function) driven from Gaussian mixture model(GMM) (Using GMM.ipynb) : The pdf represents each of pharmacological classes
2. Calculate and Visualize Kullback-Leibler Divergence of Pharmacological Class(64+64, 64+10003, 10003+10003) from gaussian mixture pdfs (Using KLD.ipynb)

(Result)

The final result form KL-divergence works as an informatical comparison method that characterizes the structures of 3D-similarity scores.
