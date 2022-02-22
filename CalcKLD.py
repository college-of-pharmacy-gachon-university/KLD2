#!/usr/bin/env python
# coding: utf-8



import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import *
import sys

def CalcKLD(P, Q):
    epsilon = 0.0001
    P = P + epsilon
    Q = Q+ epsilon
    divergence = np.sum(P*np.log(P/Q))
    return divergence

pdf_data = pd.read_csv('pdf_data.csv')
filename = sys.argv[1]
cross_pdf = pd.read_csv(filename, index_col=0)
KL = np.zeros ((len(cross_pdf.columns),len(pdf_data.columns)))
KLD = pd.DataFrame(KL).T
KLD.index = pdf_data.columns
KLD.columns = cross_pdf.columns
for i in cross_pdf.columns:
    for j in pdf_data.columns:
        KLD[i][j]=CalcKLD(cross_pdf[i].values/np.sum(cross_pdf[i].values), pdf_data[j].values/np.sum(pdf_data[j].values))
        
result = KLD.T.add_prefix(filename[10:-8] + '_KLD_')
result.to_csv(filename[10:-8] + ' KLDTable.csv')



plt.figure(figsize=(18,6))
for i in range(1, 19):
    plt.subplot(3, 6, i)
    plt.hist(result.iloc[:, i-1].values, bins = 100, range = [0, 2], label = str(i))
    #plt.x_label(str(i))
plt.title(filename[10:-8] + '_KLD_')
plt.savefig(filename[10:-8])

