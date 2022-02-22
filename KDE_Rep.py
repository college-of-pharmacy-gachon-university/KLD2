#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import numpy as np
from scipy.stats import *
import sys

filename = sys.argv[1]

result = []
[result.extend(el) for el in pd.read_csv(filename + '.zip', index_col=0).round(decimals = 4).values.tolist()] 
x = np.linspace(0,1,101)
rep_pdf = gaussian_kde(result).pdf(x)
pd.DataFrame([rep_pdf, x]).T.to_csv(filename[5:-4]+'.csv')

