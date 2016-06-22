#!/usr/bin/python

import numpy as np
import sys
from scipy import stats

def LIE(f1, f2):
	data_1 = np.genfromtxt(f1)[:,1:]
	eElec_1 = data_1[:,0]
	eVdw_1 = data_1[:,1]
	meanElec_1 = np.mean(eElec_1)
	stdErr_1e = stats.sem(eElec_1)
	meanVdw_1 = np.mean(eVdw_1)
	stdErr_1v = stats.sem(eVdw_1)

	data_2 = np.genfromtxt(f2)[:,1:]
	eElec_2 = data_2[:,0]
	eVdw_2 = data_2[:,1]
	meanElec_2 = np.mean(eElec_2)
	stdErr_2e = stats.sem(eElec_2)
	meanVdw_2 = np.mean(eVdw_2)
	stdErr_2v = stats.sem(eVdw_2)

	delElec = meanElec_1 - meanElec_2
	errE = np.sqrt(stdErr_1e**2 + stdErr_2e**2)
	delVdw = meanVdw_1 - meanVdw_2
	errV = np.sqrt(stdErr_1v**2 + stdErr_2v**2)
	dGbind = (0.16 * delVdw + 0.5 * delElec)
	error = np.sqrt((0.16*errV)**2 + (0.5 * errE)**2)

	return dGbind,error

print LIE(sys.argv[1],sys.argv[2])
