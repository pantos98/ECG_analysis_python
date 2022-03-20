import numpy as np 
import os 
import csv
import pandas as pd
import wfdb
import hr_calc_method as hr

record = wfdb.rdrecord('14046')
signals, fields = wfdb.rdsamp('14046', channels=[0], sampfrom=100, sampto=15000)
ecg = np.zeros([signals.shape[0],1])

for i in range(signals.shape[0]):
    ecg[i] = signals[i]
#wfdb.plot_wfdb(record=record, title='Example signals')
method = hr.hrMethod(128)
method.inputECG(ecg)
method.calculateHR()
print(method)
