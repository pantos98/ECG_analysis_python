import numpy as np 
import os 
import csv
import pandas as pd
import wfdb
import hr_calc_method as hr

signals, fields = wfdb.rdsamp('14046', channels=[0], sampfrom=15000, sampto=50000)
ecg = np.zeros([signals.shape[0],1])
for i in range(signals.shape[0]):
    ecg[i] = signals[i]
    
method = hr.hrMethod(128)
method.inputECG(ecg)
pos, sign, en_thres = method.calculateHR()
print(method.hrv)
method.plotDebug()