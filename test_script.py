import numpy as np 
import os 
import csv
import pandas as pd
import wfdb
import hr_calc_method as hr

external_database=1

if external_database:
    record = wfdb.rdrecord('ecg_database/14046', channels=[0], sampfrom=00, sampto=30000)
    annotation = wfdb.rdann('ecg_database/14046', 'atr', sampfrom=00, sampto=30000)
    wfdb.plot_wfdb(record=record, annotation=annotation)
    signals, fields = wfdb.rdsamp('ecg_database/14046', channels=[0], sampfrom=00, sampto=30000)
    ecg = np.zeros([signals.shape[0],1])
    for i in range(signals.shape[0]):
        ecg[i] = signals[i]
else:
    directory = "Z:\\ETH_RASMDEP_PILOT\\"
    folder = directory + "S002\\vivalnk_vv330_ecg\\20211110\\"
    df = pd.read_csv(folder + "20211110_1800.csv.gz")
    ecg = np.array([df["value.ecg"]])
    timestamps = np.array([df["value.time"]])

method = hr.hrMethod(128)
method.inputECG(ecg)
pos, sign, en_thres = method.calculateHR()
print(method.hrv)
method.plotDebug()