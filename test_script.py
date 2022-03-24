import numpy as np 
import os 
import csv
import pandas as pd
import wfdb
import hr_calc_method as hr

external_database=1

if external_database:
    start = 00000
    end = 50000
    name = 'ecg_database/14184'
    record = wfdb.rdrecord(name, channels=[0], sampfrom=start, sampto=end)
    annotation = wfdb.rdann(name, 'atr', sampfrom=start, sampto=end)
    annotation.sample = annotation.sample - start
    wfdb.plot_wfdb(record=record, annotation=annotation)
    signals, fields = wfdb.rdsamp(name, channels=[0], sampfrom=start, sampto=end)
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
method.plotDebug()