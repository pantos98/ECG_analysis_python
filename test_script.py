import numpy as np 
import os 
import csv
import pandas as pd
import wfdb

record = wfdb.rdrecord('14046')
signals, fields = wfdb.rdsamp('14046', channels=[0], sampfrom=100, sampto=15000)
wfdb.plot_wfdb(record=record, title='Example signals')

