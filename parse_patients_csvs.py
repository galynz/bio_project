# -*- coding: utf-8 -*-
"""
Created on Thu Sep 08 11:39:28 2016

@author: gal
"""

import pandas as pd
import glob
import os

def parse_patients_csvs(csv_paths):
    all_files = glob.glob(csv_paths)
    frame = pd.DataFrame()
    sum_frame = pd.DataFrame([])
    l = []
    l_sum = []
    for path in all_files:
        filename = os.path.basename(path)
        pos = filename.find('_2016')
        cancer = filename[:pos]
        df = pd.read_csv(path,index_col=None, header=0)
        df['hr_deficient'] = df.HR_mutated > 0
        df['cancer_location'] = cancer
        median = df['Mutations_Count'].median()/30.0
        mean = df['Mutations_Count'].mean()/30.0
        ix = df['hr_deficient'] == True
        deficient_median = df[ix]['Mutations_Count'].median()/30.0
        proficient_median = df[~ix]['Mutations_Count'].median()/30.0
        deficient_mean = df[ix]['Mutations_Count'].mean()/30.0
        proficient_mean = df[~ix]['Mutations_Count'].mean()/30.0
        print cancer, median, mean, deficient_median, proficient_median, deficient_mean, proficient_mean
        l.append(df)
    frame = pd.concat(l)
    return frame