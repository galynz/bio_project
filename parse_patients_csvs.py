# -*- coding: utf-8 -*-
"""
Created on Thu Sep 08 11:39:28 2016

@author: gal
"""

import pandas as pd
import glob
import os
import scipy

import matplotlib
matplotlib.use('Agg')
from plotly.offline import download_plotlyjs, init_notebook_mode,  plot
import plotly.graph_objs as go
import plotly.tools as tls   

def parse_patients_csvs(csv_paths, cancer_order):
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
        df['hr_deficient_extended'] = df.where((df.HR_mutated > 0) | (df.NER_mutated > 0) | (df.MMR_mutated  > 0), None).Tumor_Sample_Barcode.apply(lambda x: x!= None)
        #df['hr_deficient_extended'] = df.HR_mutated > 0 | df.NER_mutated > 0 | df.MMR_mutated > 0
        df['cancer_location'] = cancer
        df['cancer_barcode'] = df.Tumor_Sample_Barcode.apply(lambda x: "_".join([cancer, x]))
        median = df['Mutations_Count_per_megabase'].median()
        mean = df['Mutations_Count_per_megabase'].mean()
        if cancer not in cancer_order:
            continue
        cancer_pos = cancer_order.index(cancer)
        df['cancer_pos'] = cancer_pos
        ix = df['hr_deficient'] == True
        deficient_median = df[ix]['Mutations_Count_per_megabase'].median()
        proficient_median = df[~ix]['Mutations_Count_per_megabase'].median()
        deficient_mean = df[ix]['Mutations_Count_per_megabase'].mean()
        proficient_mean = df[~ix]['Mutations_Count_per_megabase'].mean()
#        for gene in ["ATM", "ATRX", "BRIP1", "CHEK2", "FANCA", "FANCC", "BRCA1","BRCA2","FANCD2", "FANCE", "FANCF","FANCG", "NBN", "PTEN", "U2AF1"]:
#            df[gene] = df[gene].apply(lambda x: x==True)
        df['SUB_HR'] = scipy.logical_or.reduce((df.ATM, df.ATRX, df.BRIP1, df.FANCA, df.BRCA1, df.BRCA2, df.FANCD2, df.PTEN, df.ATR))
        print cancer, median, mean, deficient_median, proficient_median, deficient_mean, proficient_mean
        l.append(df)
    frame = pd.concat(l)
    return frame
    
def plot_mutation_load_hr_box(df, cancer_order, filename_format, column_name='hr_deficient'):
    for cancer in cancer_order:
        ix = df['cancer_location'] == cancer
        df_cancer = df[ix]
        ix_HRD = df_cancer[column_name] == True
        pvalue = tls.scipy.stats.ttest_ind(df_cancer[ix_HRD]['Mutations_Count_per_megabase'], df_cancer[~ix_HRD]['Mutations_Count_per_megabase'], equal_var=False).pvalue
        print cancer, pvalue
        trace1 = go.Box(y=df_cancer[ix_HRD]['Mutations_Count_per_megabase'], name='HRD', 
                        boxpoints='all', jitter=0.5,whiskerwidth=0.2, 
                        fillcolor='rgba(93, 164, 214, 0.5)', 
                        marker=dict(size=2),line=dict(width=1, 
                        color='rgba(31, 119, 180, 0.5)'),
                        boxmean=True)
        trace2 = go.Box(y=df_cancer[~ix_HRD]['Mutations_Count_per_megabase'], name='HRP', 
                        boxpoints='all', jitter=0.5,whiskerwidth=0.2, 
                        fillcolor='rgba(44, 160, 101, 0.5)', 
                        marker=dict(size=2),line=dict(width=1, 
                        color='rgba(44, 160, 44, 0.5)'),
                        boxmean=True)
        layout = go.Layout(title = 'mutatio load HRD/HRP - %s <br>pvalue = %s' % (cancer, pvalue))
        fig = go.Figure(data=[trace1, trace2], layout=layout)
        plot(fig, auto_open=False, filename = filename_format % cancer)
             #filename=r"d:\bio_project\outputs\outputs_20161006\box_plots_HR_new\%s_mutation_load_HR.html" % cancer)
             
def plot_mutation_bar_heatmap(df, cancer_order, output_dir):
    for cancer in cancer_order:
        ix = df['cancer_location'] == cancer
        df_cancer = df[ix]
        df_sorted = df_cancer.sort_values("Mutations_Count_per_megabase")
        trace1 = go.Bar(x=df_sorted.Tumor_Sample_Barcode, y=df_sorted.Mutations_Count_per_megabase)
        trace2 = go.Heatmap(z=[df_sorted.ATM.apply(int), df_sorted.ATRX.apply(int), df_sorted.BRIP1.apply(int), df_sorted.FANCA.apply(int), df_sorted.BRCA1.apply(int), df_sorted.BRCA2.apply(int), df_sorted.FANCD2.apply(int), df_sorted.PTEN.apply(int), df_sorted.ATR.apply(int), df_sorted.CHEK2.apply(int), df_sorted.FANCC.apply(int), df_sorted.FANCE.apply(int), df_sorted.FANCF.apply(int), df_sorted.FANCG.apply(int), df_sorted.NBN.apply(int), df_sorted.U2AF1.apply(int)], y=["ATM", "ATRX", "BRIP1", "FANCA", "BRCA1", "BRCA2", "FANCD2", "PTEN", "ATR", "CHEK2","FANCC", "FANCE", "FANCF", "FANCG", "NBN", "U2AF1"],x=df_sorted.Tumor_Sample_Barcode, showscale=False, showlegend=False, yaxis = dict(ticks='',nticks=16))
        fig = tls.make_subplots(rows=3, cols=1, specs=[[{}],[{'rowspan' : 2}], [None]])
        fig.append_trace(trace1, 1, 1)
        fig.append_trace(trace2, 2, 1)
        fig['layout']['xaxis1'].update(showticklabels = False)
        fig['layout']['xaxis2'].update(showticklabels = False)
        plot(fig, auto_open=False, filename=os.path.join(output_dir, "%s_heatmap.html" % cancer))