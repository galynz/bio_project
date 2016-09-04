# -*- coding: utf-8 -*-
"""
Created on Sun Jul 24 13:33:35 2016

@author: gal
"""

import csv, datetime, sys, glob
import xml.etree.ElementTree as ET
import logging
import logging.handlers
from optparse import OptionParser
import gzip
from tqdm import tqdm
from Bio import Entrez

import matplotlib
matplotlib.use('Agg')
from plotly.offline import download_plotlyjs, init_notebook_mode,  plot
import plotly.graph_objs as go
import numpy as np
import pandas as pd
import lifelines as ll

# Plotting helpers
#from IPython.display import HTML
import matplotlib.pyplot as plt
import plotly.plotly as py
import plotly.tools as tls   
from plotly.graph_objs import *
from lifelines.statistics import logrank_test

from pylab import rcParams
rcParams['figure.figsize']=10, 5

#init_notebook_mode()

HR_DEFICIENT_GENES = ("ATM", "ATRX", "BRIP1", "CHEK2", "FANCA", "FANCC", "BRCA1","BRCA2",
                      "FANCD2", "FANCE", "FANCF","FANCG", "NBN", "PTEN", "U2AF1")

NER_DEFICIENT_GENES = ("CCNH", "CDK7", "CENT2",  "DDB1", "DDB2",
                       "ERCC1", "ERCC2","ERCC3", "ERCC4", "ERCC5", 
                       "ERCC6", "ERCC8", "LIG1", "MNAT1", "MMS19",
                       "RAD23A", "RAD23B","RPA1", "RPA2", "TFIIH",
                       "XAB2", "XPA", "XPC")
                       
MMR_DEFICIENT_GENES = ("MLH1","MLH3","MSH2","MSH3","MSH6","PMS1","PMS2")

TOP_PERCENTIL = 4 #25%

HOT_SPOT_TRESHOLD = 3

logger = logging.getLogger("mutations_summary")

class Sample(object):
    def __init__(self, patient_barcode, tumor_barcode, norm_barcode):
        self.patient_barcode = patient_barcode
        self.tumor_barcode = tumor_barcode
        self.norm_barcode = norm_barcode
        self.centers = set()
        self.mutations = {}
        self.brca1 = {}
        self.brca2 = {}
        self.hr_deficient = {}
        self.ner_deficient = {}
        self.mmr_deficient = {}
        self.survival_days = 0
        self.survival_update = None
        self.dead = False
        self.clinical_available = False
        self.top_mutation_load = False
        self.low_mutation_load = False
        self.hot_spots = 0
        self.drugs_dict = {}
        logger.debug("added sample %s", patient_barcode)
        
    def add_mutation(self,hugo_symbol, mutation_type, mutation_pos):
        self.mutations.setdefault(mutation_type, {}).setdefault(hugo_symbol, set())
        if not self.mutations[mutation_type][hugo_symbol].intersection(set(mutation_pos)):
            self.mutations[mutation_type][hugo_symbol].add(mutation_pos)
            count = sum([len(i.get(hugo_symbol, set())) for i in self.mutations.values()])
            if count == HOT_SPOT_TRESHOLD:
                self.hot_spots+=1
            if hugo_symbol == "BRCA1":
                self.brca1[mutation_type] = self.brca1.get(mutation_type, 0) + 1
            if hugo_symbol == "BRCA2":
                self.brca2[mutation_type] = self.brca2.get(mutation_type, 0) + 1
            if hugo_symbol in HR_DEFICIENT_GENES:
                self.hr_deficient[mutation_type] = self.hr_deficient.get(mutation_type, 0) + 1
            if hugo_symbol in NER_DEFICIENT_GENES:
                self.ner_deficient[mutation_type] = self.ner_deficient.get(mutation_type, 0) + 1
            if hugo_symbol in MMR_DEFICIENT_GENES:
                self.mmr_deficient[mutation_type] = self.mmr_deficient.get(mutation_type, 0) + 1
        
    def count_mutations(self, distinct=True, mutation_type=None):
        if distinct:
            if mutation_type:
                return sum([len(self.mutations.get(mut_type, {})) for mut_type in mutation_type])
            # return the number of genes that has a mutation
            return sum([len(i) for i in self.mutations.values()])
        if mutation_type:
            #  count only specific mutation types
            return sum([sum([len(j) for j in self.mutations.get(i, {}).values()]) for i in mutation_type])
        # count all the mutations, regardless of mutations type
        return sum([sum([len(j) for j in i.values()]) for i in self.mutations.values()])
        
    def add_center(self, center):
        self.centers.add(center)
        
    def get_gene_mutations(self, hugo_symbol, distinct=True, mutation_type=[]):
        if distinct:
            if mutation_type:
                for i in mutation_type:
                    if self.mutations.get(i, {}).get(hugo_symbol, set()):
                        return 1
                return 0
            if sum([len(i.get(hugo_symbol, set())) for i in self.mutations.values()]):
                return 1
            return 0
        if mutation_type:
            return sum([len(self.mutations.get(i, {}).get(hugo_symbol, set())) for i in mutation_type])
        return sum([len(i.get(hugo_symbol, set())) for i in self.mutations.values()])
        
    def update_survival(self, survival_days, update_date, dead=False):
        if ((not self.survival_update) or update_date > self.survival_update) and survival_days:
            self.survival_days = survival_days            
            self.survival_update = update_date
            self.clinical_available = True
            self.dead = dead
            logger.debug("updated patient %s survival days to %s", self.patient_barcode, survival_days)
            
    def update_top_mutation_load(self):
        self.top_mutation_load = True
        logger.info("patient %s is in top mutation load group", self.patient_barcode)
        
    def update_low_mutation_load(self):
        self.low_mutation_load = True
        logger.info("patient %s in in low mutation load group", self.patient_barcode)
        
    def get_group(self):
        #TODO: add mutation type option
        group = "HR_PROFIECIENT"
        if self.brca1:
            group = "BRCA1_MUTATED"
        elif self.brca2:
            group = "BRCA2_MUTATED"
        elif self.hr_deficient:
            group = "HR_DEFICIENT"
        elif self.ner_deficient:
            group = "NER_DEFICIENT"
        elif self.mmr_deficient:
            group =  "MMR_DEFICIENT"
        logger.debug("patient %s group: %s", self.patient_barcode, group)
        return group
        
    def check_group_deficient(self, group, mutation_type):
        has_mutation = getattr(self, group)                
        if has_mutation:                    
            if mutation_type:
            # making sure that the mutation is of type we want to consider
                has_mutation_type = False
                for mut_type in mutation_type:
                    if has_mutation.has_key(mut_type):
                        has_mutation_type = True
                        break
                if has_mutation_type:
                    return True
                else:
                    # if the mutations in this group are of wrong types, return False
                    return False
            else:
                # if we want to consider all mutation types
                return True
        else:
            return False
            
    def update_drug_data(self, drugs_dict):
        self.drugs_dict = drugs_dict
        
    def sort_drugs_chronologically(self):
        self.drugs_list_sorted = [[]] * len(self.drugs_dict)
        counter = 0
        for drug_dict in self.drugs_dict.values():
            if drug_dict['regimen_number'] and drug_dict['regimen_number'].isdigit():
                index = int(drug_dict['regimen_number']) - 1
                self.drugs_list_sorted[index].append(drug_dict)
                if index > counter:
                    counter = index+1
        self.drugs_list_sorted = self.drugs_list_sorted[:counter]
                
    def get_last_drug_name(self):
        self.sort_drugs_chronologically()
        return [i.get('drug_name', 'None') for i in self.drugs_list_sorted[-1]]

        
class Mutation(object):
    def __init__(self, hugo_code):
        self.hugo_code = hugo_code
        self.samples = set()
        logger.debug("added mutation %s", self.hugo_code)
    
    def add_sample(self, sample):
        self.samples.add(sample)
        
    def count_top_mutation_load(self, distinct=True, mutation_type=[]):
        count = 0
        for sample in self.samples:
            if sample.top_mutation_load:
                count+=sample.get_gene_mutations(self.hugo_code, distinct, mutation_type)
        logger.debug("mutation %s has %d (distinct=%s) mutations in top group (mutation_type %s)", self.hugo_code, count, distinct, mutation_type)
        return count
        
    def count_low_mutation_load(self, distinct=True, mutation_type=[]):
        count = 0
        for sample in self.samples:
            if sample.low_mutation_load:
                count+=sample.get_gene_mutations(self.hugo_code, distinct, mutation_type)
        logger.debug("mutation %s has %d (distinct=%s) mutations in low group (mutation_type %s)", self.hugo_code, count, distinct, mutation_type)
        return count
                
    def count_all_mutation_load(self, distinct=True, mutation_type=[]):
        return sum([i.get_gene_mutations(self.hugo_code, distinct, mutation_type) for i in self.samples])
        
    def count_non_top_mutation_load(self, distinct=True,  mutation_type=[]):
        return (self.count_all_mutation_load(distinct, mutation_type)-self.count_top_mutation_load(distinct, mutation_type))
        
    def calc_top_mutation_load_patients_ratio(self, top_mutation_load_patients_num, distinct=True, mutation_type=[]):
        return float(self.count_top_mutation_load(distinct,mutation_type))/top_mutation_load_patients_num
        
    def calc_low_mutation_load_patients_ratio(self, low_mutation_load_patients_num, distinct=True,mutation_type=[]):
        return float(self.count_low_mutation_load(distinct,mutation_type))/low_mutation_load_patients_num
        
    def calc_non_top_mutation_load_patients_ratio(self, non_top_mutation_load_patients_num, distinct=True,mutation_type=[]):
        return float(self.count_non_top_mutation_load(distinct,mutation_type))/non_top_mutation_load_patients_num
        
    def calc_patients_ratio(self, patients_num, distinct=True,mutation_type=[]):
        return float(self.count_all_mutation_load(distinct,mutation_type))/patients_num
            
    def calc_top_mutation_load_percent(self, mutations_num, distinct=True, mutation_type=[]):
        return float(self.count_top_mutation_load(distinct,mutation_type))/mutations_num
        
    def calc_low_mutation_load_percent(self, mutations_num, distinct=True,mutation_type=[]):
        return float(self.count_low_mutation_load(distinct,mutation_type))/mutations_num
        
    def calc_non_top_mutation_load_percent(self, mutations_num, distinct=True,mutation_type=[]):
        return float(self.count_non_top_mutation_load(distinct,mutation_type))/mutations_num                 
    
class MutationsSummary(object):
    def __init__(self, csv_paths, clinical_paths):
        self.ids_dict = {}
        self.mutations_dict = {}
        for path in csv_paths:
            self.add_csv_data(path)
            logger.info("added %s to csv files", path)
        for path in clinical_paths:            
            self.add_clinical_data(path)
            logger.info("added %s to clinical files", path)
        
    def add_csv_data(self, path):
        with gzip.open(path) as f:
            line = f.readline().lower()
            while line.startswith('hugo')==False:
                #ignoring comment lines in the begining of the file
                logger.debug("ignoring line in the begining of the file: %s", line)
                line = f.readline().lower()
            file_dict = csv.DictReader(f, dialect=csv.excel_tab, fieldnames=line.split())
            for row in file_dict:
                try:
                    tumor_barcode = row["tumor_sample_barcode"]
                    if tumor_barcode.startswith('TCGA')!=True:
                        raise Exception("Invalid file, tumor barcode doesn't start with TCGA")
                except Exception, e:
                    print row
                    print path
                    print e
                    logger.exception(e)
                    sys.exit()
                norm_barcode = row["matched_norm_sample_barcode"]
                patient_barcode = "-".join(tumor_barcode.split('-')[:3])
                mutation = row["hugo_symbol"]
                mutation_type = row["variant_classification"]
                mutation_pos = row["start_position"]
                #center = row["Center"]
                sample = self.ids_dict.setdefault(patient_barcode, Sample(patient_barcode, tumor_barcode, norm_barcode))
                #sample.add_center(center)
                sample.add_mutation(mutation, mutation_type, mutation_pos)
                
                mutation_obj = self.mutations_dict.setdefault(mutation, Mutation(mutation))
                mutation_obj.add_sample(sample)
                
    def add_clinical_data(self, path):
        tree = ET.parse(path)
        patient_barcode = tree.findtext(".//{http://tcga.nci/bcr/xml/shared/2.7}bcr_patient_barcode")
        days_to_last_followup = tree.findtext(".//{http://tcga.nci/bcr/xml/clinical/shared/2.7}days_to_last_followup")
        days_to_death = tree.findtext(".//{http://tcga.nci/bcr/xml/clinical/shared/2.7}days_to_death")
        update_day = int(tree.findtext(".//{http://tcga.nci/bcr/xml/administration/2.7}day_of_dcc_upload"))
        update_month = int(tree.findtext(".//{http://tcga.nci/bcr/xml/administration/2.7}month_of_dcc_upload"))
        update_year = int(tree.findtext(".//{http://tcga.nci/bcr/xml/administration/2.7}year_of_dcc_upload"))
        update_date = datetime.date(update_year, update_month, update_day)
        sample = self.ids_dict.get(patient_barcode, None)
        if sample:
            drugs = tree.findall('.//{http://tcga.nci/bcr/xml/clinical/pharmaceutical/2.7}drug')
            drugs_dict = {}
            for drug in drugs:
                drug_name = drug.findtext('.//{http://tcga.nci/bcr/xml/clinical/pharmaceutical/2.7}drug_name') 
                drugs_dict[drug_name] = {}
                for child in list(drug):
                    key = child.tag.split("}")[-1]
                    value = child.text
                    if list(child):
                        value = []
                        for mini_child in list(child):
                            value.append((child, child.text))
                    elif value:
                        value = value.strip()
                    drugs_dict[drug_name][key] = value
            sample.update_drug_data(drugs_dict)
            if days_to_last_followup:
                logger.debug("updating patient %s survival according to days_to_last_followup", patient_barcode)
                sample.update_survival(days_to_last_followup, update_date)
            else:
                logger.debug("updating patient %s survival according to days_to_death", patient_barcode)
                sample.update_survival(days_to_death, update_date, True)
                
    def write_output(self, output_path, cancer, mutation_type):
        logger.info("writing output file %s", output_path)
        with open(output_path, "wb") as f:
            csv_file = csv.DictWriter(f, fieldnames=["Tumor_Sample_Barcode",
                                                     "Matched_Norm_Sample_Barcode", 
                                                     "Group", "Mutations_Count", 
                                                     "Mutations_Count_distinct",
                                                     "Cancer_Site", "Survival_days",
                                                     "BRCA1_mutated", "BRCA2_mutated",
                                                     "HR_mutated", "NER_mutated", "MMR_mutated"])
            csv_file.writeheader()
            for sample in self.ids_dict.values():
                group = sample.get_group()
                if mutation_type:
                    row_dict = {"Tumor_Sample_Barcode" : sample.tumor_barcode,
                                "Matched_Norm_Sample_Barcode" : sample.norm_barcode,
                                "Group" :  group,
                                "Mutations_Count" : sample.count_mutations(False, mutation_type),
                                "Mutations_Count_distinct" : sample.count_mutations(mutation_type=mutation_type),
                                "Cancer_Site" : cancer,
                                "Survival_days" : sample.survival_days,
                                "BRCA1_mutated" : sum([sample.brca1.get(i, 0) for i in mutation_type]),
                                "BRCA2_mutated" : sum([sample.brca2.get(i, 0) for i in mutation_type]),
                                "HR_mutated" : sum([sample.hr_deficient.get(i, 0) for i in mutation_type]),
                                "NER_mutated" : sum([sample.ner_deficient.get(i, 0) for i in mutation_type]),
                                "MMR_mutated" : sum([sample.mmr_deficient.get(i, 0) for i in mutation_type])}
                else:
                    row_dict = {"Tumor_Sample_Barcode" : sample.tumor_barcode,
                                "Matched_Norm_Sample_Barcode" : sample.norm_barcode,
                                "Group" :  group,
                                "Mutations_Count" : sample.count_mutations(False, mutation_type),
                                "Mutations_Count_distinct" : sample.count_mutations(mutation_type=mutation_type),
                                "Cancer_Site" : cancer,
                                "Survival_days" : sample.survival_days,
                                "BRCA1_mutated" : sum(sample.brca1.values()),
                                "BRCA2_mutated" : sum(sample.brca2.values()),
                                "HR_mutated" : sum(sample.hr_deficient.values()),
                                "NER_mutated" : sum(sample.ner_deficient.values()),
                                "MMR_mutated" : sum(sample.mmr_deficient.values())}
                csv_file.writerow(row_dict)
                    
    def write_mutation_load_output(self, output_path,cancer, mutation_type):
        self.find_high_low_mutation_load_patients(mutation_type)
        logger.info("writing mutation load output file %s", output_path)
        with open(output_path, "wb") as f:
            csv_file = csv.DictWriter(f, fieldnames=["Hugo_code", "samples_count",
                                                     "top_mutation_load_samples_count", 
                                                     "low_mutation_load_samples_count", 
                                                     "total_patients_ratio", 
                                                     "top_mutation_load_patients_ratio", 
                                                     "low_mutation_load_patients_ratio", 
                                                     "non_top_mutation_load_patients_ratio", 
                                                     "top_patients_mutations_ratio", 
                                                     "low_patients_mutations_ratio", 
                                                     "non_top_patients_mutations_ratio",
                                                     "samples_count_distinct", 
                                                     "top_mutation_load_samples_count_distinct", 
                                                     "low_mutation_load_samples_count_distinct", 
                                                     "total_patients_ratio_distinct", 
                                                     "top_mutation_load_patients_ratio_distinct", 
                                                     "low_mutation_load_patients_ratio_distinct", 
                                                     "non_top_mutation_load_patients_ratio_distinct", 
                                                     "top_patients_mutations_ratio_distinct", 
                                                     "low_patients_mutations_ratio_distinct", 
                                                     "non_top_patients_mutations_ratio_distinct",
                                                     "cancer"])
            csv_file.writeheader()
            patients_num = len(self.ids_dict.keys())
            top_mutation_load_num = len([i for i in self.ids_dict.values() if i.top_mutation_load])
            low_mutation_load_num = len([i for i in self.ids_dict.values() if i.low_mutation_load])
            top_patients_mutations_sum = sum([i.count_mutations(mutation_type) for i in self.ids_dict.values() if i.top_mutation_load])
            low_patients_mutations_sum = sum([i.count_mutations(mutation_type) for i in self.ids_dict.values() if i.low_mutation_load])
            patients_mutations_sum = sum([i.count_mutations(mutation_type) for i in self.ids_dict.values()])
            non_top_patients_mutations_sum = patients_mutations_sum - top_patients_mutations_sum
            for mut in self.mutations_dict.values():
                row_dict = {"Hugo_code": mut.hugo_code,
                           "samples_count_distinct" : mut.count_all_mutation_load(),
                           "top_mutation_load_samples_count_distinct" : mut.count_top_mutation_load(mutation_type=mutation_type),
                           "low_mutation_load_samples_count_distinct" : mut.count_low_mutation_load(mutation_type=mutation_type),
                           "top_mutation_load_patients_ratio_distinct" : mut.calc_top_mutation_load_patients_ratio(top_mutation_load_num,mutation_type=mutation_type),
                           "low_mutation_load_patients_ratio_distinct" : mut.calc_low_mutation_load_patients_ratio(low_mutation_load_num, mutation_type=mutation_type),
                           "total_patients_ratio_distinct" : mut.calc_patients_ratio(patients_num, mutation_type=mutation_type),
                           "non_top_mutation_load_patients_ratio_distinct" : mut.calc_non_top_mutation_load_patients_ratio(patients_num - top_mutation_load_num, mutation_type=mutation_type),
                           "top_patients_mutations_ratio_distinct" : mut.calc_top_mutation_load_percent(top_patients_mutations_sum, mutation_type=mutation_type),
                           "low_patients_mutations_ratio_distinct" : mut.calc_low_mutation_load_percent(low_patients_mutations_sum, mutation_type=mutation_type),
                           "non_top_patients_mutations_ratio_distinct" : mut.calc_non_top_mutation_load_percent(non_top_patients_mutations_sum, mutation_type=mutation_type),
                           "samples_count" : mut.count_all_mutation_load(False, mutation_type=mutation_type),
                           "top_mutation_load_samples_count" : mut.count_top_mutation_load(False, mutation_type=mutation_type),
                           "low_mutation_load_samples_count" : mut.count_low_mutation_load(False, mutation_type=mutation_type),
                           "top_mutation_load_patients_ratio" : mut.calc_top_mutation_load_patients_ratio(top_mutation_load_num, False, mutation_type=mutation_type),
                           "low_mutation_load_patients_ratio" : mut.calc_low_mutation_load_patients_ratio(low_mutation_load_num, False, mutation_type=mutation_type),
                           "total_patients_ratio" : mut.calc_patients_ratio(patients_num,False,mutation_type=mutation_type),
                           "non_top_mutation_load_patients_ratio" : mut.calc_non_top_mutation_load_patients_ratio(patients_num - top_mutation_load_num,False, mutation_type=mutation_type),
                           "top_patients_mutations_ratio" : mut.calc_top_mutation_load_percent(top_patients_mutations_sum,False, mutation_type=mutation_type),
                           "low_patients_mutations_ratio" : mut.calc_low_mutation_load_percent(low_patients_mutations_sum,False, mutation_type=mutation_type),
                           "non_top_patients_mutations_ratio" : mut.calc_non_top_mutation_load_percent(non_top_patients_mutations_sum,False, mutation_type=mutation_type),
                           "cancer" : cancer}
                csv_file.writerow(row_dict)
                    
    def find_high_low_mutation_load_patients(self, mutation_type=[]):
        patients_list = sorted(self.ids_dict.values(), key=lambda x: x.count_mutations(False, mutation_type), reverse=True)
        stop_index = int(len(patients_list)/TOP_PERCENTIL)
        for patient in patients_list[:stop_index]:
            patient.update_top_mutation_load()
        for patient in patients_list[-stop_index:]:
            patient.update_low_mutation_load()
                
                
#    def write_survival_output(self, output_path, cancer):
#        logger.info("writing survival output file %s", output_path)
#        with open(output_path, "w") as f:
#            csv_file = csv.DictWriter(f,  fieldnames = ["Days", "Group", "Num", "percent out of group","Cancer"])
#            csv_file.writeheader()
#            days_dict = {}
#            group_count = {}
#            for sample in self.ids_dict.values():
#                group = sample.get_group()
#                group_list = days_dict.setdefault(group, [])
#                if sample.clinical_available:
#                    group_count[group] = group_count.get(group, 0) + 1
#                    for i in range(0,int(sample.survival_days)):
#                        if i >= len(group_list):
#                            group_list += [1]
#                        else:
#                            group_list[i] += 1
#            for group, group_list in days_dict.items():
#                for days,count in enumerate(group_list):
#                    csv_file.writerow({"Days" : days, "Group" : group, "Num" : count, "Cancer" : cancer, "percent out of group" : count/group_count[group]})
                    
    def plot_mutation_load_box(self, output_path, cancer, distinct, mutation_type):
        logger.info("plotting mutation load box plots and saving it to %s", output_path)
        count_dict ={}
        groups = ('brca1', 'brca2', 'hr_deficient', 'ner_deficient', 'mmr_deficient')
        for group in groups:
            count_dict[group] = {'deficient' : [], 'proficient' : []}
        for sample in self.ids_dict.values():
            count = sample.count_mutations(distinct, mutation_type)
            for group in groups:
                if sample.check_group_deficient(group, mutation_type):
                    count_dict[group]['deficient'].append(count)                        
                else:
                    # the patient has no mutations in this gene/pathway
                    count_dict[group]['proficient'].append(count)
        x_deficient = []
        x_proficient = []
        
        y_deficient = []
        y_proficient = []
        
        for group in groups:
#            print group, len(count_dict[group]['deficient'])
            pvalue = tls.scipy.stats.ttest_ind(count_dict[group]['deficient'], count_dict[group]['proficient'], equal_var=False).pvalue
            ks_pvalue = tls.scipy.stats.ks_2samp(count_dict[group]['deficient'], count_dict[group]['proficient']).pvalue
            x_deficient.extend(['%s<br>t_test-pvalue=%s<br>ks-pvalue=%s' % (group, pvalue, ks_pvalue)] * len(count_dict[group]['deficient']))
            y_deficient.extend(count_dict[group]['deficient'])
            x_proficient.extend(['%s<br>t_test-pvalue=%s<br>ks-pvalue=%s' % (group, pvalue, ks_pvalue)] * len(count_dict[group]['proficient']))
            y_proficient.extend(count_dict[group]['proficient'])
        deficient = go.Box(y=y_deficient, x=x_deficient, 
                           name='deficient', marker=dict(color='#3D9970'))
        proficient = go.Box(y=y_proficient, x=x_proficient, 
                           name='proficient', marker=dict(color='#FF4136'))
        data = [deficient, proficient]
        
        layout = go.Layout(yaxis=dict(title='%s mutation load by group (deficient/proficient)'% cancer,
                                      zeroline=False),
                            boxmode='group')
        fig = go.Figure(data=data, layout=layout)
        plot(fig, filename=output_path, auto_open=False)
     
     
    def plot_survival(self, output_path, cancer, mutation_type=[]):
        self.find_high_low_mutation_load_patients()
        l = [(i.patient_barcode, int(i.survival_days), i.dead, i.get_group(), i.check_group_deficient('brca1', mutation_type), i.check_group_deficient('brca2', mutation_type), i.check_group_deficient('hr_deficient', mutation_type), i.check_group_deficient('ner_deficient', mutation_type), i.check_group_deficient('mmr_deficient', mutation_type), i.top_mutation_load, i.low_mutation_load) for i in self.ids_dict.values() if i.clinical_available]
        df = pd.DataFrame(data=l, columns=["patient_barcode", "days", "dead", "group", 'BRCA1', 'BRCA2', 'HR', 'NER', 'MMR','top_mutation_load','low_mutation_load'])
#        groups = ('BRCA1', 'BRCA2', 'HR', 'NER', 'MMR')
        T = df['days']
        C = df['dead']
#        kmf = ll.KaplanMeierFitter()
#        for i, group in enumerate(groups):            
#            ax = plt.subplot(2,3,i+1)
#            ix = (df[group] == True)
#            if len(df[ix]) and len(df[~ix]):
#                kmf.fit(T[~ix], C[~ix], label='%s_proficient' % group)
#                kmf.plot(ax=ax)
#                kmf.fit(T[ix], C[ix], label='%s_deficient' % group)
#                kmf.plot(ax=ax, legend=True)
#                ax.set_title('%s survival - %s'% (group, cancer))
#            else:
#                logger.info("cancer %s has no %s deficient", cancer, group)
##            print group
##            print logrank_test(T[ix], T[~ix], C[ix], C[~ix], alpha=0.95)
#        plt.tight_layout()
#        kmf1 = plt.gcf()
#        pyplot(kmf1, output_path + '%s.html' % cancer, ci=False)
#            
#        kmf = ll.KaplanMeierFitter()
#        groups = ('BRCA1_MUTATED', 'BRCA2_MUTATED', 'HR_DEFICIENT', 'NER_DEFICIENT', 'MMR_DEFICIENT','HR_PROFIECIENT')
#        ax = plt.subplot(111)
#        for group in groups:
#            ix = (df['group'] == group)
#            if len(df[ix]) and len(df[~ix]):
#                kmf.fit(T[ix], C[ix], label=group)
#                kmf.survival_function_.plot(ax=ax)
#            #kmf.survival_function_.plot()
#        plt.title("survival by groups - %s" % cancer)
#        kmf2 = plt.gcf()
#        pyplot(kmf2, output_path + '.%s.groups_compare.html'% cancer, ci=False)
    
            
        kmf = ll.KaplanMeierFitter()
        ax = plt.subplot(111)
        ix = (df['top_mutation_load'] == True)
        kmf.fit(T[ix], C[ix], label='top_mutation_load_patients')
        kmf.survival_function_.plot(ax=ax)
        ix2 = (df['low_mutation_load'] == True)
        kmf.fit(T[ix2], C[ix2], label='low_mutation_load_patients')
        kmf.survival_function_.plot(ax=ax)
        kmf.fit(T, C, label='all_patients')
        kmf.survival_function_.plot(ax=ax)
        plt.title('top/low mutation load patients survival - %s'% cancer)
        kmf3 = plt.gcf()
        pyplot(kmf3, output_path + '.top_low_mutation_load_patietns.%s.html' % (cancer), ci=False)
        print "cancer:", cancer
        print 'top/low mutation load patients'
        print logrank_test(T[ix], T[ix2], C[ix], C[ix2], alpha=0.99)
        print 'top:', len(T[ix]), 'low:', len(T[ix2])
        print 'top/non-top mutation load patients'
        print logrank_test(T[ix], T[~ix], C[ix], C[~ix], alpha=0.99)
        print 'top:', len(T[ix]), 'non-top:', len(T[~ix])
        
        
    def plot_hot_spot_box(self, output_path, cancer, plot_type='box', mutation_type=[]):
        logger.info("plotting hot spots box plots and saving it to %s", output_path)
        count_dict ={}
        groups = ('brca1', 'brca2', 'hr_deficient', 'ner_deficient', 'mmr_deficient')
        for group in groups:
            count_dict[group] = {'deficient' : [], 'proficient' : []}
        for sample in self.ids_dict.values():
            count = sample.hot_spots
            for group in groups:
                if sample.check_group_deficient(group, mutation_type):
                    count_dict[group]['deficient'].append(count)                        
                else:
                    # the patient has no mutations in this gene/pathway
                    count_dict[group]['proficient'].append(count)
        x_deficient = []
        x_proficient = []
        
        y_deficient = []
        y_proficient = []
        
        for group in groups:
            pvalue = tls.scipy.stats.ttest_ind(count_dict[group]['deficient'], count_dict[group]['proficient'], equal_var=False).pvalue
            ks_pvalue = tls.scipy.stats.ks_2samp(count_dict[group]['deficient'], count_dict[group]['proficient']).pvalue
            x_deficient.extend(['%s<br>t_test-pvalue=%s<br>ks-pvalue=%s' % (group, pvalue, ks_pvalue)] * len(count_dict[group]['deficient']))
            y_deficient.extend(count_dict[group]['deficient'])
            x_proficient.extend(['%s<br>t_test-pvalue=%s<br>ks-pvalue=%s' % (group, pvalue, ks_pvalue)] * len(count_dict[group]['proficient']))
            y_proficient.extend(count_dict[group]['proficient'])
        if plot_type == 'bar':
            ### NOT recommended
            deficient = go.Bar(y=y_deficient, x=x_deficient, 
                           name='deficient', marker=dict(color='#199219'))
            proficient = go.Bar(y=y_proficient, x=x_proficient, 
                           name='proficient', marker=dict(color='#78c578'))
            data = [deficient, proficient]
            layout = go.Layout(yaxis=dict(title='%s hot spots load by group (deficient/proficient)'% cancer,
                                      zeroline=False),
                            barmode='group')
        elif plot_type == 'box':
            deficient = go.Box(y=y_deficient, x=x_deficient, 
                           name='deficient', marker=dict(color='#199219'))
            proficient = go.Box(y=y_proficient, x=x_proficient, 
                           name='proficient', marker=dict(color='#78c578'))
            data = [deficient, proficient]
            layout = go.Layout(yaxis=dict(title='%s hot spots load by group (deficient/proficient)'% cancer,
                                      zeroline=False),
                            boxmode='group')
        fig = go.Figure(data=data, layout=layout)
        plot(fig, filename=output_path, auto_open=False)
        
    def create_samples_mutations_dataframe(self, mutation_type):
        samples_list = sorted(self.ids_dict.keys())
        mutation_list = sorted(self.mutations_dict.keys())
        self.samples_mutations_dataframe = pd.DataFrame()
#        for mutation in tqdm(mutation_list):
        for mutation in mutation_list:
            mutation_samples = []
            for sample_id in samples_list:
                sample = self.ids_dict[sample_id]
                mutation_samples.append(sample.get_gene_mutations(mutation, False, mutation_type=mutation_type))
            self.samples_mutations_dataframe[mutation] = mutation_samples
        self.samples_mutations_dataframe['mutations_num'] = self.samples_mutations_dataframe.sum(axis=1, numeric_only=True)
        
        self.find_high_low_mutation_load_patients(mutation_type)
        self.samples_mutations_dataframe['top_mutation_load'] = [self.ids_dict[sample_id].top_mutation_load for sample_id in samples_list]
        self.samples_mutations_dataframe['low_mutation_load'] = [self.ids_dict[sample_id].low_mutation_load for sample_id in samples_list]
        
    def calc_mutation_ratio(self):
        samples_mutation_count = self.samples_mutations_dataframe.sum(axis=1, numeric_only=True)
        self.samples_mutations_ratio_dataframe = pd.DataFrame()
        for column in self.samples_mutations_dataframe.columns[:-3]:
            self.samples_mutations_ratio_dataframe[column] = self.samples_mutations_dataframe[column]/samples_mutation_count
        self.samples_mutations_ratio_dataframe['top_mutation_load'] = self.samples_mutations_dataframe['top_mutation_load']
        self.samples_mutations_ratio_dataframe['low_mutation_load'] = self.samples_mutations_dataframe['low_mutation_load']
        
    def calc_mutation_ratio_pvalue_top_non_top(self):
        top_mutation_load = (self.samples_mutations_ratio_dataframe['top_mutation_load'] == True)
        #low_mutation_load = (self.samples_mutations_ratio_dataframe['low_mutation_load'] == True)
        self.mutation_ratio_pvalue_top_low = pd.DataFrame()
        for mutation in self.samples_mutations_ratio_dataframe.columns[:-3]:
            self.mutation_ratio_pvalue_top_low[mutation] = [tls.scipy.stats.ttest_ind(self.samples_mutations_ratio_dataframe[top_mutation_load][mutation], self.samples_mutations_ratio_dataframe[~top_mutation_load][mutation]).pvalue]
            
    def create_mutations_ratio_pvalue_csv(self, mutation_type, output_path):
        self.create_samples_mutations_dataframe(mutation_type)
        self.calc_mutation_ratio()
        self.calc_mutation_ratio_pvalue_top_non_top()
        top_mutation_load = (self.samples_mutations_ratio_dataframe['top_mutation_load'] == True)
        with open(output_path, 'wb') as f:
            csv_writer = csv.DictWriter(f, ['Hugo_symbol', 'p-value', 'top_avg_ratio', 'non_top_avg_ratio'])
            csv_writer.writeheader()
            for mutation in self.mutation_ratio_pvalue_top_low.columns:
                csv_writer.writerow({'Hugo_symbol' : mutation,
                                     'p-value' : self.mutation_ratio_pvalue_top_low[mutation][0],
                                     'top_avg_ratio' : self.samples_mutations_ratio_dataframe[mutation][top_mutation_load].sum()/len(self.samples_mutations_ratio_dataframe[mutation][top_mutation_load]),
                                     'non_top_avg_ratio' : self.samples_mutations_ratio_dataframe[mutation][~top_mutation_load].sum()/len(self.samples_mutations_ratio_dataframe[mutation][~top_mutation_load])})
        
    def create_mutations_overload_pvalue_csv(self, output_path):
        Entrez.email = 'galynz@gmail.com'
        with open(output_path, 'wb') as f:
            csv_writer = csv.DictWriter(f, ['Hugo_symbol', 't-test_p-value',
                                            'ks-test_p-value',
                                            'deficient_num', 'proficient_num', 
                                            'deficient_mean', 'proficient_mean',
                                            'deficient_median', 'proficient_median',
                                            'deficient_std', 'proficient_std',
                                            'size', 'description'])
            csv_writer.writeheader()
#            for mutation in tqdm(self.samples_mutations_dataframe.columns[:-3]):
            for mutation in self.samples_mutations_dataframe.columns[:-3]:
                ix = self.samples_mutations_dataframe[mutation] > 0
                deficient = self.samples_mutations_dataframe[ix]['mutations_num']
                proficient = self.samples_mutations_dataframe[~ix]['mutations_num']
                
                #get gene info
                try:
                    handle = Entrez.esearch("gene", term="%s[Gene] AND Human[Orgn]" % mutation)
                    result = Entrez.read(handle)
                    mut_id = result['IdList'][0]
                    request = Entrez.epost('gene', id=mut_id)
                    result = Entrez.read(request)
                    webEnv = result["WebEnv"]
                    queryKey = result["QueryKey"]
                    data = Entrez.esummary(db="gene", webenv=webEnv, query_key = queryKey)
                    annotations = Entrez.read(data)
    
                    start = int(annotations['DocumentSummarySet']['DocumentSummary'][0]['GenomicInfo'][0]['ChrStart'])
                    end = int(annotations['DocumentSummarySet']['DocumentSummary'][0]['GenomicInfo'][0]['ChrStop'])                
                    size = abs(start-end)
                    desc = annotations['DocumentSummarySet']['DocumentSummary'][0]['Description']
                except Exception, e:
                    logger.exception(e)
                    logger.error("can't get info for Hugo symbol %s", mutation)
                    size = 'Failed'
                    desc = 'Failed'
                
                if len(deficient) > 20:
                    csv_writer.writerow({'Hugo_symbol' : mutation,
                                         't-test_p-value' : tls.scipy.stats.ttest_ind(deficient, proficient, equal_var=False).pvalue,
                                         'ks-test_p-value' : tls.scipy.stats.ks_2samp(deficient, proficient).pvalue,
                                         'deficient_num' : len(deficient),
                                         'proficient_num' : len(proficient),
                                         'deficient_mean' : tls.scipy.mean(deficient),
                                         'proficient_mean' : tls.scipy.mean(proficient),
                                         'deficient_median' : tls.scipy.median(deficient),
                                         'proficient_median' : tls.scipy.median(proficient),
                                         'deficient_std' : tls.scipy.std(deficient),
                                         'proficient_std' : tls.scipy.std(proficient),
                                         'size' : size, 'description' : desc})

def pyplot(fig, output_path, ci=False, legend=True):
    # Convert mpl fig obj to plotly fig obj, resize to plotly's default
    py_fig = tls.mpl_to_plotly(fig, resize=True)
    
    # Add fill property to lower limit line
    if ci == True:
        style1 = dict(fill='tonexty')
        # apply style
        py_fig['data'][2].update(style1)
        
        # Change color scheme to black
        py_fig['data'].update(dict(line=Line(color='black')))
    
    # change the default line type to 'step'
    py_fig['data'].update(dict(line=Line(shape='hv')))
    # Delete misplaced legend annotations 
    py_fig['layout'].pop('annotations', None)
    
    if legend == True:
        # Add legend, place it at the top right corner of the plot
        py_fig['layout'].update(
            showlegend=True,
            legend=Legend(
                x=1.05,
                y=1
            )
        )
        
    # Send updated figure object to Plotly, show result in notebook
    return plot(py_fig, filename=output_path, auto_open=False)
        
            
def main():
    parser = OptionParser()
    parser.add_option("-c", "--cancer", dest="cancer", help="tumor main site")
    parser.add_option("-p", "--csv_path", dest="csv_path", help="csv paths, use '' if the path contains *")
    parser.add_option("--clinical_path", dest="clinical_path", help="clinical paths, use '' if the path contains *")
    parser.add_option("--debug", default=False, action="store_false", dest="debug", help="run the script in debug mode")
    parser.add_option("-m", "--mutation_types", dest="mutation_types", default="", help="mutation types (silent, missense...) to report, a string")
    parser.add_option("-o", "--output", dest="output_path",default="output", help="where to save all the html files")
    
    (options, args) = parser.parse_args()
    
    # Logger
    formatter = logging.Formatter('%(process)d %(asctime)s %(levelname)s %(message)s')
    handler = logging.handlers.RotatingFileHandler(
              "mutations_summary.log", maxBytes=10*1024*1024, backupCount=5, mode="a")

    if options.debug:
        logging_level = logging.DEBUG
    else:
        logging_level = logging.INFO
    logger.setLevel(logging_level)
    handler.setLevel(logging_level)
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.info("starting the script")
    logger.info("cancer type: %s", options.cancer)
    
    if options.clinical_path:
        clinical_paths = glob.glob(options.clinical_path)
    else:
        clinical_paths = []
    if options.mutation_types:
        mutation_types = options.mutation_types.split(",")
    else:
        mutation_types = []
    summary = MutationsSummary(glob.glob(options.csv_path), clinical_paths)
    summary.write_mutation_load_output("%s.mutation_load.csv" % options.output_path, options.cancer, mutation_types)
    summary.write_output("%s.patients_summary.csv" % options.output_path, options.cancer, mutation_types)
#    summary.write_survival_output("survival_report_%s.csv" % options.cancer, options.cancer)
    summary.plot_mutation_load_box("%s.mutation_load" % options.output_path, options.cancer, False, mutation_types)
    summary.plot_hot_spot_box("%s.hot_spot" % options.output_path, options.cancer, mutation_type=mutation_types)
    summary.plot_survival("%s.survival" % options.output_path, options.cancer, mutation_types)
#    summary.create_samples_mutations_dataframe(mutation_types)
#    summary.create_mutations_overload_pvalue_csv('%s.mutation_load_per_mutation.csv' % options.output_path)
    
#if __name__ == "__main__":
#    main()
        
        
