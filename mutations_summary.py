# -*- coding: utf-8 -*-
"""
Created on Sun Jul 24 13:33:35 2016

@author: gal
"""

import csv, datetime
import xml.etree.ElementTree as ET

HR_DEFICIENT_GENES = ("ATM", "ATRX", "BRIP1", "CHEK2", "FANCA", "FANCC", "FANCD2", "FANCE", "FANCF",
"FANCG", "NBN", "PTEN", "U2AF1")

NER_DEFICIENT_GENES = ("CCNH", "CDK7", "CENT2",  "DDB1", "DDB2",
                       "ERCC1", "ERCC2","ERCC3", "ERCC4", "ERCC5", 
                       "ERCC6", "ERCC8", "LIG1", "MNAT1", "MMS19",
                       "RAD23A", "RAD23B","RPA1", "RPA2", "TFIIH",
                       "XAB2", "XPA", "XPC")
                       
MMR_DEFICIENT_GENES = ("MLH1","MLH3","MSH2","MSH3","MSH6","PMS1","PMS2")

TOP_PERCENTIL = 4 #25%

class Sample(object):
    def __init__(self, patient_barcode, tumor_barcode, norm_barcode):
        self.tumor_barcode = tumor_barcode
        self.norm_barcode = norm_barcode
        self.centers = set()
        self.mutations = set()
        self.brca1 = False
        self.brca2 = False
        self.hr_deficient = False
        self.ner_deficient = False
        self.mmr_deficient = False
        self.survival_days = 0
        self.survival_update = None
        self.clinical_available =False
        self.top_mutation_load = False
        
    def add_mutation(self,hugo_symbol):
        self.mutations.add(hugo_symbol)
        if hugo_symbol == "BRCA1":
            self.brca1 = True
        elif hugo_symbol == "BRCA2":
            self.brca2 = True
        elif hugo_symbol in HR_DEFICIENT_GENES:
            self.hr_deficient = True
        elif hugo_symbol in NER_DEFICIENT_GENES:
            self.ner_deficient = True
        elif hugo_symbol in MMR_DEFICIENT_GENES:
            self.mmr_deficient = True
        
    def count_mutations(self):
        return len(self.mutations)
        
    def add_center(self, center):
        self.centers.add(center)
        
    def update_survival(self, survival_days, update_date):
        if ((not self.survival_update) or update_date > self.survival_update) and survival_days:
            self.survival_days = survival_days            
            self.survival_update = update_date
            self.clinical_available = True
            
    def update_top_mutation_load(self):
        self.top_mutation_load = True
        
    def get_group(self):
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
        return group
        
class Mutation(object):
    def __init__(self, hugo_code):
        self.hugo_code = hugo_code
        self.samples = set()
    
    def add_sample(self, sample):
        self.samples.add(sample)
        
    def count_top_mutation_load(self):
        count = 0
        for sample in self.samples:
            if sample.top_mutation_load:
                count+=1
        return count
                
    def count_all_mutation_load(self):
        return len(self.samples)
        
    def count_lower_mutation_load(self):
        return (self.count_all_mutation_load()-self.count_top_mutation_load())
        
    def calc_top_mutation_load_patients_ratio(self, top_mutation_load_patients_num):
        return self.count_top_mutation_load()/top_mutation_load_patients_num
        
    def calc_non_top_mutation_load_patients_ratio(self, non_top_mutation_load_patients_num):
        return self.count_lower_mutation_load()/non_top_mutation_load_patients_num
        
    def calc_patients_ratio(self, patients_num):
        return self.count_all_mutation_load()/patients_num
            
    def calc_top_mutation_load_percent(self, mutations_num):
        return self.count_top_mutation_load()/mutations_num
        
    def calc_lower_mutation_load_percent(self, mutations_num):
        return self.count_lower_mutation_load()/mutations_num
    
class MutationsSummary(object):
    def __init__(self, csv_paths, clinical_paths):
        self.ids_dict = {}
        self.mutations_dict = {}
        for path in csv_paths:
            self.add_csv_data(path)
        for path in clinical_paths:
            self.add_clinical_data(path)
        
    def add_csv_data(self, path):
        with open(path) as f:
            file_dict = csv.DictReader(f, dialect=csv.excel_tab)
            for row in file_dict:
                tumor_barcode = row["Tumor_Sample_Barcode"]
                norm_barcode = row["Matched_Norm_Sample_Barcode"]
                patient_barcode = "-".join(tumor_barcode.split('-')[:3])
                mutation = row["Hugo_Symbol"]
                #center = row["Center"]
                sample = self.ids_dict.setdefault(patient_barcode, Sample(patient_barcode, tumor_barcode, norm_barcode))
                #sample.add_center(center)
                sample.add_mutation(mutation)
                
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
            if days_to_last_followup:
                sample.update_survival(days_to_last_followup, update_date)
            else:
                sample.update_survival(days_to_death, update_date)
        
                
    def write_output(self, output_path, cancer):
        with open(output_path, "w") as f:
            csv_file = csv.DictWriter(f, fieldnames=["Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode", "Group", "Mutations_Count", "Cancer_Site", "Survival_days"])
            csv_file.writeheader()
            for sample in self.ids_dict.values():
                group = sample.get_group()
                row_dict = {"Tumor_Sample_Barcode" : sample.tumor_barcode,
                            "Matched_Norm_Sample_Barcode" : sample.norm_barcode,
                            "Group" :  group,
                            "Mutations_Count" : sample.count_mutations(),
                            "Cancer_Site" : cancer,
                            "Survival_days" : sample.survival_days}
                csv_file.writerow(row_dict)
                    
    def write_mutation_load_output(self, output_path,cancer):
        self.find_high_mutation_load_patients()
        with open(output_path, "w") as f:
            csv_file = csv.DictWriter(f, fieldnames=["Hugo_code", "samples_count", "top_mutation_load_samples_count", "total_patients_ratio", "top_mutation_load_patients_ratio", "non_top_mutation_load_patients_ratio", "top_patients_mutations_ratio", "lower_patients_mutations_ratio","cancer"])
            csv_file.writeheader()
            patients_num = len(self.ids_dict.keys())
            top_mutation_load_num = len([i for i in self.ids_dict.values() if i.top_mutation_load])
            top_patients_mutations_sum = sum([i.count_mutations() for i in self.ids_dict.values() if i.top_mutation_load])
            patients_mutations_sum = sum([i.count_mutations() for i in self.ids_dict.values()])
            lower_patients_mutations_sum = patients_mutations_sum - top_patients_mutations_sum
            for mut in self.mutations_dict.values():
                row_dict = {"Hugo_code": mut.hugo_code,
                           "samples_count" : mut.count_all_mutation_load(),
                           "top_mutation_load_samples_count" : mut.count_top_mutation_load(),
                           "top_mutation_load_patients_ratio" : mut.calc_top_mutation_load_patients_ratio(top_mutation_load_num),
                           "total_patients_ratio" : mut.calc_patients_ratio(patients_num),
                            "non_top_mutation_load_patients_ratio" : mut.calc_non_top_mutation_load_patients_ratio(patients_num - top_mutation_load_num),
                            "top_patients_mutations_ratio" : mut.calc_top_mutation_load_percent(top_patients_mutations_sum),
                            "lower_patients_mutations_ratio" : mut.calc_lower_mutation_load_percent(lower_patients_mutations_sum),
                            "cancer" : cancer}
                csv_file.writerow(row_dict)
                    
    def find_high_mutation_load_patients(self):
        patients_list = sorted(self.ids_dict.values(), key=lambda x: x.count_mutations(), reverse=True)
        stop_index = int(len(patients_list)/TOP_PERCENTIL)
        for patient in patients_list[:stop_index]:
            patient.update_top_mutation_load()
                
                
    def write_survival_output(self, output_path, cancer):
        with open(output_path, "w") as f:
            csv_file = csv.DictWriter(f,  fieldnames = ["Days", "Group", "Num", "percent out of group","Cancer"])
            csv_file.writeheader()
            days_dict = {}
            group_count = {}
            for sample in self.ids_dict.values():
                group = sample.get_group()
                group_list = days_dict.setdefault(group, [])
                if sample.clinical_available:
                    group_count[group] = group_count.get(group, 0) + 1
                    for i in range(0,int(sample.survival_days)):
                        if i >= len(group_list):
                            group_list += [1]
                        else:
                            group_list[i] += 1
            for group, group_list in days_dict.items():
                for days,count in enumerate(group_list):
                    csv_file.writerow({"Days" : days, "Group" : group, "Num" : count, "Cancer" : cancer, "percent out of group" : count/group_count[group]})
        
            
                