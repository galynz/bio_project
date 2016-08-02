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

HR_DEFICIENT_GENES = ("ATM", "ATRX", "BRIP1", "CHEK2", "FANCA", "FANCC", "FANCD2", "FANCE", "FANCF",
"FANCG", "NBN", "PTEN", "U2AF1")

NER_DEFICIENT_GENES = ("CCNH", "CDK7", "CENT2",  "DDB1", "DDB2",
                       "ERCC1", "ERCC2","ERCC3", "ERCC4", "ERCC5", 
                       "ERCC6", "ERCC8", "LIG1", "MNAT1", "MMS19",
                       "RAD23A", "RAD23B","RPA1", "RPA2", "TFIIH",
                       "XAB2", "XPA", "XPC")
                       
MMR_DEFICIENT_GENES = ("MLH1","MLH3","MSH2","MSH3","MSH6","PMS1","PMS2")

TOP_PERCENTIL = 10 #10%

logger = logging.getLogger("mutations_summary")

class Sample(object):
    def __init__(self, patient_barcode, tumor_barcode, norm_barcode):
        self.patient_barcode = patient_barcode
        self.tumor_barcode = tumor_barcode
        self.norm_barcode = norm_barcode
        self.centers = set()
        self.mutations = {}
        self.brca1 = False
        self.brca2 = False
        self.hr_deficient = False
        self.ner_deficient = False
        self.mmr_deficient = False
        self.survival_days = 0
        self.survival_update = None
        self.clinical_available =False
        self.top_mutation_load = False
        self.least_mutation_load = False
        logger.debug("added sample %s", patient_barcode)
        
    def add_mutation(self,hugo_symbol):
        self.mutations[hugo_symbol] = self.mutations.get(hugo_symbol, 0) + 1
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
        
    def count_mutations(self, distinct=True):
        if distinct:
            return len(self.mutations)
        return sum(self.mutations.values())
        
    def add_center(self, center):
        self.centers.add(center)
        
    def get_gene_mutations(self, hugo_symbol):
        return self.mutations.get(hugo_symbol, 0)
        
    def update_survival(self, survival_days, update_date):
        if ((not self.survival_update) or update_date > self.survival_update) and survival_days:
            self.survival_days = survival_days            
            self.survival_update = update_date
            self.clinical_available = True
            logger.debug("updated patient %s survival days to %s", self.patient_barcode, survival_days)
            
    def update_top_mutation_load(self):
        self.top_mutation_load = True
        logger.info("patient %s is in top mutation load group", self.patient_barcode)
        
    def update_least_mutation_load(self):
        self.least_mutation_load = True
        logger.info("patient %s in in lower mutation load group", self.patient_barcode)
        
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
        logger.debug("patient %s group: %s", self.patient_barcode, group)
        return group
        
class Mutation(object):
    def __init__(self, hugo_code):
        self.hugo_code = hugo_code
        self.samples = set()
        logger.debug("added mutation %s", self.hugo_code)
    
    def add_sample(self, sample):
        self.samples.add(sample)
        
    def count_top_mutation_load(self, distinct=True):
        count = 0
        for sample in self.samples:
            if sample.top_mutation_load:
                if distinct:
                    count+=1
                else:
                    count+=sample.get_gene_mutations(self.hugo_code)
        logger.debug("mutation %s has %d (distinct=%s) mutations in top group", self.hugo_code, count, distinct)
        return count
        
    def count_lower_mutation_load(self, distinct=True):
        count = 0
        for sample in self.samples:
            if sample.lower_mutation_load:
                if distinct:
                    count+=1
                else:
                    count+=sample.get_gene_mutations(self.hugo_code)
        logger.debug("mutation %s has %d (distinct=%s) mutations in lower group", self.hugo_code, count, distinct)
        return count
                
    def count_all_mutation_load(self, distinct=True):
        if distinct:
            return len(self.samples)
        return sum([i.get_gene_mutations(self.hugo_code) for i in self.samples])
        
    def count_non_top_mutation_load(self, distinct=True):
        return (self.count_all_mutation_load(distinct)-self.count_top_mutation_load(distinct))
        
    def calc_top_mutation_load_patients_ratio(self, top_mutation_load_patients_num, distinct=True):
        return self.count_top_mutation_load(distinct)/top_mutation_load_patients_num
        
    def calc_lower_mutation_load_patients_ratio(self, lower_mutation_load_patients_num, distinct=True):
        return self.count_lower_mutation_load(distinct)/lower_mutation_load_patients_num
        
    def calc_non_top_mutation_load_patients_ratio(self, non_top_mutation_load_patients_num, distinct=True):
        return self.count_non_top_mutation_load(distinct)/non_top_mutation_load_patients_num
        
    def calc_patients_ratio(self, patients_num, distinct=True):
        return self.count_all_mutation_load(distinct)/patients_num
            
    def calc_top_mutation_load_percent(self, mutations_num, distinct=True):
        return self.count_top_mutation_load(distinct)/mutations_num
        
    def calc_lower_mutation_load_percent(self, mutations_num, distinct=True):
        return self.count_lower_mutation_load(distinct)/mutations_num
        
    def calc_non_top_mutation_load_percent(self, mutations_num, distinct=True):
        return self.count_non_top_mutation_load(distinct)/mutations_num
    
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
        with open(path) as f:
            while True:
                #ignoring comment lines in the begining of the file
                line = f.readline()
                if line.startswith('Hugo'):
                    break
                logger.debug("ignoring line in the begining of the file: %s", line)
            file_dict = csv.DictReader(f, dialect=csv.excel_tab, fieldnames=line.split())
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
                logger.debug("updating patient %s survival according to days_to_last_followup", patient_barcode)
                sample.update_survival(days_to_last_followup, update_date)
            else:
                logger.debug("updating patient %s survival according to days_to_death", patient_barcode)
                sample.update_survival(days_to_death, update_date)
        
                
    def write_output(self, output_path, cancer):
        logger.info("writing output file %s", output_path)
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
        self.find_high_low_mutation_load_patients()
        logger.info("writing mutation load output file %s", output_path)
        with open(output_path, "w") as f:
            csv_file = csv.DictWriter(f, fieldnames=["Hugo_code", "samples_count",
                                                     "top_mutation_load_samples_count", 
                                                     "lower_mutation_load_samples_count", 
                                                     "total_patients_ratio", 
                                                     "top_mutation_load_patients_ratio", 
                                                     "lower_mutation_load_patients_ratio", 
                                                     "non_top_mutation_load_patients_ratio", 
                                                     "top_patients_mutations_ratio", 
                                                     "lower_patients_mutations_ratio", 
                                                     "non_top_patients_mutations_ratio",
                                                     "samples_count_distinct", 
                                                     "top_mutation_load_samples_count_distinct", 
                                                     "lower_mutation_load_samples_count_distinct", 
                                                     "total_patients_ratio_distinct", 
                                                     "top_mutation_load_patients_ratio_distinct", 
                                                     "lower_mutation_load_patients_ratio_distinct", 
                                                     "non_top_mutation_load_patients_ratio_distinct", 
                                                     "top_patients_mutations_ratio_distinct", 
                                                     "lower_patients_mutations_ratio_distinct", 
                                                     "non_top_patients_mutations_ratio_distinct",
                                                     "cancer"])
            csv_file.writeheader()
            patients_num = len(self.ids_dict.keys())
            top_mutation_load_num = len([i for i in self.ids_dict.values() if i.top_mutation_load])
            lower_mutation_load_num = len([i for i in self.ids_dict.values() if i.lower_mutation_load])
            top_patients_mutations_sum = sum([i.count_mutations() for i in self.ids_dict.values() if i.top_mutation_load])
            patients_mutations_sum = sum([i.count_mutations() for i in self.ids_dict.values()])
            non_top_patients_mutations_sum = patients_mutations_sum - top_patients_mutations_sum
            for mut in self.mutations_dict.values():
                row_dict = {"Hugo_code": mut.hugo_code,
                           "samples_count_distinct" : mut.count_all_mutation_load(),
                           "top_mutation_load_samples_count_distinct" : mut.count_top_mutation_load(),
                           "lower_mutation_load_samples_count_distinct" : mut.count_lower_mutation_load(),
                           "top_mutation_load_patients_ratio_distinct" : mut.calc_top_mutation_load_patients_ratio(top_mutation_load_num),
                           "lower_mutation_load_patients_ratio_distinct" : mut.calc_lower_mutation_load_patients_ratio(lower_mutation_load_num),
                           "total_patients_ratio_distinct" : mut.calc_patients_ratio(patients_num),
                           "non_top_mutation_load_patients_ratio_distinct" : mut.calc_non_top_mutation_load_patients_ratio(patients_num - top_mutation_load_num),
                           "top_patients_mutations_ratio_distinct" : mut.calc_top_mutation_load_percent(top_patients_mutations_sum),
                           "lower_patients_mutations_ratio_distinct" : mut.calc_top_mutation_load_percent(lower_patients_mutations_sum),
                           "non_top_patients_mutations_ratio_distinct" : mut.calc_non_top_mutation_load_percent(non_top_patients_mutations_sum),
                           "samples_count" : mut.count_all_mutation_load(False),
                           "top_mutation_load_samples_count" : mut.count_top_mutation_load(False),
                           "lower_mutation_load_samples_count" : mut.count_lower_mutation_load(False),
                           "top_mutation_load_patients_ratio" : mut.calc_top_mutation_load_patients_ratio(top_mutation_load_num, False),
                           "lower_mutation_load_patients_ratio" : mut.calc_lower_mutation_load_patients_ratio(lower_mutation_load_num, False),
                           "total_patients_ratio" : mut.calc_patients_ratio(patients_num,False),
                           "non_top_mutation_load_patients_ratio" : mut.calc_non_top_mutation_load_patients_ratio(patients_num - top_mutation_load_num,False),
                           "top_patients_mutations_ratio" : mut.calc_top_mutation_load_percent(top_patients_mutations_sum,False),
                           "lower_patients_mutations_ratio" : mut.calc_lower_mutation_load_percent(lower_patients_mutations_sum,False),
                           "non_top_patients_mutations_ratio" : mut.calc_non_top_mutation_load_percent(non_top_patients_mutations_sum,False),
                           "cancer" : cancer}
                csv_file.writerow(row_dict)
                    
    def find_high_low_mutation_load_patients(self):
        patients_list = sorted(self.ids_dict.values(), key=lambda x: x.count_mutations(), reverse=True)
        stop_index = int(len(patients_list)/TOP_PERCENTIL)
        for patient in patients_list[:stop_index]:
            patient.update_top_mutation_load()
        for patient in patients_list[-stop_index:]:
            patient.update_lower_mutation_load()
                
                
    def write_survival_output(self, output_path, cancer):
        logger.info("writing survival output file %s", output_path)
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
                    
        
            
def main():
    parser = OptionParser()
    parser.add_option("-c", "--cancer", dest="cancer", help="tumor main site")
    parser.add_option("-p", "--csv_path", dest="csv_path", help="csv paths, use '' if the path contains *")
    parser.add_option("--clinical_path", dest="clinical_path", help="clinical paths, use '' if the path contains *")
    parser.add_option("--debug", default=False, action="store_false", dest="debug", help="ruu the script in debug mode")
    
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
    
    if options.clinical_path:
        clinical_paths = glob.glob(options.clinical_path)
    else:
        clinical_paths = []
    summary = MutationsSummary(options.csv_path, clinical_paths)
    summary.write_mutation_load_output("mutations_load_%s.csv" % options.cancer, options.cancer)
    
if __name__ == "__main__":
    main()
        
        
