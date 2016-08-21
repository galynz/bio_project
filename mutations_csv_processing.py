# -*- coding: utf-8 -*-
"""
Created on Sun Aug 21 10:33:37 2016

@author: gal
"""

import csv
import glob
import os

BEST_RES_NUM = 100
MIN_SAMPLES_COUNT = 25

class MutataionSummary(object):
    def __init__(self, csv_path):
        self.csv_path = csv_path
        with open(self.csv_path) as f:
            self.csv_dict = csv.DictReader(f)
            self.create_mutations_dict()
    
    def create_mutations_dict(self):
        self.mutation_dict = {}
        for row in self.csv_dict:
            self.mutation_dict[row['Hugo_code']] = row
            
    def get_highest_top_mutation_load_count(self):
        mut_count_row_list = [(hugo_code, row['top_mutation_load_samples_count'], row) for hugo_code, row in self.mutation_dict.items()]
        mut_count_row_list.sort(key=lambda x: int(x[1]), reverse=True)
        return mut_count_row_list[:BEST_RES_NUM]
        
    def get_highest_top_non_top_mutations_ratio(self):
        mut_ratio_row_list = []
        for hugo_code, row in self.mutation_dict.items():
            top_patients_mutations_ratio = float(row['top_patients_mutations_ratio'])
            non_top_patients_mutations_ratio = float(row['non_top_patients_mutations_ratio'])
            if non_top_patients_mutations_ratio == 0.0:
#                non_top_patients_mutations_ratio = 0.00001
                ratio = 'only top patients (%s)' % top_patients_mutations_ratio
            else:
                ratio = top_patients_mutations_ratio/non_top_patients_mutations_ratio
#            ratio = top_patients_mutations_ratio/non_top_patients_mutations_ratio
            mut_ratio_row_list += [(hugo_code, ratio, row['samples_count'], row)]
#            if int(row['samples_count']) >= MIN_SAMPLES_COUNT and ratio>1:
#                mut_ratio_row_list += [(hugo_code, ratio, row['samples_count'], row)]
        return sorted(mut_ratio_row_list, key=lambda x: x[1], reverse=True)
        
        
class CancerMutations(object):
    def __init__(self, csv_dir):
        self.csv_dir = csv_dir
        self.csv_files_dict = {}
        search_path = os.path.join(csv_dir, '*_20160821_without_silent_mut.mutation_load.csv')
        for path in glob.glob(search_path):
            filename = os.path.basename(path)
            end_pos = filename.find("_2")
            cancer = filename[:end_pos]
            self.csv_files_dict[cancer] = path
        self.mutations_summary_dict = {cancer: MutataionSummary(path) for cancer, path in self.csv_files_dict.items()}
        
    def create_top_non_top_ratio(self, output_path):
        mutation_cancer_ratio_dict = {}
        for cancer, mutation_summary in self.mutations_summary_dict.items():
            for mut_info in mutation_summary.get_highest_top_non_top_mutations_ratio():
                mutation = mut_info[0]
                mutation_cancer_ratio_dict.setdefault(mutation, {'Hugo_symbol' : mutation})
                ratio = mut_info[1]
                mutation_cancer_ratio_dict[mutation][cancer] = ratio
        #return mutation_cancer_ratio_dict
        headers = ['Hugo_symbol'] + self.mutations_summary_dict.keys()
        with open(output_path, 'w') as f:
            writer = csv.DictWriter(f, headers)
            writer.writeheader()
            for cancers_dict in mutation_cancer_ratio_dict.values():
                writer.writerow(cancers_dict)