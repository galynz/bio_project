# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 13:24:33 2017

@author: gal
"""

import glob, urllib2, re, csv

#Getting ids
l = [i.replace("\n","") for i in open("/groups/nshomron/galyankovitz/storage_1_root/WXS/Ovarian_Serous_Cystadenocarcinoma/ids_list.txt", "rb")]

with open("output.csv", "wb") as csv_output:
    writer = csv.DictWriter(csv_output, fieldnames=["tcga_id", "A1", "A2", "B1", "B2", "C1", "C2"])
    writer.writeheader()
    for i in l:
        if glob.glob("/groups/nshomron/galyankovitz/storage_1_root/WXS/Ovarian_Serous_Cystadenocarcinoma/" + i + "/2017_03_12*/*result.tsv"):
            res_path = glob.glob("/groups/nshomron/galyankovitz/storage_1_root/WXS/Ovarian_Serous_Cystadenocarcinoma/" + i + "/2017_03_12*/*result.tsv")[-1]
            res_con = csv.DictReader(open(res_path, "rb"), delimiter="\t").next()
            response = urllib2.urlopen("https://gdc-api.nci.nih.gov/files/"+ i +"?fields=cases.submitter_id&pretty=true").read()
            tcga_id = re.compile('.*"(TCGA-.*)".*').findall(response)[0]
            writer.writerow({"tcga_id" : tcga_id, "A1" : res_con["A1"], "A2" : res_con["A2"], "B1" : res_con["B1"], "B2" : res_con["B2"], "C1" : res_con["C1"], "C2" : res_con["C2"]})