# -*- coding: utf-8 -*-
"""
Created on Sun Mar 05 10:37:43 2017

@author: gal
"""
import gzip, csv
from varcode import Variant
from pyensembl import ensembl_grch38
import subprocess
import sys
import re


def read_maf(path):
    with gzip.open(path) as f:
        line = f.readline().lower()
        while line.startswith('hugo')==False:
            #ignoring comment lines in the begining of the file
            line = f.readline().lower()
        file_dict = csv.DictReader(f, dialect=csv.excel_tab, fieldnames=line.split())
    return file_dict
    
def get_sample_mutations(file_dict, tcga_id):
    return [line for line in file_dict if line["tumor_sample_barcode"].startswith(tcga_id)]
    
def get_variants(mutations):
    var_list = []
    for mut in mutations:
        contig = mut['chromosome'].replace("chr", "")
        if contig.isdigit():
            contig = int(contig)
        start = int(mut['start_position'])
        ref=mut['reference_allele']
        alt1=mut['tumor_seq_allele1']
        alt2=mut['tumor_seq_allele2']
        if alt1 != ref:
            var_list.append(Variant(contig=contig, start=start, ref=ref, alt=alt1, ensembl=ensembl_grch38))
        if alt2 != ref:
            var_list.append(Variant(contig=contig, start=start, ref=ref, alt=alt2, ensembl=ensembl_grch38))
    return var_list
    
def get_top_varirants(var_list):
    return [var.effects().top_priority_effect() for var in var_list]
    
def get_aa_frame(effect_list):
    aa_frame_list = []
    for eff in effect_list:
        if not isinstance(eff, varcode.effects.effect_classes.Silent):
            variantLocation = eff.aa_mutation_start_offset
            if variantLocation:
#                for j in xrange(9,13):
#                    for i in xrange(-j, 1):
#                        aa_frame_list.append(eff.mutant_protein_sequence[variantLocation+i:variantLocation+j+i])
                aa_frame_list.append(eff.mutant_protein_sequence[variantLocation-12:variantLocation+12])
    return aa_frame_list
    
def write_aa_frames(aa_frame_list, path):
    with open(path, "wb") as f:
        for aa_frame in aa_frame_list:
            f.write(">\n" + aa_frame + "\n")
            
                    
def count_neoantigen(aa_frame_list, hla_list, tmp_output):
    counter = 0
    for aa_frame in aa_frame_list:
        with open(tmp_output, "wb") as f:
            f.write(">\n" + aa_frame)
        out = subprocess.check_output(["/groups/nshomron/galyankovitz/bio_project/installations/netMHCpan-3.0/netMHCpan",
                                       "-a", ",".join(hla_list), 
                                       "-f", tmp_output])
        if re.compile(".*Number of high binders [1-9].*").findall(out):
            counter += 1
    return counter
    
def main():
    hla_path = sys.argv[1]
    maf_path = sys.argv[2]
    output_path = sys.argv[3]
    file_dict = read_maf(maf_path)
    with open(output_path, "wb") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=["tcga_id", "neoantigen_count"])
        writer.writeheader()
        with open(hla_path, "rb") as hla_file:
            csv_reader = csv.DictReader(hla_file)
            for info in csv_reader:
                tcga_id = info["tcga_id"]
                hla_list = ["HLA-" + i.replace("*", "") for i in [info["A1"], info["A2"], info["B1"], info["B2"], info["C1"], info["C2"]]]
                mutations = get_sample_mutations(file_dict, tcga_id)
                var_list = get_variants(mutations)
                effects = get_top_varirants(var_list)
                aa_frame_list = get_aa_frame(effects)
                count = count_neoantigen(aa_frame_list, hla_list, "tmp.txt")
                writer.writerow({"tcga_id" : tcga_id, "neoantigen_count" : count})
            
if __name__ == "__main__":
    main()