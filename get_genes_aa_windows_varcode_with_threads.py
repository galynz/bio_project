# -*- coding: utf-8 -*-
"""
Created on Sun Mar 05 10:37:43 2017

@author: gal
"""
import gzip, csv
import varcode
from varcode import Variant
from pyensembl import ensembl_grch38
import subprocess
import sys
import re
import threading
import logging

LOG_FILENAME = 'count_neoantigen.log'
logging.basicConfig(level=logging.DEBUG, filename=LOG_FILENAME,
                    format='[%(levelname)s] (%(threadName)-10s) %(message)s',
                    )


class ActiveMHCpanThreads(object):
    def __init__(self):
        self.active = []
        self.lock = threading.Lock()
    def makeActive(self, name):
        with self.lock:
            self.active.append(name)
            logging.debug('Running: %s', self.active)
    def makeInactive(self):
        with self.lock:
            self.active.remove(name)
            logging.debug('Running %s', self.active)
            

def read_maf(path):
    with gzip.open(path) as f:
        line = f.readline().lower()
        while line.startswith('hugo')==False:
            #ignoring comment lines in the begining of the file
            line = f.readline().lower()
        file_dict = [i for i in csv.DictReader(f, dialect=csv.excel_tab, fieldnames=line.split())]
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
            
                    
def count_neoantigen(aa_frame_list, hla_list, tmp_output, pool, s, lock, writer, output_file, tcga_id):
    counter = 0
    for aa_frame in aa_frame_list:
        with s:
            name = threading.currentThread().getName()
            pool.makeActive(name)
            with open(tmp_output, "wb") as f:
                f.write(">\n" + aa_frame)
                out = subprocess.check_output(["/groups/nshomron/galyankovitz/bio_project/installations/netMHCpan-3.0/netMHCpan",
                                       "-a", ",".join(hla_list), 
                                       "-f", tmp_output])
                if re.compile(".*Number of high binders [1-9].*").findall(out):
                counter += 1
            pool.makeInactive(name)
    logging.debug("Number of neoantigens: %d", count)
    writer.writerow({"tcga_id" : tcga_id, "neoantigen_count" : count})
    logging.debug("Wrote output")
    output_file.flush()
    logging.debug("Done!", tcga_id)


def main():
    print "Started the script"
    hla_path = sys.argv[1]
    print hla_path
    maf_path = sys.argv[2]
    print maf_path
    output_path = sys.argv[3]
    print output_path
    file_dict = read_maf(maf_path)
    with open(output_path, "wb") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=["tcga_id", "neoantigen_count"])
        writer.writeheader()
        print "Started writing output"
        with open(hla_path, "rb") as hla_file:
            csv_reader = csv.DictReader(hla_file)
            threads = []
            lock = threading.Lock()
            pool = ActiveMHCpanThreads()
            s = threading.Semaphore(10)
            for hla_info in csv_reader:
                tcga_id = hla_info["tcga_id"]
                hla_list = list(set(["HLA-" + i.replace("*", "") for i in [hla_info["A1"], hla_info["A2"], hla_info["B1"], hla_info["B2"], hla_info["C1"], hla_info["C2"]]]))
                logging.debug("Starting, tgca_id: %s, hla_list: %s", tcga_id, ",".join(hla_list))
                mutations = get_sample_mutations(file_dict, tcga_id)
                logging.debug("Found %d mutations for tcga_id %s", len(mutations), tcga_id)
                var_list = get_variants(mutations)
                logging.debug("Found %d vars for tcga_id %s", len(var_list), tcga_id)
                effects = get_top_varirants(var_list)
                logging.debug("Found %d top effects for tcga_id %s", len(effects), tcga_id)
                aa_frame_list = get_aa_frame(effects)
                logging.debug("AA frame list len: %d for tcga_id %s", len(aa_frame_list), tcga_id)
                t = threading.Thread(target=count_neoantigen, args=(aa_frame_list, hla_list, "tmp2.txt", pool, s, lock, writer, output_file, tcga_id), name=tcga_id)
                threads.append(t)
                t.start()                
                #count = count_neoantigen(aa_frame_list, hla_list, "tmp2.txt", pool, s, lock, writer, output_file)
                                
                
                
            
if __name__ == "__main__":
    main()
