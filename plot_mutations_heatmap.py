# -*- coding: utf-8 -*-
"""
Created on Sun Nov 27 13:00:41 2016

@author: gal
"""

import gzip, re, csv, datetime, sys, glob, os
import xml.etree.ElementTree as ET
import logging
import logging.handlers
from optparse import OptionParser
import matplotlib
matplotlib.use('Agg')
from plotly.offline import download_plotlyjs, init_notebook_mode,  plot
import plotly.graph_objs as go
import numpy as np
import pandas as pd
import lifelines as ll
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import cophenet
from scipy.spatial.distance import pdist
import pylab

# Plotting helpers
#from IPython.display import HTML
import matplotlib.pyplot as plt
import plotly.plotly as py
import plotly.tools as tls   
from plotly.graph_objs import *
from lifelines.statistics import logrank_test

logger = logging.getLogger("mutations_heatmap")

FRACTION = 0.15

BIOTYPE_PRIORITY = {'protein_coding' : 1, # Contains an open reading frame (ORF)
                    'LRG_gene' : 2, # Gene in a "Locus Reference Genomic" region known to have disease-related sequence variations
                    'IG_C_gene' : 2, # Immunoglobulin (Ig) variable chain genes imported or annotated according to the IMGT
                    'IG_D_gene' : 2, # Immunoglobulin (Ig) variable chain genes imported or annotated according to the IMGT
                    'IG_J_gene' : 2, # Immunoglobulin (Ig) variable chain genes imported or annotated according to the IMGT
                    'IG_LV_gene' : 2, # Immunoglobulin (Ig) variable chain genes imported or annotated according to the IMGT
                    'IG_V_gene' : 2, # Immunoglobulin (Ig) variable chain genes imported or annotated according to the IMGT
                    'TR_C_gene' : 2, # T-cell receptor (TcR) genes imported or annotated according to the IMGT
                    'TR_D_gene' : 2, # T-cell receptor (TcR) genes imported or annotated according to the IMGT
                    'TR_J_gene' : 2, # T-cell receptor (TcR) genes imported or annotated according to the IMGT
                    'TR_V_gene' : 2, # T-cell receptor (TcR) genes imported or annotated according to the IMGT
                    'miRNA' : 3, # Non-coding RNA predicted using sequences from RFAM and miRBase
                    'snRNA' : 3, # Non-coding RNA predicted using sequences from RFAM and miRBase
                    'snoRNA' : 3, # Non-coding RNA predicted using sequences from RFAM and miRBase
                    'ribozyme' : 3, # Non-coding RNA predicted using sequences from RFAM and miRBase
                    'tRNA' : 3, #Added by Y. Boursin
                    'sRNA' : 3, # Non-coding RNA predicted using sequences from RFAM and miRBase
                    'scaRNA' : 3, # Non-coding RNA predicted using sequences from RFAM and miRBase
                    'rRNA' : 3, # Non-coding RNA predicted using sequences from RFAM and miRBase
                    'lincRNA' : 3, # Long, intervening noncoding (linc) RNAs, that can be found in evolutionarily conserved, intergenic regions
                    'bidirectional_promoter_lncrna' : 3, # A non-coding locus that originates from within the promoter region of a protein-coding gene, with transcription proceeding in the opposite direction on the other strand
                    'bidirectional_promoter_lncRNA' : 3, # A non-coding locus that originates from within the promoter region of a protein-coding gene, with transcription proceeding in the opposite direction on the other strand
                    'known_ncrna' : 4,
                    'vaultRNA' : 4, # Short non coding RNA genes that form part of the vault ribonucleoprotein complex
                    'macro_lncRNA' : 4, # unspliced lncRNAs that are several kb in size
                    'Mt_tRNA' : 4, # Non-coding RNA predicted using sequences from RFAM and miRBase
                    'Mt_rRNA' : 4, # Non-coding RNA predicted using sequences from RFAM and miRBase
                    'antisense' : 5, # Has transcripts that overlap the genomic span (i.e. exon or introns) of a protein-coding locus on the opposite strand
                    'antisense_RNA' : 5, # Alias for antisense (Y. Boursin)
                    'sense_intronic' : 5, # Long non-coding transcript in introns of a coding gene that does not overlap any exons
                    'sense_overlapping' : 5, # Long non-coding transcript that contains a coding gene in its intron on the same strand
                    '3prime_overlapping_ncrna' : 5, # Transcripts where ditag and/or published experimental data strongly supports the existence of short non-coding transcripts transcribed from the 3'UTR
                    '3prime_overlapping_ncRNA' : 5, # Transcripts where ditag and/or published experimental data strongly supports the existence of short non-coding transcripts transcribed from the 3'UTR
                    'misc_RNA' : 5, # Non-coding RNA predicted using sequences from RFAM and miRBase
                    'non_coding' : 5, # Transcript which is known from the literature to not be protein coding
                    'regulatory_region' : 6, # A region of sequence that is involved in the control of a biological process
                    'disrupted_domain' : 6, # Otherwise viable coding region omitted from this alternatively spliced transcript because the splice variation affects a region coding for a protein domain
                    'processed_transcript' : 6, # Doesn't contain an ORF
                    'TEC' : 6, # To be Experimentally Confirmed. This is used for non-spliced EST clusters that have polyA features. This category has been specifically created for the ENCODE project to highlight regions that could indicate the presence of protein coding genes that require experimental validation, either by 5' RACE or RT-PCR to extend the transcripts, or by confirming expression of the putatively-encoded peptide with specific antibodies
                    'TF_binding_site' : 7, # A region of a nucleotide molecule that binds a Transcription Factor or Transcription Factor complex
                    'CTCF_binding_site' :7, # A transcription factor binding site with consensus sequence CCGCGNGGNGGCAG, bound by CCCTF-binding factor
                    'promoter_flanking_region' : 7, # A region immediately adjacent to a promoter which may or may not contain transcription factor binding sites
                    'enhancer' : 7, # A cis-acting sequence that increases the utilization of (some) eukaryotic promoters, and can function in either orientation and in any location (upstream or downstream) relative to the promoter
                    'promoter' : 7, # A regulatory_region composed of the TSS(s) and binding sites for TF_complexes of the basal transcription machinery
                    'open_chromatin_region' : 7, # A DNA sequence that in the normal state of the chromosome corresponds to an unfolded, un-complexed stretch of double-stranded DNA
                    'retained_intron' : 7, # Alternatively spliced transcript believed to contain intronic sequence relative to other, coding, variants
                    'nonsense_mediated_decay' : 7, # If the coding sequence (following the appropriate reference) of a transcript finishes >50bp from a downstream splice site then it is tagged as NMD. If the variant does not cover the full reference coding sequence then it is annotated as NMD if NMD is unavoidable i.e. no matter what the exon structure of the missing portion is the transcript will be subject to NMD
                    'non_stop_decay' : 7, # Transcripts that have polyA features (including signal) without a prior stop codon in the CDS, i.e. a non-genomic polyA tail attached directly to the CDS without 3' UTR. These transcripts are subject to degradation
                    'ambiguous_orf' : 7, # Transcript believed to be protein coding, but with more than one possible open reading frame
                    'pseudogene' : 8, # Have homology to proteins but generally suffer from a disrupted coding sequence and an active homologous gene can be found at another locus. Sometimes these entries have an intact coding sequence or an open but truncated ORF, in which case there is other evidence used (for example genomic polyA stretches at the 3' end) to classify them as a pseudogene. Can be further classified as one of the following
                    'processed_pseudogene' : 8, # Pseudogene that lack introns and is thought to arise from reverse transcription of mRNA followed by reinsertion of DNA into the genome
                    'polymorphic_pseudogene' : 8, # Pseudogene owing to a SNP/DIP but in other individuals/haplotypes/strains the gene is translated
                    'retrotransposed' : 8, # Pseudogene owing to a reverse transcribed and re-inserted sequence
                    'translated_processed_pseudogene' : 8, # Pseudogenes that have mass spec data suggesting that they are also translated
                    'translated_unprocessed_pseudogene' : 8, # Pseudogenes that have mass spec data suggesting that they are also translated
                    'transcribed_processed_pseudogene' : 8, # Pseudogene where protein homology or genomic structure indicates a pseudogene, but the presence of locus-specific transcripts indicates expression
                    'transcribed_unprocessed_pseudogene' : 8, # Pseudogene where protein homology or genomic structure indicates a pseudogene, but the presence of locus-specific transcripts indicates expression
                    'transcribed_unitary_pseudogene' : 8, #Pseudogene where protein homology or genomic structure indicates a pseudogene, but the presence of locus-specific transcripts indicates expression
                    'unitary_pseudogene' : 8, # A species specific unprocessed pseudogene without a parent gene, as it has an active orthologue in another species
                    'unprocessed_pseudogene' : 8, # Pseudogene that can contain introns since produced by gene duplication
                    'Mt_tRNA_pseudogene' : 8, # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
                    'tRNA_pseudogene' : 8, # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
                    'snoRNA_pseudogene' : 8, # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
                    'snRNA_pseudogene' : 8, # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
                    'scRNA_pseudogene' : 8, # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
                    'rRNA_pseudogene' : 8, # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
                    'misc_RNA_pseudogene' : 8, # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
                    'miRNA_pseudogene' : 8, # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
                    'IG_C_pseudogene' : 8, # Inactivated immunoglobulin gene
                    'IG_D_pseudogene' : 8, # Inactivated immunoglobulin gene
                    'IG_J_pseudogene' : 8, # Inactivated immunoglobulin gene
                    'IG_V_pseudogene' : 8, # Inactivated immunoglobulin gene
                    'TR_J_pseudogene' : 8, # Inactivated immunoglobulin gene
                    'TR_V_pseudogene' : 8, # Inactivated immunoglobulin gene
                    'artifact' : 9, # Used to tag mistakes in the public databases (Ensembl/SwissProt/Trembl)
    }

class Sample(object):
    def __init__(self, patient_id):
        self.patient_id = patient_id
        self.germline_mutations = dict()
        self.somatic_mutations = set()
        
    def add_germline_mutation(self, gene, mutation_priority):
        self.germline_mutations[gene] = {'severity' : min(mutation_priority, self.germline_mutations.get(gene, {}).get('severity', 10)), 
                                        'count' : self.germline_mutations.get(gene, {}).get('count', 0) + 1}
        
    def add_somatic_mutation(self, gene):
        self.somatic_mutations.add(gene)
        
    def count_somatic_mutations(self):
        return len(self.somatic_mutations)
    
    def count_germline_mutations(self, gene):
        return self.germline_mutations.get(gene, {}).get('count', 0)
        
        if self.germline_mutations.get(gene, {}).get('count', 0) == 0:
            return 0
        return 1
        
        
def parse_vcf(vcf_path, samples_dict):#, vcf_list_file):
    try:
        with gzip.open(vcf_path) as f:
            fields = []
            logger.info("going over %s", vcf_path)
            n = 0
            for line in f:
                n += 1
                if line.startswith("chr"):
                    if line.find("germline_risk") > -1:
                        info = [i.split("|") for i in line[line.find("CSQ=")+len("CSQ="):line.find(" ", line.find("CSQ=")+len("CSQ="))].split(",")]
                        for i in info:
                            if len(i) == len(fields):
                                hugo_code = i[fields.index("SYMBOL")]
                                biotype = i[fields.index("BIOTYPE")]
                                if BIOTYPE_PRIORITY.has_key(biotype):
                                    priority = BIOTYPE_PRIORITY[biotype]
                                    samples_dict[sample].add_germline_mutation(hugo_code, priority)
                                    logger.debug("added gene %s with priority %d to sample %s", hugo_code, priority, sample)
                                else:
                                    logger.debug("didn't add gene %s with biotype %s to sample %s because the biotype is not recognized", hugo_code, biotype, sample)
                elif line.startswith("##INDIVIDUAL"):
                    sample = re.compile(r'INDIVIDUAL=<NAME=(TCGA-[A-Z0-9-]*),').findall(line)[0]
                    if not samples_dict.has_key(sample):
                        samples_dict.setdefault(sample, Sample(sample))
                        logger.warn("%s - the id isn't in the ids list", vcf_path)
                elif line.startswith("##INFO"):
                    fields = re.compile("([\w_]{1,})[\|>\"]").findall(line)
                elif line.startswith("#CHROM"):
                    if not fields:
                        logger.warn("failed to parse %s", vcf_path)
                        return
                elif line.startswith("##gdcWorkflow") and line.find("mutect2") == -1:
                    return
#            vcf_list_file.write(vcf_path + "\n")
            logger.debug("went over %d lines in file %s", n, vcf_path)
    except IOError:
        logger.exception("can't read %s", vcf_path)
            
def add_csv_data(path, samples_dict):
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
                patient_barcode = "-".join(tumor_barcode.split('-')[:3])
                gene = row["hugo_symbol"]
                mutation_type = row["variant_classification"]
                if mutation_type != "Silent":
                    sample = samples_dict.setdefault(patient_barcode, Sample(patient_barcode))
                    sample.add_somatic_mutation(gene)
                    
def create_df(samples_dict, z_param):
    # Creating a set of all the genes    
    all_genes_dict = {}
    for sample in samples_dict.values():
        for gene in sample.germline_mutations:
            all_genes_dict[gene] = all_genes_dict.get(gene, 0) + 1
    all_genes = [i[0] for i in sorted(all_genes_dict.items(), key=lambda x: x[1], reverse=True)] #Sorted, incase we'll want it  sorted in the future (taking the 500 most communly mutated genes didn't work)
    logger.info("germline genes num: %d", len(all_genes))

    # Creating a df            
    l = []
    for sample in samples_dict.values():
        if z_param == 'severity':
            l_sample = [sample.count_somatic_mutations(), sample.patient_id] + [-1*sample.germline_mutations.get(i, 10) for i in all_genes]
        elif z_param == 'germline_count':
            l_sample = [sample.count_somatic_mutations(), sample.patient_id] + [sample.count_germline_mutations(i) for i in all_genes]
        elif z_param == 'germline_binary':
            func = lambda x : int(x>0)
            l_sample = [sample.count_somatic_mutations(), sample.patient_id] + [func(sample.count_germline_mutations(i)) for i in all_genes]
        l.append(l_sample)
    tmp_df = pd.DataFrame(data=l, columns=["somatic_mutations_count", "patient_id"]+all_genes)
    df = tmp_df.sort_values("somatic_mutations_count")

    return df, all_genes

def plot_heatmap_var(samples_dict, output_path, cancer, df, all_genes):
    #Creating a list of top var genes
    genes_var_list = []
    for gene in all_genes:
        gene_var = df[gene].var()
        genes_var_list.append((gene, gene_var))
    genes_var_list.sort(key=lambda x: x[1], reverse=True)
    top_genes = [i[0] for i in genes_var_list][:500]
    plot_heatmap(samples_dict, output_path + ".top_var_genes", cancer, df, top_genes)
    
    
def plot_heatmap_top_low_unique(samples_dict, output_path, cancer, df, all_genes):
    #top mutation load patients genes
    top_mutation_load_patients_ix = (df["somatic_mutations_count"] >= df["somatic_mutations_count"].quantile(1-FRACTION))
    top_mutation_load_df = df[top_mutation_load_patients_ix]
    low_mutation_load_patients_ix = (df["somatic_mutations_count"] <= df["somatic_mutations_count"].quantile(FRACTION))
    low_mutation_load_df = df[low_mutation_load_patients_ix]
    top_genes = []
    low_genes = []
    for gene in all_genes:
        if top_mutation_load_df[gene].sum() > 0 and low_mutation_load_df[gene].sum() == 0:
            top_genes.append(gene)
        elif top_mutation_load_df[gene].sum() == 0 and low_mutation_load_df[gene].sum() > 0:
            low_genes.append(gene)
    logger.info("top_genes: %d, low_genes: %d, all_genes: %d", len(top_genes), len(low_genes), len(all_genes))
    plot_heatmap(samples_dict, output_path + ".top_low_genes", cancer, df, top_genes+low_genes)
        
def plot_heatmap(samples_dict, output_path, cancer, df, genes):
    # Plotting the heatmap
    heatmap_trace = go.Heatmap(z=[df[i] for i in genes], y=genes, x=df.patient_id, showscale=False, colorscale=[[0, "rgb(111, 168, 220)"], [1, "rgb(5, 10, 172)"]])
    mutation_load_trace = go.Bar(x=df.patient_id, y=df.somatic_mutations_count/30.0)
    fig = tls.make_subplots(rows=29, cols=1, specs=[[{'rowspan':5, 'colspan' : 1}]] + [[None]] * 4 + [[{'rowspan' : 24, 'colspan' : 1}]] + [[None]] * 23)
    fig.append_trace(mutation_load_trace, 1, 1)
    fig.append_trace(heatmap_trace, 6, 1)
    fig['layout']['xaxis1'].update(showticklabels = False)
    fig['layout']['xaxis1'].update(zeroline = False, showgrid=False)
    fig['layout']['yaxis1'].update(zeroline = False, showgrid = False, tickfont=dict(family='Arial', size=4))
    fig['layout']['xaxis2'].update(showticklabels = False)
    fig['layout']['xaxis2'].update(zeroline = False, showgrid=False)
    fig['layout']['yaxis2'].update(zeroline = False, showgrid = False, tickfont=dict(family='Arial', size=4))
    plot(fig, auto_open=False, filename="%s_%s_heatmap.html" % (output_path, cancer))
    
    plot_clustered_heatmap(df, genes, cancer, output_path)
    
def plot_clustered_heatmap(df, genes_list, cancer, output_path):
    # Build nxm matrix (n samples, m genes)
    X = df[genes_list].as_matrix().transpose()
    
    Z = linkage(X, method='complete', metric='hamming')
    c, coph_dists = cophenet(Z, pdist(X))
    print "Cophenetic Correlation Coefficient:", c
    
    #layout = go.Layout(yaxis=dict(title='%s germline mutations (ordered by samples somatic mutation load)'% cancer, zeroline=False))    
#    fig = pylab.figure(figsize=(8,8))
#    ax1 = fig.add_axes([0.09,0.1,0.2,0.6])
#    ax1.set_xticks([])
#    ax1.set_yticks([])
#    axmatrix = fig.add_axes([0.3,0.1,0.6,0.6])
    den = dendrogram(Z, orientation='left')
    idx = den['leaves']
    X = X[idx,:]
    print "X shape:", X.shape
    genes_ordered = [genes_list[i] for i in idx]
    logger.info("ordered genes: %s", str(genes_ordered))
    
#    im = axmatrix.matshow(X, aspect='auto', origin='lower', cmap=pylab.cm.YlGnBu)
#    axmatrix.set_xticks([])
#    axmatrix.set_yticks([])
#    # Plot colorbar.
#    axcolor = fig.add_axes([0.91,0.1,0.02,0.6])
#    pylab.colorbar(im, cax=axcolor)
#    fig.savefig(output_path)
    
    # Plotting the heatmap (without the hirarchy)
    heatmap_trace = go.Heatmap(z=X.tolist(), x=df.patient_id, y=genes_ordered, showscale=False, colorscale=[[0, "rgb(111, 168, 220)"], [1, "rgb(5, 10, 172)"]])
    mutation_load_trace = go.Bar(x=df.patient_id, y=df.somatic_mutations_count/30.0)
    fig = tls.make_subplots(rows=29, cols=1, specs=[[{'rowspan':5, 'colspan' : 1}]] + [[None]] * 4 + [[{'rowspan' : 24, 'colspan' : 1}]] + [[None]] * 23)
    fig.append_trace(mutation_load_trace, 1, 1)
    fig.append_trace(heatmap_trace, 6, 1)
    fig['layout']['xaxis1'].update(showticklabels = False)
    fig['layout']['xaxis1'].update(zeroline = False, showgrid=False)
    fig['layout']['yaxis1'].update(zeroline = False, showgrid = False, tickfont=dict(family='Arial', size=4))
    fig['layout']['xaxis2'].update(showticklabels = False)
    fig['layout']['xaxis2'].update(zeroline = False, showgrid=False)
    fig['layout']['yaxis2'].update(zeroline = False, showgrid = False, tickfont=dict(family='Arial', size=4))
    plot(fig, auto_open=False, filename="%s_%s_heatmap_clustered.html" % (output_path, cancer))
    
    
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
    parser.add_option("--vcf_path", dest="vcf_path", help="clinical paths, use '' if the path contains *")
    parser.add_option("--debug", default=False, action="store_true", dest="debug", help="run the script in debug mode")
    parser.add_option("-o", "--output", dest="output_path",default="output", help="where to save all the html files")
    parser.add_option("--vcf_list", dest="vcf_list", help="a txt file with a list of vcf paths to use")
    
    (options, args) = parser.parse_args()
    
    # Logger
    formatter = logging.Formatter('%(process)d %(asctime)s %(levelname)s %(message)s')
    handler = logging.handlers.RotatingFileHandler(
              "plot_mutations_heatmap.log", maxBytes=10*1024*1024, backupCount=5, mode="a")

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
    
    if options.vcf_list:
        vcf_paths = []        
        with open(options.vcf_list, "rb") as f:
            for line in f:
                vcf_paths.append(line.replace("\n", ""))
        
    elif options.vcf_path:
        vcf_paths = glob.glob(options.vcf_path)
        logger.info("going over %d vcf files", len(vcf_paths))
    else:
        vcf_paths = []

    csv_paths = glob.glob(options.csv_path)
    samples_dict = {}
    logger.info("parsing csv files")
    for csv_path in csv_paths:
        add_csv_data(csv_path, samples_dict)
        
    logger.info("parsing vcf files")
#    vcf_list_file = open("vcf_files.txt", "wb")
    for vcf_path in vcf_paths:
#        parse_vcf(vcf_path, samples_dict, vcf_list_file)
        parse_vcf(vcf_path, samples_dict)

    df, all_genes = create_df(samples_dict, 'germline_binary')
    print len(all_genes)
    plot_heatmap_top_low_unique(samples_dict, options.output_path, options.cancer, df, all_genes)
#    plot_heatmap_var(samples_dict, options.output_path, options.cancer, df, all_genes)
    
if __name__ == "__main__":
    main()
        
                
