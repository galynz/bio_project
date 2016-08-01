# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 16:08:53 2016

@author: gal
"""

#Imports
from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, String
from sqlalchemy.orm import sessionmaker
import csv

ENGINE = create_engine('sqlite:///:memory:', echo=True)
Base = declarative_base()


class Maf(Base):
    __tablename__ = 'maf_files'
    
    id = Column(Integer, primary_key=True)
    Hugo_Symbol = Column(String)
    Entrez_Gene_Id = Column(String)
    Center = Column(String)
    NCBI_Build = Column(String)
    Chromosome = Column(String)
    Start_Position = Column(String)
    End_Position = Column(String)
    Strand = Column(String)
    Variant_Classification = Column(String)
    Variant_Type = Column(String)
    Reference_Allele = Column(String)
    Tumor_Seq_Allele1 = Column(String)
    Tumor_Seq_Allele2 = Column(String)
    dbSNP_RS = Column(String)
    dbSNP_Val_Status = Column(String)
    Tumor_Sample_Barcode = Column(String)
    Matched_Norm_Sample_Barcode = Column(String)
    Match_Norm_Seq_Allele1 = Column(String)
    Match_Norm_Seq_Allele2 = Column(String)
    Tumor_Validation_Allele1 = Column(String)
    Tumor_Validation_Allele2 = Column(String)
    Match_Norm_Validation_Allele1 = Column(String)
    Match_Norm_Validation_Allele2 = Column(String)
    Verification_Status = Column(String)
    Validation_Status = Column(String)
    Mutation_Status = Column(String)
    Sequencing_Phase = Column(String)
    Sequence_Source = Column(String)
    Validation_Method = Column(String)
    Score = Column(String)
    BAM_File = Column(String)
    Sequencer = Column(String)
    Tumor_Sample_UUID = Column(String)
    Matched_Norm_Sample_UUID = Column(String)
    chromosome_name_WU = Column(String)
    start_WU = Column(String)
    stop_WU = Column(String)
    reference_WU = Column(String)
    variant_WU = Column(String)
    type_WU = Column(String)
    gene_name_WU = Column(String)
    transcript_name_WU = Column(String)
    transcript_species_WU = Column(String)
    transcript_source_WU = Column(String)
    transcript_version_WU = Column(String)
    strand_WU = Column(String)
    transcript_status_WU = Column(String)
    trv_type_WU = Column(String)
    c_position_WU = Column(String)
    amino_acid_change_WU = Column(String)
    ucsc_cons_WU = Column(String)
    domain_WU = Column(String)
    all_domains_WU = Column(String)
    deletion_substructures_WU = Column(String)
    transcript_error = Column(String)
    
Base.metadata.create_all(ENGINE)

def load_csv(csv_path):
    #Create the session
    session = sessionmaker()
    session.configure(bind=ENGINE)
    s = session()
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        try:
            for row in reader:
                new_line = Maf(**row)
                s.add(new_line)
            s.commit()
        except:
            s.rollback()