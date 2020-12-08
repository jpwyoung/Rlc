#!/usr/bin/env python3
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import subprocess
import numpy as np

'''
Author: Peter Young (peter.young@york.ac.uk).

This script takes a set of single-gene alignments (created by find_genes.py) and 
concatenates them to produce a file with one long concatenated sequence per strain.
'''

source_folder = "../output/find_genes/genes_alignments/"
genome_file_list = "../data/genome_file_list.txt"
output_file = "../output/find_genes/concatenated_genes.fas"
alignment_file_list = "../output/find_genes/genes_alignments/alignment_file_list.txt"

#Exclude genomes that have too many missing genes
missing_genes_allowed = 20

#Make a list of all the alignment files
subprocess.call("ls " + source_folder + "*.fas > " + alignment_file_list, shell=True)

strain_list = []
with open(genome_file_list) as genome_list:
    for genome in genome_list:
        strain = genome.split("/")[-1][:-5]
        strain_list.append(strain)

#Initialise concatenation_dict with an empty sequence record for each strain
concatenation_dict = {}
missing_gene_count = {}
for strain in strain_list:
    concatenation_dict[strain] = SeqRecord(Seq(""), id=strain, description="")
    missing_gene_count[strain] = 0
    
#Define a pad_seq to separate the genes in the concatenated sequences
pad_seq = Seq("NNN")

with open(alignment_file_list) as list_file:

    for alignment_file in list_file:
        strains_to_use = strain_list[:]
        with open( alignment_file.strip("\n")) as alignment:
            #append the aligned sequence to the end of the concatenated sequence for each strain        
            for seq_record in SeqIO.parse(alignment, "fasta"):
                strain = seq_record.id
                if strain in strains_to_use:                
                    concatenation_dict[strain].seq += pad_seq + seq_record.seq
                    strains_to_use.remove(strain)
                    gene_length = len(seq_record.seq)
        #append a string of gaps for strains that do not have the gene        
        missing_gene_seq = Seq(gene_length*"-")
        for strain in strains_to_use:
            print("Strain " + strain + " missing from " + alignment_file)
            concatenation_dict[strain].seq += pad_seq + missing_gene_seq
            missing_gene_count[strain] += 1

#Remove strains that have too many missing genes
for strain in strain_list:
    if missing_gene_count[strain] > missing_genes_allowed:
        concatenation_dict.pop(strain)
            
#Write out the concatenated alignment
with open(output_file, "w") as outfile:
    SeqIO.write(concatenation_dict.values(), outfile, "fasta")
        