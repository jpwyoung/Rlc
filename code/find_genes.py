#!/usr/bin/env python3
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import subprocess
import numpy as np
import datetime


'''
Author: Peter Young (peter.young@york.ac.uk). November 2018.

This script searches a set of genomes (all .fna files in data_folder, or those listed in 
genome_file_list.txt) for matches to a set of protein or DNA sequences 
(target_standards_file) and extracts a set of full-length genes corresponding
to each sequence.  It then aligns each set and constructs phylogenetic trees from them.

The type of query sequence is specified by setting q_type = "protein" or "dna"..

'''

#Define the input files and files for intermediate results, as well as folders for the 
#sets of (unaligned) genes, alignments, and trees.

q_type = "protein"
#q_type = "dna"

data_folder = "../data"
output_folder = "../output/find_genes"
reference_folder = "../reference"
target_standards_file =  reference_folder + "/NodC_stds.fas"

blast_file = output_folder + "/genes_blast.tab"
log_file = output_folder + "/find_genes_log.txt"

os.mkdir(output_folder)
with open(log_file, "w") as log_name:
    log_name.write("Starting at " + str(datetime.datetime.now()) + "\n")    

os.mkdir(output_folder + "/genes_sequences")
os.mkdir(output_folder + "/genes_alignments")
os.mkdir(output_folder + "/genes_trees")
os.mkdir(output_folder + "/top_genes_hits")

type_dict = {"protein":(3,"tblastn"),"dna":(1,"blastn")}
try: 
    (bp_unit,blast_type) = type_dict[q_type]
except:
    print("q_type must be protein or dna")

target_dict = {} #will be used to collect all the gene sequences

target_standards_dict = SeqIO.to_dict(SeqIO.parse(target_standards_file, "fasta"))
for target_gene in target_standards_dict.keys():
    target_dict[target_gene] = []

#Make a list of the genome .fas files if a list does not exist already
if not os.path.exists(data_folder + "/genome_file_list.txt"):
    subprocess.call("ls " + data_folder + "/*.fna > " + data_folder +"/genome_file_list.txt", shell=True)

with open(data_folder + "/genome_file_list.txt") as genome_list:

    for genome_file in genome_list:
    
        genome_file = genome_file.rstrip("\n")
        strain_name = genome_file.split("/")[-1][:-4]
        
        #make a blast database for the genome
        subprocess.call("makeblastdb -in " + genome_file + " -input_type fasta -dbtype nucl", shell=True)
        
        #Run a blastn or tblastn search of all scaffolds with the standards file as query

        subprocess.call(blast_type + " -db " + genome_file + " -query " + target_standards_file +" -out " + 
        blast_file +" -outfmt 6 -max_hsps 1 -evalue 1E-10 -num_threads 6", shell=True)

        #Read the blast hits into a dict with queries as keys.
        #Keep only the best hit for each query in each strain.

        hit_dict = {}

        with open(blast_file) as blast_output:
            for hit in blast_output:
                hit_bits = hit.rstrip("\n").split("\t")
                query = hit_bits[0]
                if query not in hit_dict:
                    hit_dict[query] = hit_bits
                else:                
                    if float(hit_bits[11]) > float(hit_dict[query][11]):
                        #new hit is better, so replace the existing
                        hit_dict[query] = hit_bits 

        #Transfer these best hits to a dict indexed by scaffold, so that the actual gene sequences
        #can be extracted from the scaffolds file.

        top_hit_dict = {} 

        for query in hit_dict:
            hit = hit_dict[query]
            scaffold = hit[1]
            if scaffold not in top_hit_dict:
                top_hit_dict[scaffold] = [hit]
            else:
                top_hit_dict[scaffold].append(hit)

        #Write the top hits to a tab-delimited file arranged like the Blast output, sorted by scaffold and query
        
        top_hits_file = output_folder + "/top_genes_hits/" + strain_name + "_hits.tab"  
        hit_count = 0

        with open(top_hits_file, "w") as outfile:
            sorted_top_hits = sorted(top_hit_dict)
            for scaffold in sorted_top_hits:
                sorted_by_query = sorted(top_hit_dict[scaffold], key = lambda hit: hit[0])
                for hit in sorted_by_query:
                    for hit_bit in hit:
                        outfile.write(str(hit_bit) + "\t")
                    outfile.write("\n")
                    hit_count +=1

        with open(log_file, "a") as log_name:
            log_name.write(str(datetime.datetime.now()))
            log_name.write(" Found " + str(hit_count) + " genes in " + strain_name + "\n")

        if hit_count >= 1: #only include genomes with at least 1 gene
        
            #for each top hit, extract the gene sequence from the scaffold and add it to a list for that 
            #query
                
            all_scaffolds_dict = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))

            for scaffold in top_hit_dict:
                for hit in top_hit_dict[scaffold]:
                    target_class = hit[0]
                    q_start = int(hit[6])
                    q_end = int(hit[7])
                    s_start = int(hit[8])
                    s_end = int(hit[9])
                    g_dir = np.sign(s_end - s_start)
                    g_len = bp_unit*len(target_standards_dict[target_class])
                    scaf_len = len(all_scaffolds_dict[scaffold])
    
                    g_start = s_start - bp_unit*(q_start-1)*g_dir 
                    g_end = s_end + (g_len - bp_unit*q_end)*g_dir
    
                    if g_dir == 1:
                        g_left = g_start
                        g_right = g_end
                        dir = ""
                    elif g_dir == -1:
                        g_left = g_end
                        g_right = g_start 
                        dir = "c"
            
                    l_flag = r_flag = ""    
                    if g_left < 1: 
                        g_left = 1
                        l_flag = "*"
                    if g_right > scaf_len: 
                        g_right = scaf_len
                        r_flag = "*"
                    g_seq = all_scaffolds_dict[scaffold][g_left-1:g_right]
                    if g_dir == -1:
                        g_seq = g_seq.reverse_complement()
                    g_seq.id = strain_name
                    g_seq.description = scaffold + "_" + dir + l_flag + str(g_left) + "-" + str(g_right) + r_flag
    
                    target_dict[target_class].append(g_seq)
                       

for target_class in target_dict:
    target_sequence_file = output_folder + "/genes_sequences/" + target_class + "_raw.fas"
    target_alignment_file = output_folder + "/genes_alignments/" + target_class + ".fas"
    target_tree_file = output_folder + "/genes_trees/" + target_class + ".tre"

    #Save all the gene sequences in a separate file for each target gene    
    with open(target_sequence_file, "w") as outfile:
        SeqIO.write(target_dict[target_class], outfile, "fasta")

    #Make an alignment and a phylogeny if there are at least two sequences    
    num_of_seqs = len(target_dict[target_class])    
    if num_of_seqs < 2:
        if num_of_seqs == 0: file_flag = "_empty.fas"
        elif num_of_seqs == 1: file_flag = "_single.fas" 
        target_alignment_file = target_alignment_file[:-4] + file_flag
        with open(target_alignment_file, "w") as outfile:
            SeqIO.write(target_dict[target_class], outfile, "fasta")
            
    else: 
        command_string = "clustalo -i " + target_sequence_file + " -o " + target_alignment_file + " --threads=8"    
        subprocess.call(command_string, shell=True)
        with open(log_file, "a") as log_name:
            log_name.write(str(datetime.datetime.now()))
            log_name.write(" Finished alignment of " + target_class + "\n")
    
        command_string = "fasttree -quiet -nt " + target_alignment_file + " > " + target_tree_file    
        subprocess.call(command_string, shell=True)
        with open(log_file, "a") as log_name:
            log_name.write(str(datetime.datetime.now()))
            log_name.write(" Finished tree of " + target_class + "\n")
    
    

        
        
      
                
        