 #!/usr/bin/env python3
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os

'''
This script sorts scaffolds into chromosomal and nonchromosomal sets, based on whether they
carry any of the 3215 'orthocore' genes that we have found are normally encoded on the
chromosomal backbone in genomes of the Rhizobium leguminosarum species cluster.

It takes the output from find_genes.py, run using the orthocore genes as query, in hits_folder.
It requires access to the corresponding genome files in genomes_folder.
It creates two fasta files from each genome, with chromosomal scaffolds in chr_folder and 
other scaffolds in non_folder. Filenames are the same as the input genome filenames.
'''

genomes_folder = "../data/"

hits_folder = "../output/find_genes/top_genes_hits/"
chr_folder = "../output/chrom_nonchrom/chr/"
non_folder = "../output/chrom_nonchrom/non/"

os.mkdir("../output/chrom_nonchrom")
os.mkdir("../output/chrom_nonchrom/chr")
os.mkdir("../output/chrom_nonchrom/non")

#prepare to output a table with scaffold counts and sizes for each genome
summary_file = "../output/chrom_nonchrom/chr_non_summary.tab"
summary_list = []

for file in os.listdir(hits_folder):
    hits_file = os.path.join(hits_folder, file)
    
    #Make a dict of hits from the chr blast output
    #The values are the number of chromosomal genes identified in that scaffold
    chr_scaff = {}
    with open(hits_file) as blast_matches:
        for hit in blast_matches:
            hit_data = hit.rstrip("\n").split("\t")
            scaff = hit_data[1]
            if scaff in chr_scaff:
                chr_scaff[scaff] += 1
            else: 
                chr_scaff[scaff] = 1

    #find the corresponding genome file    
    genome_file = os.path.basename(hits_file).replace("_hits.tab",".fna")
    fasta_in = genomes_folder + genome_file
    
    if os.path.exists(fasta_in):
    
        #open the genome file and sort scaffolds into chr and non output files 
        all_scaffs_dict = SeqIO.to_dict(SeqIO.parse(fasta_in, "fasta"))
        chr_list = []
        non_list = []
        for scaff in all_scaffs_dict:
            if scaff in chr_scaff:
                chr_list.append(all_scaffs_dict[scaff])
            else:
                non_list.append(all_scaffs_dict[scaff])
        
        with open(chr_folder+genome_file, "w") as outfile:
            SeqIO.write(chr_list, outfile, "fasta")

        with open(non_folder+genome_file, "w") as outfile:
            SeqIO.write(non_list, outfile, "fasta")

        #count the chr and non scaffolds               
        chr_num_scaff = len(chr_list)
        non_num_scaff = len(non_list)
        
        #get the total length of chr and non scaffolds
        chr_length  = non_length = 0
        for scaff in chr_list:
            chr_length += len(scaff)        
        for scaff in non_list:
            non_length += len(scaff)        
        
        summary_list += [(genome_file, chr_num_scaff, non_num_scaff, chr_length, non_length)]
        
#write the summary file
with open(summary_file, "w") as outfile:
        
    for genome_data in summary_list:
        for item in genome_data:            
            outfile.write(str(item) + "\t")
        outfile.write("\n")
       
                    

     
