#!/usr/bin/env python3
import pandas as pd
import numpy as np
from matplotlib import cm as cm
from matplotlib import pyplot as plt
import seaborn as sns
import matplotlib
import os

'''
This script takes files called Orthogroups.GeneCount.tsv that are generated by
Orthofinder and extracts a count for the number of orthogroups shared by each pair of
genomes. Then it calculates a gene sharing index for each pair which is the fraction
of orthogroups shared, averaged for the two genomes. If the same genome pair appears in
more than one input file, the index is averaged.

INPUT FILES (in the data folder)
One or more Orthogroups.GeneCount.tsv files, with unique prefixes

REFERENCE FILES (in the reference folder)
names_in_itol_order.txt A list of the genomes, in the desired output order
strain_info_440.csv A file that specifies the genospecies of each genome
gs_colour_list.csv  A file that specifies the colour to use for each genospecies

OUTPUT FILES (will appear in the output folder)
gene_sharing_pairwise_list.csv  List of pairwise genome combinations and their gene sharing index
gene_sharing_square_matrix.csv Square pivot table of gene sharing index
gene_sharing_square_plot.pdf   Square plot of gene sharing index
'''

#output files
pairwise_file = "../output/gene_sharing_pairwise_list.csv"
square_file = "../output/gene_sharing_square_matrix.csv"
figure_file = "../output/gene_sharing_square_plot.pdf"

#reference files
itol_order_file = "../reference/names_in_itol_order.txt"
strain_info_file = "../reference/strain_info_440.csv"
gs_colour_file = "../reference/gs_colour_list.csv"

#get the list of input files
def list_files(directory, name_end):
    return (directory +"/"+ f for f in os.listdir(directory) if f.endswith(name_end))

input_file_list =  list(list_files("../data","GeneCount.tsv"))

#create df_all to hold concatenated results from all inputs
df_all  = pd.DataFrame()

for infile in input_file_list:
    #open the data in a Pandas dataframe and drop the first and last columns
    df = pd.read_table(infile)
    df.pop('Orthogroup')
    df.pop('Total')

    #filter out orthogroups that have more than two copies in any genome
    df = df[df.max(axis=1) <= 2]

    #create a list with pairwise shared gene counts [[strain1, strain2, count],...]
    pair_list = []
    for col1 in df:
        for col2 in df:
            #add a temporary column 'Shared' that is 1 if both have the gene, 0 otherwise
            df['Shared'] = np.sign(df[col1]*df[col2])
            shared = df['Shared'].sum()
            #transfer the number of shared genes to pair_list
            if col2 != 'Shared':
                pair_list.append([col1,col2,shared])

    #make pair_list into a dataframe
    df2 = pd.DataFrame(pair_list, columns = ['Str1', 'Str2', 'Shared'])
    df2 = df2.replace(0, np.nan)

    #make a dictionary of total genes per strain
    df2_self = df2[df2['Str1'] == df2['Str2']]
    df2_dict = df2_self.set_index('Str1').to_dict()
    tot_dict = df2_dict['Shared']

    def tot_gene(x):
        return tot_dict[x]
    
    df2 = df2[df2.Str1.isin(tot_dict) & df2.Str2.isin(tot_dict)]

    #add a column to df2 that is Shared / mean total genes for each strain pair
    df2['Str1_tot'] = df2['Str1'].apply(lambda x : tot_gene(x))
    df2['Str2_tot'] = df2['Str2'].apply(lambda x : tot_gene(x))
    df2['Shared_norm'] = df2['Shared']/(df2['Str1_tot']*2) + df2['Shared']/(df2['Str2_tot']*2)
    
    df_all = df_all.append(df2)

#if a strain pair appears more than once, take the mean
df_all = df_all.groupby(['Str1','Str2'], as_index = False)[['Shared_norm']].mean()

#Make a sorter dict to convert strain name to rank in the sorted order
with open(itol_order_file) as tree_order:
    sorter = list(tree_order)
sorter = list(map(lambda x: x.strip(), sorter))

sorter_index = dict(zip(sorter, range(len(sorter))))
inverse_sorter_index  = dict([(value, key) for key, value in sorter_index.items()])

def sort_order(x):
    return sorter_index[x]
    
#Add columns of ranks for sorting, then sort
df_all['Strain1'] = df_all['Str1'].apply(lambda x : sort_order(x))
df_all['Strain2'] = df_all['Str2'].apply(lambda x : sort_order(x))

df_all = df_all.sort_values(by = ['Strain1','Strain2'])

#output the pairwise list
df_all.to_csv(pairwise_file)

#pivot the table from long form to square
#NB pivot sorts by the labels, so it is essential to have ranks not names here
dfp_rank = df_all.pivot(index='Strain1', columns='Strain2', values='Shared_norm')

#put strain names back in place of ranks
dfp = dfp_rank.rename(columns=inverse_sorter_index, index=inverse_sorter_index)

#output the square table
dfp.to_csv(square_file)

#get the sort order from dfp 
sorted_strains = list(dfp)

#Make a dict to map strain to genospecies
with open(strain_info_file) as strain_info:
    strain2gs = {}
    for strain_gs in strain_info:
        (strain, gs) = strain_gs.strip().split(',')
        strain2gs[strain] = gs

#Make a dict to map genospecies to colour        
with open(gs_colour_file) as gs_colours:
    gs2colour = {}
    for gs_col in gs_colours:
        (gs,colour) = gs_col.strip().split(',')
        gs2colour[gs] = colour

#define the order of colours in the genospecies stripe        
colour_list = [gs2colour[strain2gs[strain]] for strain in sorted_strains]       

#make the background grey (who knows why this works)
sns.set(font_scale=1)

#plot the ANI table with a blue-yellow-red scheme
cmap = cm.get_cmap('jet', 300)

ax = plt.figure(figsize=(10, 10))

ax = sns.clustermap(
    dfp, 
    cmap = cmap, 
    row_cluster=False, col_cluster=False,
    row_colors = colour_list,
    xticklabels=False,
    yticklabels=False   
)

#add a key to genospecies
patches1 = [
    matplotlib.patches.Patch(color=color, label=label)
    for label, color in zip(gs2colour.keys(),gs2colour.values())]
legend1 = plt.legend(handles=patches1, loc=(-.2,-3))
ax = plt.gca().add_artist(legend1)

plt.savefig(figure_file)

