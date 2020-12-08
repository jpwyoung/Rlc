#!/usr/bin/env python3
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

'''
This script takes the gene_sharing_pairwise_list.csv generated by gene_sharing_tables+plot.py
and the 445_ANI.tab generated for 445* genomes by ANI_square_plot.py and draws a scatterplot
of gene sharing v. ANI using all pairwise values of gene sharing.

*all 440 Rlc+anhuiense plus a few duplicates
'''

figure_file = "../output/gene_sharing_v_ANI.pdf"

#read in the ANI table
dfANI = pd.read_csv('../data/445_ANI.tab', names = ['Query_name', 'Reference_name', 'ANI', 'Fragments', 'Matches'], sep="\t")

dfANI['Query_name'] = dfANI['Query_name'].str.strip('.fna')
dfANI['Reference_name'] = dfANI['Reference_name'].str.strip('.fna')

#read in the shared genes table
df2 = pd.read_csv('../output/gene_sharing_pairwise_list.csv')

#merge the orthogroup and ANI dataframes
df_merged = pd.merge(df2, dfANI, how='left', left_on=['Str1', 'Str2'], right_on = ['Query_name','Reference_name'])

df_merged = df_merged.rename(columns = {'Shared_norm':'Gene sharing index'}) 

sns.set(font_scale=2)
sns.set_style("white")
sns.set_style("ticks")

ax = plt.figure(figsize=(10, 10))

ax = sns.scatterplot(
    data = df_merged,
    x='ANI', 
    y='Gene sharing index', 
    marker='.', 
    alpha = 0.3) 

plt.savefig(figure_file)
