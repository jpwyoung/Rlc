#!/usr/bin/env python3
import pandas as pd
import numpy as np
from matplotlib import cm as cm
from matplotlib import pyplot as plt
import seaborn as sns
import matplotlib

'''
This script calculates ANI separately for chromosomal and nonchromosomal scaffolds.
It takes output from fastANI run on the separated chr and non sequence data files that 
were created by chr_non_orthocore.py, using the genospecies representatives as reference.
Output is a scatterplot and a facetgrid of plots of non v. chr ANI with all strains as queries and 
genospecies representatives as reference.

N.B. When I upgraded to seaborn 0.11.0, it stopped plotting in the last facet. To get the 
full set of 18 Rlc plots, I added R. anhuiense as a 19th plot (which is displayed as an 
empty facet that can be edited out of the figure).
'''

figure_file = "../output/chrom_nonchrom/chr_v_nonchr_ANI_plot.pdf"
grid_file = "../output/chrom_nonchrom/chr_v_nonchr_ANI_plot.pdf"

#read in the chr_ANI table
dfANIC = pd.read_csv('../output/chrom_nonchrom/chr/all_type_strain_ANI.tab', names = ['Query_name_C', 'Reference_name_C', 'ANI_C', 'Fragments_C', 'Matches_C'], sep="\t")

#read in the non_ANI table
dfANIN = pd.read_csv('../output/chrom_nonchrom/non/all_type_strain_ANI.tab', names = ['Query_name_N', 'Reference_name_N', 'ANI_N', 'Fragments_N', 'Matches_N'], sep="\t")

#merge the two ANI dataframes
df_merged = pd.merge(dfANIC, dfANIN, how='left', left_on=['Query_name_C','Reference_name_C'], right_on = ['Query_name_N','Reference_name_N'])

###
#plot a scatterplot with all data
###

sns.set(font_scale=2)
sns.set_style("white")
sns.set_style("ticks")

X_plot = np.linspace(91, 100, 100)
Y_plot = X_plot

ax = plt.figure(figsize=(10, 10))

ax = sns.scatterplot(
    data = df_merged,
    x='ANI_C', 
    y='ANI_N', 
    marker='.', 
    alpha = 0.3) 

plt.plot(X_plot, Y_plot, color='gray')

plt.savefig(figure_file)

###
#plot a grid of scatterplots, one for each genospecies representative
###

df_merged['Query_name_C'] = df_merged['Query_name_C'].str.strip('.fna')
df_merged['Reference_name_C'] = df_merged['Reference_name_C'].str.strip('.fna')

#Get lists of the type strains and the genospecies, in the standard order    
with open('../reference/Type_strain_list_anh.tab') as type_list:    
    ref_list = []
    gs_list = []
    for ref_gs in type_list:
        ref_list.append(ref_gs[:ref_gs.find("\t")])
        gs_list.append(ref_gs[ref_gs.find("\t")+1:-1])
        
#Make a dict to map strain to genospecies
with open('../reference/strain_info_440.csv') as strain_info:
    strain2gs = {}
    for strain_gs in strain_info:
        (strain, gs) = strain_gs.strip().split(',')
        strain2gs[strain] = gs
        
#Make a dict to map genospecies to colour
with open('../reference/gs_colour_list.csv') as gs_colours:
    gs2colour = {}
    for gs_col in gs_colours:
        (gs,colour) = gs_col.strip().split(',')
        gs2colour[gs] = colour

def get_gs(strain):
    return strain2gs[strain]

#Add columns with the genospecies of the reference ('gs') and the query ('Genospecies'),
#and whether these are the same or not ('in_gs')
#These column names will appear in the figure
df_merged['gs'] = df_merged['Reference_name_C'].apply(lambda x: get_gs(x))
df_merged['Genospecies'] = df_merged['Query_name_C'].apply(lambda x: get_gs(x))
df_merged['in_gs'] = df_merged['Genospecies'] == df_merged['gs']

#Use Seaborn relplot to make the plot

sns.set(font_scale=1.5)
sns.set_style("white")
sns.set_style("ticks")


ax = plt.figure(figsize=(10, 10))

ax = sns.relplot(
    data=df_merged, 
    x="ANI_C", 
    y="ANI_N", 
    hue="Genospecies", 
    palette = gs2colour,
    edgecolor = None, 
    style = "in_gs",
    markers = {True:"o", False: "o"},
    aspect = 1,
    col="gs", 
    col_order = gs_list,
    col_wrap=6,
    legend = False)
    
#ax.map(plt.axvline, x=96, ls='--', c='gray')
#ax.map(plt.axhline, y=96, ls='--', c='gray')

plt.plot(X_plot, Y_plot, color='gray')
X_plot = np.linspace(91, 100, 100)
Y_plot = X_plot

plt.savefig(grid_file)     
