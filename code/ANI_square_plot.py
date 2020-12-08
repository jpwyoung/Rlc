import pandas as pd
import numpy as np

import pylab as pl
import time
from matplotlib import cm as cm
from matplotlib import pyplot as plt
import seaborn as sns
import matplotlib

'''
This script plots a square matrix of all-by-all ANI values, sorted in a defined order, 
with a coloured band showing genospecies down the side.
'''

#read in the ANI table
df = pd.read_csv('../data/445_ANI.tab', names = ['Query_name', 'Reference_name', 'ANI', 'Fragments', 'Matches'], sep="\t")

df['Query_name'] = df['Query_name'].str.strip('.fna')
df['Reference_name'] = df['Reference_name'].str.strip('.fna')

#Make a sorter dict to convert strain name to rank in the sorted order
with open('../reference/names_in_itol_order.txt') as tree_order:
    sorter = list(tree_order)
sorter = list(map(lambda x: x.strip(), sorter))

sorterIndex = dict(zip(sorter, range(len(sorter))))
inverse_sorterIndex  = dict([(value, key) for key, value in sorterIndex.items()])

def sort_order(x):
    return sorterIndex[x]

#Restrict the DataFrame to the strains of interest    
df = df[df.Query_name.isin(sorter) & df.Reference_name.isin(sorter)]

#Add columns of ranks for sorting, then sort
df['Query'] = df['Query_name'].apply(lambda x : sort_order(x))
df['Reference'] = df['Reference_name'].apply(lambda x : sort_order(x))

dfs = df.sort_values(by = ['Reference','Query'])

#pivot the table from long form to square
dfp = dfs.pivot(index='Query', columns='Reference', values='ANI')
dfp = dfp.rename(columns=inverse_sorterIndex, index=inverse_sorterIndex)

#output the square table
#dfp.to_csv('/Users/peter/Documents/Manuscripts/R_leg_genospecies/data/leg-anh_2020-09-18/ANI/ANI_pivot_01.csv')

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

#define the order of colours in the genospecies stripe        
colour_list = [gs2colour[strain2gs[strain]] for strain in sorter]       

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

plt.savefig('../output/ANI_square.pdf')
plt.savefig('../output/ANI_square_dpi150.png',dpi=150)


#recode the ANI values as -4, -2, +4 for <95, 95-96, >96
dfp_cut = dfp.apply(lambda x : (3*np.sign(x-96) + np.sign(x-95)))

#plot the thresholded ANI table with pale-grey-black
ax = plt.figure(figsize=(10, 10))

ax = sns.clustermap(
    dfp_cut, 
    cmap = sns.light_palette("#000000", as_cmap=True),
    row_cluster=False, col_cluster=False, 
    row_colors = colour_list, 
    xticklabels=False,
    yticklabels=False
)

patches1 = [
    matplotlib.patches.Patch(color=color, label=label)
    for label, color in zip(gs2colour.keys(),gs2colour.values())]
legend1 = plt.legend(handles=patches1, loc=(-.2,-3))
ax = plt.gca().add_artist(legend1)

plt.savefig('../output/ANI_square_b+w.pdf')
plt.savefig('../output/ANI_square_b+w_dpi150.png',dpi=150)