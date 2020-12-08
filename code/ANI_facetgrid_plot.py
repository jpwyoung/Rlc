#!/usr/bin/env python3
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns


#Create a DataFrame with the pairwise ANI values
df = pd.read_csv('../data/445_ANI.tab', names = ['Query_name', 'Reference_name', 'ANI', 'Fragments', 'Matches'], sep="\t")

df['Query_name'] = df['Query_name'].str.strip('.fna')
df['Reference_name'] = df['Reference_name'].str.strip('.fna')

#Get a list of the query names
with open('../reference/names_in_itol_order.txt') as tree_order:
    query_list = list(tree_order)
    query_list = list(map(lambda x: x.strip(), query_list))

#Get lists of the type strains and the genospecies, in the standard order    
with open('../reference/Type_strain_list.tab') as type_list:    
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

#Edit the DataFrame to just the values we are going to use  
df = df[df.Query_name.isin(query_list) & df.Reference_name.isin(ref_list)]

#Add a column 'rank' that ranks the queries by ANI for each reference separately
df['rank'] = df.groupby('Reference_name')['ANI'].rank(method='first',ascending=False)

#Add columns with the genospecies of the reference ('gs') and the query ('Genospecies'),
#and whether these are the same or not ('in_gs')
#These column names will appear in the figure
df['gs'] = df['Reference_name'].apply(lambda x: get_gs(x))
df['Genospecies'] = df['Query_name'].apply(lambda x: get_gs(x))
df['in_gs'] = df['Genospecies'] == df['gs']

#Use Seaborn relplot to make the plot

sns.set(font_scale=1.5)
sns.set_style("white")
sns.set_style("ticks")

ax = plt.figure(figsize=(10, 10))

ax = sns.relplot(
    data=df, 
    x="rank", 
    y="ANI", 
    hue="Genospecies", 
    palette = gs2colour,
    edgecolor = None, 
    style = "in_gs",
    markers = {True:"o", False: "."},
    aspect = 0.5,
    col="gs", 
    col_order = gs_list,
    col_wrap=6,
    legend = False)
    
ax.map(plt.axhline, y=96, ls='--', c='gray')

plt.savefig('../output/ANI_grid.pdf')     
plt.savefig('../output/ANI_grid.png')     

