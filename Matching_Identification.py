#! /usr/bin/env python
#packages that come with python
import pandas as pd
import numpy as np
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline
import argparse
import os
import sys
import string
import re

#Methods for Results
def ref_insertion_number(insert_str):
    return list(re.findall(r'\d+', insert_str))

def mini_rename(mini):
    progress = 0; MINI_LENGTH = len(mini)
    for i in range(len(mini)):
        old_name = mini.index[i]
        if old_name[-2:] != '.0':
            new_name = old_name[:-1] + '.0\n'
            mini = mini.rename(index={old_name: new_name})
        progress+=1
        percentage = (progress/MINI_LENGTH)*100
        print('Minimap2 Results Renaming Progress [%.5f%%]\r'%percentage, end="")
    print('\n')
    return mini

def error_type_matching(small_cell, large_cell):
    MAX_DIFFERENCE = 1
    matches = []
    print('jr')
    for i in range(len(large_cell)):
        if len(large_cell) > 0:
            large_error = int(large_cell[i])
            for j in range(len(small_cell)): 
                small_error = int(small_cell[j])
                if abs(large_error - small_error) <= MAX_DIFFERENCE: #area around each other 
                    if small_error not in matches: #no duplicates
                        matches.append(small_error)
    return matches  

def common_names(small,large):
    common = []
    indexs = large.index
    progress = 0; SMALL_LENGTH = len(small)
    for i in range(len(small)):
        name = small.index[i]
        if name in indexs:
            common.append(name)
        progress+=1
        percentage = (progress/SMALL_LENGTH)*100
        print('Finding Common Strands Progress [%.5f%%]\r'%percentage, end="")
    print('\n')
    return common

def matches(small,large,col):
    col_matches = []
    indexs = large.index
    for i in range(len(small)):
        name = small.index[i]
        if name in indexs:
            small_cell = ref_insertion_number(small.loc[name][col])
            large_cell = ref_insertion_number(large.iloc[i][col])
            col_matches.append(error_type_matching(small_cell,large_cell))
    return col_matches

def mismatches(all_matches):
    mismatch = []
    for rowNum in range(len(all_matches[5])):
        mismatch.append(len(all_matches[5][rowNum]) + len(all_matches[6][rowNum]) + len(all_matches[7][rowNum]) + len(all_matches[8][rowNum]) + len(all_matches[9][rowNum]) + len(all_matches[10][rowNum]))
    return mismatch
def all_counts(all_matches):
    counts = []
    nonHP_insertion = []; nonHP_deletion = []; subs = []; hp_insertion = []; hp_deletion = []; long_del = []
    for row in range(len(all_matches[5])):
        nonHP_insertion.append(len(all_matches[5][row])); nonHP_deletion.append(len(all_matches[6][row]))
        subs.append(len(all_matches[7][row])); hp_insertion.append(len(all_matches[8][row])); hp_deletion.append(len(all_matches[9][row])); long_del.append(len(all_matches[10][row]))
    counts.append(nonHP_insertion); counts.append(nonHP_deletion); counts.append(subs); counts.append(hp_insertion); counts.append(hp_deletion); counts.append(long_del)
    return counts
#******************************* Arguments **************************************
#1. add consensus positions and counts
#2. combine two inputs and put this ouput in back?
parser = argparse.ArgumentParser()
parser.add_argument("Minimap2ResultFile1", help="File path for CSV with the minimap2 Results using Minimap2_Error_Detection.py.")
parser.add_argument("Minimap2ResultFile2", help="File path for CSV with the minimap2 Results using Minimap2_Error_Detection.py. The results are from pbdagcon consensus generation")
parser.add_argument("OutputCsv", nargs='?', const = 1, default="Match_Combined_Results.csv", help="The msa sites compiled into a csv file.")

args = parser.parse_args()
mini1 = args.Minimap2ResultFile1
mini2 = args.Minimap2ResultFile2
output = args.OutputCsv

#******************************* New Index **************************************
# NEED faster renaming scheme
# aligned_file = AlignIO.read(msa,'clustal')
mini1_df = pd.read_csv(mini1)
mini1_df = mini1_df.set_index("Name")
#mini1_df = mini_rename(mini1_df)
mini2_df = pd.read_csv(mini2)
mini2_df = mini2_df.set_index("Name")
#mini2_df = mini_rename(mini2_df)

smaller = ""; larger = ""
if len(mini2_df) >= len(mini1_df):
    smaller = mini1_df
    larger = mini2_df
else:
    smaller = mini2_df
    larger = mini1_df
#+- a certain position with all the rows

#******************************* Same Strand IDs **************************************
names = common_names(smaller,larger)
all_matches = []
columns = ["Consensus NonHP Insertion Positions", "Consensus NonHP Deletion Positions", "Consensus Substitution Positions", "Consensus HP Insertion Positions", "Consensus HP Deletion Positions", "Reference NonHP Insertion Positions", "Reference NonHP Deletion Positions", "Reference Substitution Positions", "Reference HP Insertion Positions", "Reference HP Deletion Positions", "Reference Long Deletion Sites"]

#******************************* Same Strand Information **************************************
progress = 0; COL_LENGTH = len(columns)
for col in columns: #for each column in output
    all_matches.append(matches(smaller,larger,col))
    progress+=1
    percentage = (progress/COL_LENGTH)*100
    print('Matching Progress [%.5f%%]\r'%percentage, end="")

#******************************* Same Strand Summarize Information **************************************
mismatch = mismatches(all_matches)
counts = all_counts(all_matches)

#******************************* Add Summary Information **************************************
#df.insert(column position, column name, data,True)
info = {"Name":names, 
        "Mistmatches Similar from Reference":mismatch,
        "NonHP Insertion Similar from Reference":counts[0],
        "NonHP Deletion Similar from Reference":counts[1],
        "Substitution Similar from Reference":counts[2],
        "HP Insertion Similar from Reference":counts[3],
        "HP Deletion Similar from Reference":counts[4],
        "Long Deletion Similar from Reference":counts[5],
        "Consensus NonHP Insertion Positions":all_matches[0], 
        "Consensus NonHP Deletion Positions":all_matches[1],
        "Consensus Substitution Positions":all_matches[2], 
        "Consensus HP Insertion Positions":all_matches[3],
        "Consensus HP Deletion Positions":all_matches[4],
        "Reference NonHP Insertion Positions":all_matches[5], 
        "Reference NonHP Deletion Positions":all_matches[6],
        "Reference Substitution Positions":all_matches[7], 
        "Reference HP Insertion Positions":all_matches[8],
        "Reference HP Deletion Positions":all_matches[9],
        "Reference Long Deletion Positions":all_matches[10]}

#******************************* Concat Old Information **************************************
df = pd.DataFrame(info)
df = df.set_index("Name")

identify = 0; identity = ['_sparc','_pbd']
for comp_df in [smaller, larger]:
    del comp_df['Unnamed: 0']
    #comp_df = comp_df.reset_index()
    for column in comp_df:
        comp_df.rename(columns={column:column+identity[identify]}, inplace=True)
    identify+=1

#inner join based on names
df_final = pd.concat([df,smaller,larger], axis=1, join = "inner")

#delete unecessarry columns and reset indexs to row number
df_final = df_final.reset_index()
df_final.dropna(how='all', axis=1)
df_final.insert(19,'Name_sparc',names,allow_duplicates = False)
df_final.insert(44,'Name_pbd',names,allow_duplicates = False)

print(df_final)
df_final.to_csv(output, index=False)