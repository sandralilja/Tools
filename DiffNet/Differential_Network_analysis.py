#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 09:01:24 2019

@author: sanli71
"""

import argparse
import pandas as pd
import os
import shutil
import numpy as np
import statsmodels.stats.multitest as multitest
from datetime import datetime

#os.chdir(os.getcwd() + '/script')
def define_dirs(perm_min, perm, timepoint):
    dir_in = os.getcwd() + '/../DGE_data/CellPhoneDB_indata/'
    cp_indir = os.getcwd() + '/../results/CellPhoneDB_out/tmp_in' + str(perm_min) + 'to' + str(perm) + '_' + timepoint + '/'
    cp_outdir = os.getcwd() + '/../results/CellPhoneDB_out/tmp_out' + str(perm_min) + 'to' + str(perm) + '_' + timepoint + '/'
    dir_out = os.getcwd() + '/../results/CellPhoneDB_DiffNet_out/'
    if os.path.exists(cp_indir) == False:
        os.mkdir(cp_indir)
    if os.path.exists(cp_outdir) == False:
        os.mkdir(cp_outdir)
    if os.path.exists(dir_out) == False:
        os.mkdir(dir_out)
    return dir_in, cp_indir, cp_outdir, dir_out

def get_data(samp1, samp2, dir_in):
    if samp1.split('_')[2] == samp2.split('_')[2]:
        group = samp1.split('_')[3]
        timepoint = samp1.split('_')[4]
        timepoint = timepoint.split('.')[0]
        Group1_data = pd.read_table(dir_in + samp1, index_col = 0)
        Group2_data = pd.read_table(dir_in + samp2, index_col = 0)
    return Group1_data, Group2_data, group, timepoint 

def make_meta(samp):
    Cell = samp.columns
    meta = pd.DataFrame(columns=['Cell', 'cell_type'])
    for row in range(0, len(Cell)):
        df = pd.DataFrame([[Cell[row], 
                          Cell[row].split('_')[0]]], 
                          columns=['Cell', 'cell_type'])
        meta = meta.append(df)
    return meta

def cellphone_analysis(data, cp_indir, cp_outdir, save=False, sampname='random', permute=True):
    meta = make_meta(data)
    # save the files for cellphone analysis
    data.to_csv(cp_indir + 'expression_matrix.txt', sep = '\t')
    meta.to_csv(cp_indir + 'meta_data.txt', sep = '\t', index = False)
    # Run CellPhoneDB    
    if permute == True:
        iterations = 1000
    if permute == False:
        iterations = 1
    os.system('./CellPhoneDB_analysis.sh %s %s %s'%(cp_indir, cp_outdir, iterations))
    # Read in significant interactions
    if permute == True:
        interactions = pd.read_table(cp_outdir + 'significant_means.txt')
    if permute == False:
        interactions = pd.read_table(cp_outdir + 'means.txt')    
    if save == True:
        savedir = cp_outdir + '../real_networks/' + sampname
        shutil.copytree(cp_outdir, savedir)         
    filesToRemove = [os.path.join(cp_outdir,f) for f in os.listdir(cp_outdir)]
    for f in filesToRemove:
        os.remove(f) 
    filesToRemove = [os.path.join(cp_indir,f) for f in os.listdir(cp_indir)]
    for f in filesToRemove:
        os.remove(f) 
    return interactions       

def subset_interaction_tables(net):
    intcols = [col for col in net.columns if '|' in col]
    int_net = net[intcols]
    int_net.index = net['id_cp_interaction']
    return int_net    

def calculate_scores(net, data):
    # Count total number of cells
    Nall = len(data.columns)
    # Count number of cells per cell type
    uniqueCells = []
    allCells = []
    for cell in list(data.columns):
        ct = cell.split('_')[0]
        allCells.append(ct)
        if ct not in uniqueCells:
            uniqueCells.append(ct)
    NperCellType = pd.DataFrame(columns = ['celltype', 'count'], index = range(0, len(uniqueCells)))
    for cellu in range(0, len(uniqueCells)): 
        NperCellType.loc[cellu] = uniqueCells[cellu], allCells.count(uniqueCells[cellu]) 
    del uniqueCells, allCells
    # calculate the scores based on (N * mij) / |Nall|
    int_net = subset_interaction_tables(net)
    int_net = int_net.fillna(0)
    for col in range(len(int_net.columns)):
        intcells = int_net.columns[col]
        intcell1 = intcells.split('|')[0]
        intcell2 = intcells.split('|')[1]
        if intcell1 == intcell2:
            N = NperCellType.loc[NperCellType['celltype'] == intcell1].iloc[0]['count']
        elif intcell1 != intcell2:
            n1 = NperCellType.loc[NperCellType['celltype'] == intcell1].iloc[0]['count']
            n2 = NperCellType.loc[NperCellType['celltype'] == intcell2].iloc[0]['count']
            N = n1 + n2
        for row in range(len(int_net)):
            int_net.iloc[row,col] = (int_net.iloc[row,col] * N)/Nall
    return int_net

def net_diff(score_net1, score_net2):
    ''' 
    The function calculates the difference between the networks as;
    score_net1 - score_net2
    
    The reason for filling the networks is so that they include the same rows 
    and columns, so that the difference will be calculated correctly, not 
    producing NaNs.
    '''
    # Add missing columns to make network dimentions equal
    score_net1 = pd.concat([score_net1, pd.DataFrame(columns = list(set(score_net2.columns).difference(set(score_net1.columns))))], sort = True)
    score_net2 = pd.concat([score_net2, pd.DataFrame(columns = list(set(score_net1.columns).difference(set(score_net2.columns))))], sort = True)
    # Add missing rows to make network dimentions equal
    score_net1 = pd.concat([score_net1, pd.DataFrame(index = list(set(score_net2.index).difference(set(score_net1.index))))], sort = True)
    score_net2 = pd.concat([score_net2, pd.DataFrame(index = list(set(score_net1.index).difference(set(score_net2.index))))], sort = True)
    # Change nan to zeros
    score_net1 = score_net1.fillna(0)
    score_net2 = score_net2.fillna(0)
    # calculate the network difference (score_net1 - score_net2)
    DiffNet = score_net1.sub(score_net2)
    return DiffNet

def compare_DiffNets(diffnet_1, diffnet_2, corr_df):    
    # Add any potential missing rows and columns to diffnet_2 so that all interactions in diffnet_1 can be compared
    diffnet_2 = pd.concat([diffnet_2, pd.DataFrame(index = list(set(diffnet_1.index).difference(set(diffnet_2.index))))], sort = True)
    diffnet_2 = pd.concat([diffnet_2, pd.DataFrame(columns = list(set(diffnet_1.columns).difference(set(diffnet_2.columns))))], sort = True)
    # Change nan to zeros in diffnet_2
    diffnet_2  = diffnet_2.fillna(0)
    # sort and subset diffnet_2 so that is has the same dimentions and order of rows and comulns as diffnet_1
    diffnet_2 = diffnet_2.loc[list(diffnet_1.index)]
    diffnet_2 = diffnet_2[list(diffnet_1.columns)]
    if (diffnet_1.index.equals(diffnet_2.index) == False):
        raise ValueError('The DiffNet indexes are not equal')
    if (diffnet_1.columns.equals(diffnet_2.columns) == False):
        raise ValueError('The DiffNet columns are not equal')
    if (diffnet_1.index.equals(corr_df.index) == False):
        raise ValueError('The corr_df indexes are not equal to the diffnet indexes')
    if (diffnet_1.columns.equals(corr_df.columns) == False):
        raise ValueError('The corr_df columns are not equal to the diffnet columns')
    # Check for which interactions |diffnet1[i,j]| > |diffnet2[i,j]|
    for col in range(len(diffnet_1.columns)):
        for row in range(len(diffnet_1.index)):
            if abs(diffnet_1.iloc[row,col]) > abs(diffnet_2.iloc[row,col]):
                corr_df.iloc[row, col] += 1
    return corr_df
            
def DiffNet_permutation(data, groupsize_1, groupsize_2, seed, cp_indir, cp_outdir, real_diffnet, corr_df):
    '''
    1. Randomly assign the cells into the Group1 and Group2.
    2. Do CellPhone and DiffNet analysis as before
    3. Test if the |difference| is larger for the real network then for 
    this permutation-based network
    '''
    if groupsize_1 + groupsize_2 != len(data.columns):
        raise ValueError('The merging of the samples is incorrect, gives the wrong amount of cells')
    group1 = data.sample(n=groupsize_1, random_state=seed, axis = 1)
    group2 = data[list(set(data.columns).difference(set(group1.columns)))]
    if len(group1.columns) + len(group2.columns) != len(data.columns):
        raise ValueError('Something wrong with the subsampling, gives the wrong amount of cells into the different groups')
    # CellPhone analysis
    group1_interactions = cellphone_analysis(group1, cp_indir, cp_outdir, permute=False)
    group2_interactions = cellphone_analysis(group2, cp_indir, cp_outdir, permute=False)
    # Calculate the scores for each interaction network
    group1_scores = calculate_scores(group1_interactions, group1)
    group2_scores = calculate_scores(group2_interactions, group2)
    # Calculate the difference between the networks 
    # (group1 - group2)
    DiffNet = net_diff(group1_scores, group2_scores)
    corr_df = compare_DiffNets(real_diffnet, DiffNet, corr_df)
    return corr_df
            
def FDR_correction(pval_df):
    pval_array = np.array(pval_df)    
    mask = np.isfinite(pval_array)
    pval_corrected = np.full(pval_array.shape, np.nan)
    pval_corrected[mask] = multitest.multipletests(pval_array[mask], method='fdr_bh')[1]
    pval_corrected_df = pd.DataFrame(pval_corrected, index=pval_df.index, columns=pval_df.columns)
    return pval_corrected_df   

def merge_CellPhoneout_DiffNet(diffnet, CellPhone_net_1, CellPhone_net_2):
    if [col for col in CellPhone_net_1.columns if '|' not in col and 'rank' not in col] != [col for col in CellPhone_net_2.columns if '|' not in col and 'rank' not in col]:
        raise ValueError('Columns in CellPhone significant means are not the same between samples')
    CP_cols = [col for col in CellPhone_net_1.columns if '|' not in col and 'rank' not in col]
    CP_info_1 = CellPhone_net_1[CP_cols]
    CP_info_2 = CellPhone_net_2[CP_cols]
    CP_info = CP_info_1.merge(CP_info_2, how='outer')
    if len(set(set(CP_info_1.iloc[:,1]).union(set(CP_info_2.iloc[:,1]))).union(set(CP_info.iloc[:,1]))) != len(CP_info):
        raise ValueError('Merging issues: the id_cp_interaction union is not complete after merging')
    if len(CP_info) != len(diffnet):
        raise ValueError('Dataframes to merge are of different row lengths')
    diffnet['id_cp_interaction'] = diffnet.index
    df_out = CP_info.merge(diffnet, how='outer')    
    return df_out

    
def main():
    parser = argparse.ArgumentParser(description='Cunduct DiffNet analysis between group1 and group2 using CellPhoneDB analysis for identification of cellular interactions')
    parser.add_argument('--seed_from', help='the first seed to use for permutation analysis', type=int)
    parser.add_argument('--seed_to', help='the last seed to use for permutation analysis', type=int)
    parser.add_argument('--timepoint', help='the sampel timepoint to analyze', type=str)
    parser.add_argument('--group1', help='group1 specific string. eg AllergenChallenged. Must be unique part of the input file names.', type=str)
    parser.add_argument('--group2', help='group2 specific string. eg NonChallenged. Must be unique part of the input file names.', type=str)
    parser.add_argument('--skip_to_stat', dest='skip_to_stat', action='store_true', default=False, help='Use if the DiffNet basic output is already calculated and stored in their respective folders.')
    args = parser.parse_args()
    
    print('Initiate analysis')
    starttime = datetime.now()
    perm_min = args.seed_from
    perm = args.seed_to
    timepoint = args.timepoint
    group1 = args.group1
    group2 = args.group2    
#    timepoint = '0h'
#    group1 = 'HA'
#    group2 = 'HC'
#    perm_min = 1
#    perm = 2
    dir_in, cp_indir, cp_outdir, dir_out = define_dirs(perm_min, perm, timepoint)        
    samps = os.listdir(dir_in)   
    samps = list(filter(lambda x:timepoint in x, samps))
    if len(samps) != 2:
        raise ValueError('Number of samples should be two')        
    if samps[0].split('_')[0] == samps[1].split('_')[0]:
        sampname = samps[0].split('_')[0] + '_' + samps[0].split('_')[4]
        sampname  = sampname.split('.')[0]            
    elif samps[0].split('_')[0] != samps[1].split('_')[0]:
        sampname = [samps[1].split('_')[0], samps[0].split('_')[0]]
        sampname.sort()
        sampname = sampname[0] + sampname[1] + '_' + samps[0].split('_')[4]
        sampname  = sampname.split('.')[0]
    G1_sampname = samps[0].split('_')[0] + '_' + samps[0].split('_')[3] + '_' + samps[0].split('_')[4]
    G1_sampname  = G1_sampname.split('.')[0]            
    G2_sampname = samps[1].split('_')[0] + '_' + samps[1].split('_')[3] + '_' + samps[1].split('_')[4]
    G2_sampname  = G2_sampname.split('.')[0]            

    if args.skip_to_stat == False:
        G1s_mat = list(filter(lambda x:str(group1) in x, samps))
        G2s_mat = list(filter(lambda x:str(group2) in x, samps))
        if len(G1s_mat) != len(G2s_mat) != 1:
            raise ValueError('List of group1 and/or group2 files are not 1')
    
        samp = 0
        print('Load and prepare the data')
        Group1_data, Group2_data, group, timepoint = get_data(G1s_mat[samp], G2s_mat[samp], dir_in)
        print('Analysing samp from timepoint ' + timepoint)
        print('CellPhone analysis') 
        if perm_min == 1:    
            Group1_interactions = cellphone_analysis(Group1_data, cp_indir, cp_outdir, save=True, sampname=G1_sampname)
            Group2_interactions = cellphone_analysis(Group2_data, cp_indir, cp_outdir, save=True, sampname=G2_sampname)
        if perm_min != 1:    
            Group1_interactions = cellphone_analysis(Group1_data, cp_indir, cp_outdir)
            Group2_interactions = cellphone_analysis(Group2_data, cp_indir, cp_outdir)
        print('Calculate the scores for each interaction network')
        Group1_scores = calculate_scores(Group1_interactions, Group1_data)
        Group2_scores = calculate_scores(Group2_interactions, Group2_data)
        print('Calculate the difference between the real networks') 
        # (Group1 - Group2)
        AC_min_NC_DiffNet = net_diff(Group1_scores, Group2_scores)
        
        print('Permutation test to find the significant network differences')
        # Merge all data into one dataframe
        AllData = pd.merge(Group1_data, Group2_data, left_index=True, right_index=True, how='outer')
        AllData = AllData.fillna(0)   
        # Permutation test
        corr_df = pd.DataFrame(0, index=AC_min_NC_DiffNet.index, columns=AC_min_NC_DiffNet.columns)
        print('Total number of permutations: ' + str(perm))
        for z in range(perm_min,perm+1):
            print('Permutation seed: ' + str(z))
            corr_df = DiffNet_permutation(AllData, 
                                      groupsize_1=len(Group1_data.columns), 
                                      groupsize_2=len(Group2_data.columns), 
                                      seed=z, 
                                      cp_indir=cp_indir, 
                                      cp_outdir=cp_outdir, 
                                      real_diffnet=AC_min_NC_DiffNet, 
                                      corr_df=corr_df)
    
        print('Write basic output')
        corr_df.to_csv(dir_out + sampname + '_DiffNet_' + group1 + '_vs_' + group2 + '_Ntimes_realNet_higherthan_permNet_permutations_' + str(perm_min) + 'to' + str(perm) + '.csv', index=True)
        AC_min_NC_DiffNet.to_csv(dir_out + sampname + '_DiffNet_' + group1 + '_vs_' + group2 + '_s1_min_s2_realNets_permutations_' + str(perm_min) + 'to' + str(perm) + '.csv')

    if args.skip_to_stat == True:
        corr_df = pd.read_csv(dir_out + sampname + '_DiffNet_' + group1 + '_vs_' + group2 + '_Ntimes_realNet_higherthan_permNet_permutations_' + str(perm_min) + 'to' + str(perm) + '.csv', index_col=0)
        Group1_interactions = pd.read_table(cp_outdir + '../real_networks/' + G1_sampname + '/significant_means.txt')
        Group2_interactions = pd.read_table(cp_outdir + '../real_networks/' + G2_sampname + '/significant_means.txt')

    print('Calculate P-values')
    pval_df = pd.DataFrame(0, index=corr_df.index, columns=corr_df.columns)
    for col in range(len(corr_df.columns)):
        for row in range(len(corr_df.index)):
            pval_df.iloc[row,col] = (perm - corr_df.iloc[row,col]) / perm
    print('Perform FDR correction')
    pval_FDR_df = FDR_correction(pval_df)
    
    print('Produce the stat output dataframes')
    Nsign_preFDR = (pval_df < 0.05).astype(int).sum(axis=0)
    Nsign_postFDR = (pval_FDR_df < 0.05).astype(int).sum(axis=0)            
    pval_df = merge_CellPhoneout_DiffNet(pval_df, Group1_interactions, Group2_interactions)
    pval_FDR_df = merge_CellPhoneout_DiffNet(pval_FDR_df, Group1_interactions, Group2_interactions)
    print('Write stat output')
    pval_df.to_csv(dir_out + sampname + '_DiffNet_' + group1 + '_vs_' + group2 + '_pval_permutations_' + str(perm_min) + 'to' + str(perm) + '.csv', index=False)
    pval_FDR_df.to_csv(dir_out + sampname + '_DiffNet_' + group1 + '_vs_' + group2 + '_qval_permutations_' + str(perm_min) + 'to' + str(perm) + '.csv', index=False)
    Nsign_preFDR.to_csv(dir_out + sampname + '_DiffNet_' + group1 + '_vs_' + group2 + '_Nsignificant_interactions_preFDR_permutations_' + str(perm_min) + 'to' + str(perm) + '.csv', header = False)
    Nsign_postFDR.to_csv(dir_out + sampname + '_DiffNet_' + group1 + '_vs_' + group2 + '_Nsignificant_interactions_postFDR_permutations_' + str(perm_min) + 'to' + str(perm) + '.csv', header = False)

    os.rmdir(cp_indir)
    os.rmdir(cp_outdir)
    print(datetime.now() - starttime)

#    sign_preFDR = pval_df[pval_df < 0.05]
#    sign_preFDR.boxplot()
#    Nsign_preFDR = (pval_df < 0.05).astype(int).sum(axis=0)
#    Nsign_nonzero_preFDR = (sign_preFDR > 0).astype(int).sum(axis=0)
#    Nsign_zero_preFDR = (sign_preFDR == 0).astype(int).sum(axis=0)
#    sum(list(Nsign_nonzero_preFDR))
#    sum(list(Nsign_zero_preFDR))
#    

if __name__ == '__main__':
    main()           
            
         
