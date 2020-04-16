import numpy as np
import pandas as pd
import pickle
import os.path
from subprocess import call


def ReadEssential_NonEssential_Genes():
    essential_genes = pd.read_csv('./Hart/Essential_Genes.txt', sep='\t', header=0).drop(columns=['HGNC_ID'])
    nonessential_genes = pd.read_csv('./Hart/Nonessential_Genes.txt', sep='\t', header=0).drop(columns=['HGNC_ID', 'ENTREZ_ID'])
    pickle.dump( (essential_genes, nonessential_genes), open( "./pickled/essential_nonessential_genes.p", "wb" ) )
    return


def returnNonZeros(df):
    checkZero = df.iloc[:,2:].copy()
    checkZero['sum'] = checkZero.sum(axis=1)
    df = df[checkZero.loc[:,'sum']>0]
    return(df)

def ProcessTKODataForDESEQ():
    # Code to compare the simulation data to the Hart TKO data

    # Process: readcount-HCT116_1-lib1
    filename = 'HCT116_1-lib1'
    print('starting to read %s'%filename)
    read_counts = pd.read_csv('./data/readcount-%s'%filename, sep='\t', header=0)
    read_counts = read_counts.loc[:,['GENE', 'GENE_CLONE', 'LIB1_T0','LIB1_T18_A', 'LIB1_T18_B']]
    read_counts.rename(columns={'LIB1_T0':'T0', 'LIB1_T18_A':'R1_F', 'LIB1_T18_B':'R2_F'}, inplace=True)
    read_counts = returnNonZeros(read_counts)
    read_counts.to_csv('./input/%s.csv'%filename, sep=',', index=False)
    call('rscript ./R_Scripts/DESEQ2_TKO_2Rep.R %s' %filename, shell = True)

    # Process: readcount-HCT116_2-lib1
    filename = 'HCT116_2-lib1'
    print('starting to read %s'%filename)
    read_counts = pd.read_csv('./data/readcount-%s'%filename, sep='\t', header=0)
    read_counts = read_counts.loc[:,['GENE', 'GENE_CLONE', 'T0','T18A', 'T18B', 'T18C']]
    read_counts.rename(columns={'T18A':'R1_F', 'T18B':'R2_F', 'T18C':'R3_F'}, inplace=True)
    read_counts = returnNonZeros(read_counts)
    read_counts.to_csv('./input/%s.csv'%filename, sep=',', index=False)
    call('rscript ./R_Scripts/DESEQ2_TKO_3Rep.R %s' %filename, shell = True)

    # Process: readcount-DLD1-lib1
    filename = 'DLD1-lib1'
    print('starting to read %s'%filename)
    read_counts = pd.read_csv('./data/readcount-%s'%filename, sep='\t', header=0)
    read_counts = read_counts.loc[:,['GENE', 'GENE_CLONE', 'DLD_T0','DLD_ETOH_R1', 'DLD_ETOH_R2', 'DLD_ETOH_R3']]
    read_counts.rename(columns={'DLD_T0':'T0', 'DLD_ETOH_R1':'R1_F', 'DLD_ETOH_R2':'R2_F', 'DLD_ETOH_R3':'R3_F'}, inplace=True)
    read_counts = returnNonZeros(read_counts)
    read_counts.to_csv('./input/%s.csv'%filename, sep=',', index=False)
    call('rscript ./R_Scripts/DESEQ2_TKO_3Rep.R %s' %filename, shell = True)

    # Process: readcount-DLD1-lib1
    filename = 'HeLa-lib1'
    print('starting to read %s'%filename)
    read_counts = pd.read_csv('./data/readcount-%s'%filename, sep='\t', header=0)
    read_counts = read_counts.loc[:,['GENE', 'GENE_CLONE', 'T0','T18A', 'T18B', 'T18C']]
    read_counts.rename(columns={'T18A':'R1_F', 'T18B':'R2_F', 'T18C':'R3_F'}, inplace=True)
    read_counts = returnNonZeros(read_counts)
    read_counts.to_csv('./input/%s.csv'%filename, sep=',', index=False)
    call('rscript ./R_Scripts/DESEQ2_TKO_3Rep.R %s' %filename, shell = True)

    # Process: readcount-RPE1-lib1
    filename = 'RPE1-lib1'
    print('starting to read %s'%filename)
    read_counts = pd.read_csv('./data/readcount-%s'%filename, sep='\t', header=0)
    read_counts = read_counts.loc[:,['GENE', 'GENE_CLONE', 'T0','T18A', 'T18B']]
    read_counts.rename(columns={'T18A':'R1_F', 'T18B':'R2_F'}, inplace=True)
    read_counts = returnNonZeros(read_counts)
    read_counts.to_csv('./input/%s.csv'%filename, sep=',', index=False)
    call('rscript ./R_Scripts/DESEQ2_TKO_2Rep.R %s' %filename, shell = True)

    # Process: readcount-GBM-lib1
    filename = 'GBM-lib1'
    print('starting to read %s'%filename)
    read_counts = pd.read_csv('./data/readcount-%s'%filename, sep='\t', header=0)
    read_counts = read_counts.loc[:,['GENE', 'GENE_CLONE', 'T0','T21A', 'T21B']]
    read_counts.rename(columns={'T21A':'R1_F', 'T21B':'R2_F'}, inplace=True)
    read_counts = returnNonZeros(read_counts)
    read_counts.to_csv('./input/%s.csv'%filename, sep=',', index=False)
    call('rscript ./R_Scripts/DESEQ2_TKO_2Rep.R %s' %filename, shell = True)


    # essential_genes = pd.read_csv('./Hart/Essential_Genes.txt', sep='\t', header=0).drop(columns=['HGNC_ID'])
    # genesInData= set(read_counts.loc[:,'GENE'])
    # essGenes = set(essential_genes.loc[:,'Gene'])
    # print(len(essGenes))
    # print(essGenes.difference(genesInData))









