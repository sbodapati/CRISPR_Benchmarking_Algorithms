import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import seaborn as seas
import pandas as pd
from subprocess import call
# from Analysis import AnalyzeSim as AnalyzeSim

def RunSimulation(essential_genes_mean, non_essential_variance, essential_variance, num_essential_genes,
                  num_nonessential_genes, num_guides_per_gene, guide_efficiency_mean, guide_efficiency_variance,
                  num_control_guides, depth_factor, binding_efficiency, rep_num):

    # Creating Effect per gene
    nonEssential_effect = np.random.randn(num_nonessential_genes).reshape(-1,1)* non_essential_variance
    nonEssential_effect[nonEssential_effect<0] = 0
    nonEssential_effect[nonEssential_effect>1] = 1

    Essential_effect = np.random.randn(num_essential_genes).reshape(-1,1) * essential_variance + essential_genes_mean
    Essential_effect[Essential_effect<0.7] = 0.7
    # Essential_effect[Essential_effect>1] = 1 #Need to remove this to get wider distribution of the essential gene variance



    # seas.distplot(nonEssential_effect, label = 'nonEssential Effect', color = 'orange')
    # seas.distplot(Essential_effect, label = 'Essential Effect', color = 'blue')
    # plt.legend()
    # plt.title('Non-Essential Gene Effect vs Essential Gene Effect')


    # Creating variable number of guides per gene (uniform distribution between 5 and 16 guides)
    num_guides = np.ones(num_essential_genes+num_nonessential_genes)*num_guides_per_gene
    num_guides = num_guides.astype(int)
    # num_guides = vector of length total genes, but each element contains the number of guides the gene has.


    # Creating the various binding efficiencies per guide
    guide_efficiency = np.random.randn(np.sum(num_guides)) * guide_efficiency_variance + guide_efficiency_mean
    guide_efficiency[guide_efficiency>1] = 1
    guide_efficiency[guide_efficiency<0] = 0
    plt.figure()
    seas.distplot(guide_efficiency, label = 'guide efficiency', color = 'orange')
    plt.legend()
    plt.title('Distribution of guide efficiency')


    # guide_efficiency = vector of the total number of guides, with each element being the percent binding rate.


    numEssentialGuides = np.sum(num_guides[:num_essential_genes]) #number of guides for just essential genes.



    # Creating the total effect per guide (using gene effect and guide efficiency)
    # plt.figure()

    sg_Count_Names = range(0,sum(num_guides))
    sg_Count_Names = ['sg_'+ str(i+1) for i in sg_Count_Names]
    gene_Names = range(0,sum(num_guides))
    gene_Names = ['gene_'+ str(int(i/num_guides_per_gene) + 1) for i in gene_Names]


    essential_guide_efficiency_effect = guide_efficiency[:numEssentialGuides]
    index = 0
    for i in range(num_essential_genes):
        essential_guide_efficiency_effect[index:index + num_guides[i]] = Essential_effect[i] * essential_guide_efficiency_effect[index:index + num_guides[i]]
        index = index + num_guides[i]

    # plt.figure()
    # plt.hist(essential_guide_efficiency_effect,bins=10)
    essential_mask = np.random.uniform(0.0,1.0,(essential_guide_efficiency_effect.shape[0],1))
    essential_mask[essential_mask>=(1-binding_efficiency)] = 1
    essential_mask[essential_mask<(1-binding_efficiency)] = 0
    # print(essential_guide_efficiency_effect.shape)
    # print(essential_mask.shape)

    essential_guide_efficiency_effect = essential_guide_efficiency_effect[:,np.newaxis] * essential_mask
    # print(essential_guide_efficiency_effect.shape)
    # wait = input()
    # plt.figure()
    # plt.hist(essential_guide_efficiency_effect[:,0], bins=10)
    # plt.show()


    # print(essential_mask)
    # plt.figure()
    # plt.hist(essential_mask, range=(0.0,1.0))
    # plt.show()
    # wait = input()

    nonessential_guide_efficiency_effect = guide_efficiency[numEssentialGuides:]
    index = 0
    for i in range(num_nonessential_genes):
        nonessential_guide_efficiency_effect[index:index + num_guides[num_essential_genes + i]] = nonEssential_effect[i] * nonessential_guide_efficiency_effect[index:index + num_guides[num_essential_genes + i]]
        index = index + num_guides[i+num_essential_genes]
    nonessential_guide_efficiency_effect = nonessential_guide_efficiency_effect[:,np.newaxis]


    # seas.distplot(essential_guide_efficiency_effect, label = 'essential effect', color = 'blue')
    # seas.distplot(nonessential_guide_efficiency_effect, label = 'nonessential effect', color = 'orange')
    # plt.title('Final Effect with nonessential and essential genes')
    # plt.legend()


    # Creating Logfold change per guide without sequencing depth
    lfc_nonessential = nonessential_guide_efficiency_effect
    lfc_essential = essential_guide_efficiency_effect * -1



    # plt.figure()
    # seas.distplot(lfc_nonessential, label = 'nonessential effect', color = 'orange')
    # seas.distplot(lfc_essential, label = 'essential effect', color = 'blue')
    # plt.title('Logfold effect change with nonessential and essential genes')
    # plt.legend()


    # new sequencing

    # mean normalizing

    onesAndZeros = np.concatenate((np.ones(lfc_essential.shape[0]),np.zeros(lfc_nonessential.shape[0])), axis=0)[:,np.newaxis]

    fullGuideVec = np.concatenate((lfc_essential,lfc_nonessential), axis=0)
    mean_of_fullGuideVec = np.mean(fullGuideVec)

    lfc_essential_mean_normalized = np.exp(lfc_essential-mean_of_fullGuideVec)#/mean_of_fullGuideVec
    lfc_nonessential_mean_normalized = np.exp(lfc_nonessential-mean_of_fullGuideVec)#/mean_of_fullGuideVec




    lfc_merged = np.concatenate((lfc_essential_mean_normalized, lfc_nonessential_mean_normalized), axis = 0)
    gamma = np.random.gamma(1,1,lfc_merged.size)[:,np.newaxis]



    # plt.figure()
    # seas.distplot(gamma)

    gammaNoise = np.random.gamma(1, 1, lfc_merged.size)[:,np.newaxis]
    final_counts_1 = np.random.poisson(depth_factor*gamma*lfc_merged * gammaNoise)
    gammaNoise = np.random.gamma(1, 1, lfc_merged.size)[:,np.newaxis]
    start_counts_1 = np.random.poisson(depth_factor*gamma * gammaNoise)


    # plt.figure()
    # seas.distplot(final_counts_1[:numEssentialGuides])
    # seas.distplot(final_counts_1[numEssentialGuides:])


    gammaNoise = np.random.gamma(1, 1, lfc_merged.size)[:,np.newaxis]
    final_counts_2 = np.random.poisson(depth_factor*gamma*lfc_merged * gammaNoise)
    gammaNoise = np.random.gamma(1, 1, lfc_merged.size)[:,np.newaxis]
    start_counts_2 = np.random.poisson(depth_factor*gamma*gammaNoise)



    # print(start_counts_1.shape)
    # print(start_counts_2.shape)
    # print(final_counts_1.shape)
    # print(final_counts_2.shape)
    # print(onesAndZeros.shape)
    # wait = input()

    finalData = np.concatenate((start_counts_1, start_counts_2, final_counts_1, final_counts_2, onesAndZeros), axis = 1)
    # np.savetxt("./input/finalSimData.csv", finalData, delimiter=",") #incase I wanted a NP version
    # creating new output file
    df_data= pd.DataFrame(data=finalData)
    df_data['guide_id'] = sg_Count_Names
    df_data['gene_id'] = gene_Names
    df_data = df_data[['gene_id', 'guide_id', 0, 1, 2, 3, 4]]
    df_data.columns = ['gene_id', 'guide_id', 'R1_T0', 'R2_T0', 'R1_F', 'R2_F', 'Essentiality']
    checkZero = df_data.copy().loc[:,['R1_T0', 'R2_T0', 'R1_F', 'R2_F']]
    checkZero['sum'] = checkZero.sum(axis=1)
    df_data = df_data[checkZero.loc[:,'sum']>0]

    with open('./input/Simulation.csv', "w") as myfile:
        df_data.to_csv(myfile, index=False, sep =',')

    call('/usr/local/bin/Rscript ./R_Scripts/SimToLFC_R.R', shell = True)



    # Creating input files for Controls

    checkZero = checkZero.iloc[-num_control_guides:,:]
    controlZeros = checkZero[checkZero.loc[:,'sum']<=0].shape[0]
    num_control_guides = num_control_guides - controlZeros

    control_l2fc = pd.read_csv('./input/Simulation_l2fc.csv', sep = ',')
    control_l2fc = control_l2fc.iloc[-num_control_guides:,1]
    control_guide_names = df_data.loc[:,'guide_id']
    control_guide_names = control_guide_names[-num_control_guides:]
    with open('./controls/Simulation_guide_names_control.csv', "w") as myfile:
        control_guide_names.to_csv(myfile, index=False, sep =',')
    with open('./controls/Simulation_l2fc_control.csv', "w") as myfile:
        control_l2fc.to_csv(myfile, index=False, header = 'x', sep =',')

    # CRISPhieRmix Specific File Creations
    CRISPhieRmix_Simulation = df_data.iloc[0:-num_control_guides,:]
    with open('./input/CRISPhieRmix_Simulation.csv', "w") as myfile:
        CRISPhieRmix_Simulation.to_csv(myfile, index=False, sep =',')
    CRISPhieRmix_Simulation_l2fc = pd.read_csv('./input/Simulation_l2fc.csv', sep = ',')
    CRISPhieRmix_Simulation_l2fc = CRISPhieRmix_Simulation_l2fc.iloc[0:-num_control_guides,:]
    with open('./input/CRISPhieRmix_Simulation_l2fc.csv', "w") as myfile:
        CRISPhieRmix_Simulation_l2fc.to_csv(myfile, sep = ',')

    #AnalyzeSim(rep_num)
    # plt.show()

    return()

def CombineDataIntoFiles(variables, total_reps):
    algorithms = ['RSA', 'RRA_MAGeCK', 'MLE_MAGeCK', 'CRISPhieRmix']
    for j in range(len(algorithms)):
        for i in range(total_reps):
            rep_num = i + 1
            if i == 0:
                alg_ROC_data = np.genfromtxt("./results/%s_Simulation_results_ROC_AUC_%s.txt"%(algorithms[j], str(rep_num)), dtype='str', skip_header=0)
                alg_ROC_data = np.array([alg_ROC_data])[:, np.newaxis]
                alg_PRC_data = np.genfromtxt("./results/%s_Simulation_results_PRC_AUC_%s.txt"%(algorithms[j], str(rep_num)), dtype='str', skip_header=0)
                alg_PRC_data = np.array([alg_PRC_data])[:, np.newaxis]
            else:
                temp_ROC_array = np.genfromtxt("./results/%s_Simulation_results_ROC_AUC_%s.txt"%(algorithms[j], str(rep_num)), dtype='str', skip_header=0)
                temp_ROC_array = temp_ROC_array[:,np.newaxis]
                alg_ROC_data = np.concatenate((alg_ROC_data,temp_ROC_array), axis = 1)

                temp_PRC_array = np.genfromtxt("./results/%s_Simulation_results_PRC_AUC_%s.txt"%(algorithms[j], str(rep_num)), dtype='str', skip_header=0)
                temp_PRC_array = temp_PRC_array[:,np.newaxis]
                alg_PRC_data = np.concatenate((alg_PRC_data,temp_PRC_array), axis = 1)

        alg_ROC_data = alg_ROC_data.astype(np.float)
        alg_PRC_data = alg_PRC_data.astype(np.float)


        x_axis = variables

        print(x_axis)
        print(x_axis.shape)
        print(alg_ROC_data)
        print(alg_ROC_data.shape)
        alg_ROC_data = np.concatenate((x_axis,alg_ROC_data), axis = 1)
        alg_PRC_data = np.concatenate((x_axis, alg_PRC_data), axis = 1)
        with open("./results/%s_ROC_Final.txt"%(algorithms[j]+'_results'), "w") as myfile:
            np.savetxt(myfile, alg_ROC_data, delimiter = '\t' )

        with open("./results/%s_PRC_Final.txt"%(algorithms[j]+'_results'), "w") as myfile:
            np.savetxt(myfile, alg_PRC_data, delimiter = '\t' )



    return()

##### hyperparameters ######
essential_genes_mean = 1.2 # 1.2
non_essential_variance = 0 # 0
essential_variance = 5 #5 seems to be similar to TKO
num_essential_genes = 60 # 100
num_nonessential_genes = 1700 # 1000
num_guides_per_gene = 5.0 # 5
guide_efficiency_mean = 0.8 # 0.8
guide_efficiency_variance = .1 # 0.1
depth_factor = 100 #sequencing depths starting at 100
total_reps = 1
num_control_guides = 500
binding_efficiency = 1 # 1
np.random.seed(0) #fixing for repeatability

variables = np.arange(1,100,5)
# variables = np.asarray([])
variables = np.array([100])[:,np.newaxis]


for i in range(total_reps):
    for j in range(variables.shape[0]):
        rep_num = i+1
        depth_factor = variables[j,0]
        print(depth_factor)
        np.random.seed(0)
        RunSimulation(essential_genes_mean, non_essential_variance, essential_variance, num_essential_genes,
                      num_nonessential_genes, num_guides_per_gene, guide_efficiency_mean, guide_efficiency_variance, num_control_guides, depth_factor, binding_efficiency, rep_num)

# plt.show()
# CombineDataIntoFiles(variables, total_reps)
# plt.show()


