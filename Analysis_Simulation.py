import Read_Data as rd
import numpy as np
import pandas as pd
import pickle
import os.path
from subprocess import call
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn import metrics
from scipy import stats
import statsmodels.stats.multitest as smt

import sys


ff_name = "TKO_out.txt"
# variable = str(sys.argv[1])
# Sherlock_filename = str(sys.argv[2])

with open(ff_name, "w") as ff:
  print("opening TKO_out.txt", file = ff)

controlFile = 'NonHuman_TKO_Control_Guides.csv'
fdr_cutoff = 0.1

def ConvertTKOForRSA(filename):
    print('converting %s' % (filename))
    countData = pd.read_csv('./input/%s.csv' % (filename), sep=',', header=0)
    lfc = pd.read_csv('./input/%s_l2fc.csv' % (filename), sep=',', header=0)
    inputRSA = countData[['GENE', 'GENE_CLONE']]
    inputRSA['Score'] = lfc['x']
    inputRSA.columns = ['Gene_ID', 'Well_ID', 'Score']
    with open('./input/RSA_%s_Input.csv' % (filename), "w") as myfile:
        inputRSA.to_csv(myfile, index=False, sep=',')
    return inputRSA

def ConvertSimForRSA():
    simDatadf = pd.read_csv('./input/Simulation.csv', sep=',', header=0)
    simlfc = pd.read_csv('./input/Simulation_l2fc.csv', sep=',',header=0)
    inputRSA = simDatadf[['gene_id', 'guide_id']]
    inputRSA['Score'] = simlfc['x']
    inputRSA.columns = ['Gene_ID', 'Well_ID', 'Score']
    with open('./input/RSA_Simulation_Input.csv', "w") as myfile:
        inputRSA.to_csv(myfile, index=False, sep=',')

    return inputRSA

def RunRSA(filename, isSim, rep_num = 1):
    # https://admin-ext.gnf.org/publications/RSA/
    print('============RUNNING RSA===========')
    call("python ./Algorithms/RSA/RSA.py -o ./output/%s_output.csv ./input/%s_Input.csv" %(filename,filename), shell=True)#TKT: Not shared in the git
    output = pd.read_csv('./output/%s_output.csv' %(filename))

    if isSim:
        rawData = pd.read_csv('./input/Simulation.csv', header=0, sep=',')
        totalGenes = len(set(rawData['gene_id']))
        essentialGenes = rawData[rawData['Essentiality'] == 1.0]
        essentialGenes = set(essentialGenes['gene_id'])
    else:
        essential_genes = pd.read_csv('./Hart/Essential_Genes.txt', sep='\t', header=0).drop(columns=['HGNC_ID'])
        genesInData= set(output.loc[:,'Gene_ID'])
        totalGenes=len(genesInData)
        essGenes = set(essential_genes.loc[:,'Gene'])
        essentialGenes = genesInData.intersection(essGenes)

    AnalyzeAlgorithmOutput(output, essentialGenes, totalGenes, filename, rep_num)
    return

def ConvertSimForMAGeCK():
    simDatadf = pd.read_csv('./input/Simulation.csv', sep=',', header=0)
    inputMAGeCK = simDatadf[['guide_id', 'gene_id','R1_T0', 'R2_T0', 'R1_F', 'R2_F']]
    inputMAGeCK.columns = ['sgRNA', 'gene_id', 'R1_T0', 'R2_T0', 'R1_F', 'R2_F']
    with open('./input/MAGeCK_Simulation_Input.csv', "w") as myfile:
        inputMAGeCK.to_csv(myfile, index=False, sep=',')
    return inputMAGeCK

def ConvertTKOForMAGeCK(filename):
    data_df = pd.read_csv('./input/%s.csv'%filename, sep=',', header=0)

    if data_df.shape[1] == 5:
        inputMAGeCK = data_df.loc[:,['GENE_CLONE', 'GENE', 'T0', 'R1_F', 'R2_F']]
    else:
        inputMAGeCK = data_df.loc[:,['GENE_CLONE', 'GENE', 'T0', 'R1_F', 'R2_F', 'R3_F']]

    inputMAGeCK.rename(columns={'GENE':'gene_id', 'GENE_CLONE':'sgRNA'}, inplace=True)
    with open('./input/MAGeCK_%s_Input.csv' %filename, "w") as myfile:
        inputMAGeCK.to_csv(myfile, index=False, sep=',')
    return

def RunMAGeCKRRA(name, isSim, rep_num = 1):
    print('============RUNNING RRA===========')
    filename = name + str('_Input.csv')
    inputFile = pd.read_csv('./input/%s'%filename, sep=',', header=0)
    # Call the file to build the ranked list
    if isSim is False:
        if inputFile.shape[1] == 5:
            call('mageck test -k ./input/%s -t R1_F,R2_F -c T0 --control-sgrna ./controls/%s -n %s_RRA_Output' %(filename, controlFile, name), shell = True)
        else:
            call('mageck test -k ./input/%s -t R1_F,R2_F,R3_F -c T0 --control-sgrna ./controls/%s -n %s_RRA_Output' %(filename, controlFile, name), shell = True)
    else:
        call('mageck test -k ./input/%s -t R1_F,R2_F -c R1_T0,R2_T0  --control-sgrna ./controls/Simulation_guide_names_control.csv -n %s_RRA_Output' %(filename,name), shell = True) #--control-sgrna ./input/Simulation_guide_names_control.csv

    # Calling commands to move the file to the appropriate location and remove extra files
    call('mv -f %s_RRA_Output.gene_summary.txt ./output/' %(name), shell=True)
    call('rm %s_RRA_Output.log' %(name), shell=True)
    call('rm %s_RRA_Output.R' %(name), shell=True)
    call('rm %s_RRA_Output_summary.Rnw' %(name), shell=True)
    call('rm %s_RRA_Output.sgrna_summary.txt' %(name), shell=True)
    if isSim:
        rawData = pd.read_csv('./input/Simulation.csv', header=0, sep=',')
        totalGenes = len(set(rawData['gene_id']))
        essentialGenes = rawData[rawData['Essentiality'] == 1.0]
        essentialGenes = set(essentialGenes['gene_id'])
    else:
        essential_genes = pd.read_csv('./Hart/Essential_Genes.txt', sep='\t', header=0).drop(columns=['HGNC_ID'])
        genesInData= set(inputFile.loc[:,'gene_id'])
        totalGenes=len(genesInData)
        essGenes = set(essential_genes.loc[:,'Gene'])
        essentialGenes = genesInData.intersection(essGenes)


    output = pd.read_csv('./output/%s_RRA_Output.gene_summary.txt' %(name), sep='\t')
    output['Gene_ID'] = output['id']
    AnalyzeAlgorithmOutput(output, essentialGenes, totalGenes, 'RRA_'+name, rep_num)
    return

def RunMAGeCKMLE(name, isSim, rep_num = 1):
    print('============RUNNING MLE===========')
    filename = name + str('_Input.csv')
    inputFile = pd.read_csv('./input/%s'%filename, sep=',', header=0)
    if isSim is False:
        if inputFile.shape[1] == 5 and isSim is False:
        # Call the file to build the ranked list
            call('mageck mle -k ./input/%s -d ./DesignMats/MAGeCK_2Rep_DesignMat.txt --control-sgrna ./controls/%s  -n %s_MLE_Output' %(filename, controlFile, name), shell = True)
        else:
            call('mageck mle -k ./input/%s -d ./DesignMats/MAGeCK_3Rep_DesignMat.txt --control-sgrna ./controls/%s -n %s_MLE_Output' %(filename, controlFile, name), shell = True)
    else:
        call('mageck mle -k ./input/%s -d ./DesignMats/%s_DesignMat.txt --control-sgrna ./controls/Simulation_guide_names_control.csv -n %s_MLE_Output' %(filename,name, name), shell = True)
    # # Calling commands to move the file to the appropriate location and remove extra files
    call('mv -f %s_MLE_Output.gene_summary.txt ./output/' %(name), shell=True)
    call('rm %s_MLE_Output.log' %(name), shell=True)
    call('rm %s_MLE_Output.sgrna_summary.txt' %(name), shell=True)
    if isSim:
        rawData = pd.read_csv('./input/Simulation.csv', header=0, sep=',')
        totalGenes = len(set(rawData['gene_id']))
        essentialGenes = rawData[rawData['Essentiality'] == 1.0]
        essentialGenes = set(essentialGenes['gene_id'])
    else:
        essential_genes = pd.read_csv('./Hart/Essential_Genes.txt', sep='\t', header=0).drop(columns=['HGNC_ID'])
        genesInData= set(inputFile.loc[:,'gene_id'])
        totalGenes=len(genesInData)
        essGenes = set(essential_genes.loc[:,'Gene'])
        essentialGenes = genesInData.intersection(essGenes)

    output = pd.read_csv('./output/%s_MLE_Output.gene_summary.txt' %(name), sep='\t')
    output['Gene_ID'] = output['Gene']
    output.sort_values(by='Sample|fdr', ascending=True, inplace=True, axis=0)
    AnalyzeAlgorithmOutput(output, essentialGenes, totalGenes, 'MLE_'+ name, rep_num)
    return

def ConvertTKOForCRISPhieRmix(name):
    # needed to get TKO control guides
    l2f = pd.read_csv('./input/%s_l2fc.csv'%name, sep=',')
    origData = pd.read_csv('./input/%s.csv'%name, sep=',')
    controlGuides = pd.read_csv('./controls/%s' %controlFile, header=None)
    controlGuides = set(controlGuides.iloc[:,0])

    # print(controlGuides)
    control2fc = l2f[origData.loc[:,'GENE_CLONE'].isin(controlGuides)]
    with open('./controls/CRISPhieRmix_TKO_CG.csv', "w") as myfile:
        control2fc.to_csv(myfile, sep=',')

    x = l2f[~origData.loc[:,'GENE_CLONE'].isin(controlGuides)]
    with open('./input/CRISPhieRmix_%s_l2fc.csv'%name, "w") as myfile:
        x.to_csv(myfile, sep=',')
    guideNames = origData[~origData.loc[:,'GENE_CLONE'].isin(controlGuides)]
    guideNames.rename(columns={'GENE':'gene_id', 'GENE_CLONE':'sgRNA'}, inplace=True)
    guideNames.to_csv('./input/CRISPhieRmix_%s.csv'%name, sep=',')

    return()

def RunCRISPhieRmix(name, isSim, rep_num = 1):

    print('============RUNNING CRISPhieRmix===========')
    dataset = None
    if name != 'Simulation':
        dataset = 'TKO'
        controlGuides = 'CRISPhieRmix_TKO_CG.csv'

    else:
        controlGuides = 'Simulation_l2fc_control.csv'
        dataset = 'Simulation'

    print('Checkpoint0')
    call('/usr/local/bin/Rscript ./R_Scripts/RunCRISPhieRmix.R %s %s %s' %(name, dataset, controlGuides), shell = True)
    output = pd.read_csv('./output/CRISPhieRmix_%s_Output.csv' %(name))
    output = output.loc[:,['genes', 'ranks']]
    print(output.head())
    output.sort_values(by=['ranks'], inplace=True)
    print(output.head())
    output['Gene_ID'] = output['genes']
    if isSim:
        rawData = pd.read_csv('./input/Simulation.csv', header=0, sep=',')
        totalGenes = len(set(rawData['gene_id']))
        essentialGenes = rawData[rawData['Essentiality'] == 1.0]
        essentialGenes = set(essentialGenes['gene_id'])
    else:
        essential_genes = pd.read_csv('./Hart/Essential_Genes.txt', sep='\t', header=0).drop(columns=['HGNC_ID'])
        genesInData= set(output.loc[:,'Gene_ID'])
        totalGenes=len(genesInData)
        essGenes = set(essential_genes.loc[:,'Gene'])
        essentialGenes = genesInData.intersection(essGenes)

    AnalyzeAlgorithmOutput(output, essentialGenes, totalGenes, 'CRISPhieRmix_'+name, rep_num)
    return

def RunTTest(name, isSim, rep_num = 1):
    if isSim:
        SimOriginal = pd.read_csv('./input/Simulation.csv')
        SimL2fc = pd.read_csv('./input/Simulation_l2fc.csv')
        controlGuideL2fc = pd.read_csv('./controls/Simulation_l2fc_control.csv').values
        controlGuideGeneNames = pd.read_csv('./controls/Simulation_guide_names_control.csv', header=None)
        targetingGenes = SimOriginal[~SimOriginal.loc[:,'guide_id'].isin(controlGuideGeneNames.loc[:,0])]
        targetingGenes = set(targetingGenes.loc[:,'gene_id'])
        targetingGenes = list(targetingGenes)


        fdrs = pd.DataFrame(columns=['Gene_ID', 'fdrs'])

        geneList = []
        fdrList = []

        for i in range(len(targetingGenes)):
            tempGene = targetingGenes[i]
            geneGuides = SimL2fc[SimOriginal.loc[:,'gene_id'].isin([tempGene])].loc[:, 'x'].values[:,np.newaxis]
            _ , tempProb = stats.ttest_ind(geneGuides, controlGuideL2fc)
            geneList.append(tempGene)
            fdrList.append(tempProb[0])

        fdrs['Gene_ID'] = geneList
        fdrs['fdrs'] = fdrList
        fdrs.sort_values(by='fdrs', ascending=True, inplace=True, axis=0)
        rawData = pd.read_csv('./input/Simulation.csv', header=0, sep=',')
        totalGenes = len(set(rawData['gene_id']))
        essentialGenes = rawData[rawData['Essentiality'] == 1.0]
        essentialGenes = set(essentialGenes['gene_id'])
        fdrs.to_csv(path_or_buf='./output/TTest_%s_Output.csv' %(name), sep=',', header=True, index=False)

    else:
        print('T test real here')
        # get a list of the control guides
        controlGuideL2fc = pd.read_csv('./controls/CRISPhieRmix_TKO_CG.csv', header=0, sep = ',')
        # get a list of the genes as well as their respective log fold changes
        genel2fc = pd.read_csv('./input/CRISPhieRmix_%s_l2fc.csv'%name, header=0, sep=',')
        geneIds_df = pd.read_csv('./input/CRISPhieRmix_%s.csv'%name, header=0, sep=',')
        targetingGenes = list(set(geneIds_df.loc[:,'gene_id']))


        fdrs = pd.DataFrame(columns=['Gene_ID', 'fdrs'])

        geneList = []
        fdrList = []

        for i in range(len(targetingGenes)):
            tempGene = targetingGenes[i]
            geneGuides = genel2fc[geneIds_df.loc[:,'gene_id'].isin([tempGene])].loc[:, 'x'].values[:,np.newaxis]
            _ , tempProb = stats.ttest_ind(geneGuides, controlGuideL2fc)
            geneList.append(tempGene)
            fdrList.append(tempProb[0])

        fdrs['Gene_ID'] = geneList
        fdrs['fdrs'] = fdrList
        fdrs.sort_values(by='fdrs', ascending=True, inplace=True, axis=0)

        essential_genes = pd.read_csv('./Hart/Essential_Genes.txt', sep='\t', header=0).drop(columns=['HGNC_ID'])
        totalGenes=len(targetingGenes)
        essGenes = set(essential_genes.loc[:,'Gene'])
        essentialGenes = set(targetingGenes).intersection(essGenes)


    AnalyzeAlgorithmOutput(fdrs, essentialGenes, totalGenes, 'TTest_'+name, rep_num)

    return


# Takes a ranked list of genes. Takes a set of essential genes. Takes the total number of overall genes, takes the name of the algorithm
def AnalyzeAlgorithmOutput(output, essential, totalGenes, name, rep_num = 1):
    totalEssentialGenes = len(essential)
    totalRows = output.shape[0]
    essentialCounted = 0
    genesSeen = set()
    essentialCountedSet = set()
    print('total essential gene count = %d' %(totalEssentialGenes))
    EssentialROCArray = np.zeros(totalGenes)
    TPR = np.zeros(totalGenes)
    FDR = np.zeros(totalGenes)
    trueLabel = np.zeros(totalGenes)
    confidence = np.linspace(1,0,totalGenes)

    for i in range(totalRows):
        gene = output.iloc[i,:].loc['Gene_ID']

        if gene not in genesSeen:
            genesSeen.add(gene)
            if gene in essential:
                essentialCounted = essentialCounted + 1
                essentialCountedSet.add(gene)
                trueLabel[len(genesSeen)-1] = 1
            else:
                trueLabel[len(genesSeen)-1] = 0
        # EssentialROCArray[len(genesSeen)-1] = essentialCounted/totalEssentialGenes
        # TPR[len(genesSeen)-1] = essentialCounted*1.0/totalEssentialGenes
        # nonEssentialCounted = len(genesSeen)-essentialCounted
        # FDR[len(genesSeen)-1] = (nonEssentialCounted)*1.0/(totalGenes-totalEssentialGenes)

    fpr, tpr, _ = metrics.roc_curve(trueLabel, confidence)
    print(trueLabel.shape)
    print(confidence.shape)

    areaROC = round(metrics.roc_auc_score(trueLabel, confidence),5)
    sns.lineplot(fpr, tpr, label = '%s_ROC:%s'%(name, str(areaROC)))
    precision, recall, _ = metrics.precision_recall_curve(trueLabel, confidence)
    pr_AUC = round(metrics.auc(recall, precision), 5)

    # sns.lineplot(recall, precision, label = '%s_Precision_Recall:%s'%(name, str(pr_AUC)))

    with open("./results/%s_ROC_AUC_%s.txt"%(name+'_results', str(rep_num)), "a") as myfile:
        myfile.write(str(areaROC))
        myfile.write('\n')

    with open("./results/%s_PRC_AUC_%s.txt"%(name+'_results', str(rep_num)), "a") as myfile:
        myfile.write(str(pr_AUC))
        myfile.write('\n')


    fpr = fpr[:, np.newaxis]
    tpr = tpr[:, np.newaxis]
    roc_data = np.concatenate((fpr, tpr), axis = 1)

    with open("./results/%s_ROC_%s.npy"%(name+'_results', str(rep_num)), "wb") as myfile:
        np.save(file=myfile, arr=roc_data)

    precision = precision[:, np.newaxis]
    recall = recall[:, np.newaxis]
    pr_data = np.concatenate((recall, precision), axis = 1)

    with open("./results/%s_PR_%s.npy"%(name+'_results', str(rep_num)), "wb") as myfile:
        np.save(file=myfile, arr=pr_data)

    return

def RunFDRARRA(name, isSim, rep_num = 1):
    # FDR for RRA:
    print('============RUNNING RRA FDR===========')
    filename = name + str('_Input.csv')
    inputFile = pd.read_csv('./input/%s'%filename, sep=',', header=0)
    if isSim:
        rawData = pd.read_csv('./input/Simulation.csv', header=0, sep=',')
        totalGenes = len(set(rawData['gene_id']))
        essentialGenes = rawData[rawData['Essentiality'] == 1.0]
        essentialGenes = set(essentialGenes['gene_id'])
    else:
        essential_genes = pd.read_csv('./Hart/Essential_Genes.txt', sep='\t', header=0).drop(columns=['HGNC_ID'])
        genesInData= set(inputFile.loc[:,'gene_id'])
        totalGenes=len(genesInData)
        essGenes = set(essential_genes.loc[:,'Gene'])
        essentialGenes = genesInData.intersection(essGenes)


    output = pd.read_csv('./output/%s_RRA_Output.gene_summary.txt' %(name), sep='\t')
    output['Gene_ID'] = output['id']

    # Checking FDR rate accuracy
    CheckFDR(output, output.loc[:,'neg|fdr'].values, essentialGenes, totalGenes, 'RRA_' + name, rep_num)
    return()

def RunFDRMLE(name, isSim, rep_num=1):
    print('============RUNNING MLE FDR===========')
    filename = name + str('_Input.csv')
    inputFile = pd.read_csv('./input/%s'%filename, sep=',', header=0)

    if isSim:
        rawData = pd.read_csv('./input/Simulation.csv', header=0, sep=',')
        totalGenes = len(set(rawData['gene_id']))
        essentialGenes = rawData[rawData['Essentiality'] == 1.0]
        essentialGenes = set(essentialGenes['gene_id'])
    else:
        essential_genes = pd.read_csv('./Hart/Essential_Genes.txt', sep='\t', header=0).drop(columns=['HGNC_ID'])
        genesInData= set(inputFile.loc[:,'gene_id'])
        totalGenes=len(genesInData)
        essGenes = set(essential_genes.loc[:,'Gene'])
        essentialGenes = genesInData.intersection(essGenes)

    output = pd.read_csv('./output/%s_MLE_Output.gene_summary.txt' %(name), sep='\t')
    output['Gene_ID'] = output['Gene']
    output.sort_values(by='Sample|fdr', ascending=True, inplace=True, axis=0)
    CheckFDR(output, output.loc[:,'Sample|fdr'].values, essentialGenes, totalGenes, 'MLE_' + name, rep_num)

    return()

def RunFDRCRISPhieRmix(name, isSim, rep_num = 1):
    print('============RUNNING CRISPhieRmix FDR ===========')
    dataset = None

    output = pd.read_csv('./output/CRISPhieRmix_%s_Output.csv' %(name))
    output = output.loc[:,['genes', 'ranks']]
    output.sort_values(by=['ranks'], inplace=True)
    output['Gene_ID'] = output['genes']
    if isSim:
        rawData = pd.read_csv('./input/Simulation.csv', header=0, sep=',')
        totalGenes = len(set(rawData['gene_id']))
        essentialGenes = rawData[rawData['Essentiality'] == 1.0]
        essentialGenes = set(essentialGenes['gene_id'])
    else:
        essential_genes = pd.read_csv('./Hart/Essential_Genes.txt', sep='\t', header=0).drop(columns=['HGNC_ID'])
        genesInData= set(output.loc[:,'Gene_ID'])
        totalGenes=len(genesInData)
        essGenes = set(essential_genes.loc[:,'Gene'])
        essentialGenes = genesInData.intersection(essGenes)


    CheckFDR(output, output.loc[:,'ranks'].values, essentialGenes, totalGenes, 'CRISPhieRmix_' + name, rep_num)
    return()

def RunFDRTTest(name, isSim, rep_num = 1):
    print('============RUNNING TTest FDR ===========')
    dataset = None

    output = pd.read_csv('./output/TTest_%s_Output.csv' %(name))

    if isSim:
        rawData = pd.read_csv('./input/Simulation.csv', header=0, sep=',')
        totalGenes = len(set(rawData['gene_id']))
        essentialGenes = rawData[rawData['Essentiality'] == 1.0]
        essentialGenes = set(essentialGenes['gene_id'])
        _, BHCorrected, _, _ = smt.multipletests(pvals=output.loc[:,'fdrs'], alpha=fdr_cutoff, method='fdr_bh')
        output['BHCorrected'] = BHCorrected
        CheckFDR(output, output.loc[:,'fdrs'].values, essentialGenes, totalGenes, 'TTest_' + name, rep_num)

    else:
        return()




    return()

def CheckFDR(output, total_fdrs, essential, totalGenesNum, name, rep_num = 1):

    totalEssentialGenes = len(essential)
    totalRows = output.shape[0]
    essentialCounted = 0
    genesSeen = set()
    essentialCountedSet = set()
    print('total essential gene count = %d' %(totalEssentialGenes))

    FDR_actual = np.zeros(totalGenesNum)
    FDR_predicted = np.zeros(totalGenesNum)


    # going through the cutoffs



    cutoffs = [0.05, 0.1, 0.15, 0.2]

    for i in range(len(cutoffs)):
        if 'TTest' in name:
            _, total_fdrs_bh, _, _ = smt.multipletests(pvals=total_fdrs, alpha=cutoffs[i], method='fdr_bh')
            fdrs = total_fdrs_bh[:-600]
        else:
            fdrs = total_fdrs[:-600]

        total_genes = output.loc[:,'Gene_ID']
        genes = total_genes[:-600]


        genes = genes[fdrs<cutoffs[i]]
        essentialsInGenes = genes[genes.isin(essential)]
        if genes.shape[0] == 0:
            fdr = 0
            precision = 0
        else:
            fdr = 1 - essentialsInGenes.shape[0] / genes.shape[0]
            precision = essentialsInGenes.shape[0] / genes.shape[0]
        recall = essentialsInGenes.shape[0] / totalEssentialGenes
        if precision == 0.0 and recall == 0.0:
            F1 = 0
        else:
            F1 = 2*precision*recall/(precision + recall)


        with open("./results/%s_FDR_F1_%s_Cutoff_%s.txt"%(name+'_results', str(rep_num), str(cutoffs[i])), "a") as myfile:
            myfile.write(str(F1))
            myfile.write('\n')
        with open("./results/%s_FDR_%s_Cutoff_%s.txt"%(name+'_results', str(rep_num), str(cutoffs[i])), "a") as myfile:
            myfile.write(str(fdr))
            myfile.write('\n')


    with open("./results/%s_FDR_%s.npy"%(name+'_results', str(rep_num)), "wb") as myfile:
        np.save(file=myfile, arr=(total_fdrs))

    return()

def RunAlgorithms(name, algo, rep_num=1):
    # fig = plt.figure()
    # fig.suptitle('ROC curve')

    if name == 'Simulation':
        isSim = True
    else:
        isSim = False

    # Creating appropriate input files
    if isSim:
        ConvertSimForRSA()
        ConvertSimForMAGeCK()
    else:
        ConvertTKOForRSA(name)
        ConvertTKOForMAGeCK(name)
        ConvertTKOForCRISPhieRmix(name)

    # Running ROC/ AR Curve analysis
    #if algo == 'RSA' or algo == 'All': #TKT: Algo not shared
    # Run RSA Analysis
    #    RunRSA('RSA_%s'%name, isSim, rep_num)
    #
    #
    if algo == 'RRA' or algo == 'All':
    #     Run MAGeCK RRA Analysis
        RunMAGeCKRRA('MAGeCK_%s'%name, isSim, rep_num)
    #
    if algo == 'MLE' or algo == 'All':
    #     Run MAGeCK MLE Analysis
        RunMAGeCKMLE('MAGeCK_%s'%name, isSim, rep_num)

    if algo == "CRISPhieRmix" or algo == 'All':
        RunCRISPhieRmix(name, isSim, rep_num)

    if algo == "TTest" or algo == "All":
        RunTTest(name, isSim, rep_num)

    # Running FDR Analysis here:
    # fig = plt.figure()
    # fig.suptitle('FDR: predicted vs actual')
    if algo == 'RRA' or algo == 'All':
        # Run MAGeCK RRA Analysis
        RunFDRARRA('MAGeCK_%s'%name, isSim, rep_num)
    #
    if algo == 'MLE' or algo == 'All':
        # Run MAGeCK MLE Analysis
        RunFDRMLE('MAGeCK_%s'%name, isSim, rep_num)

    if algo == "CRISPhieRmix" or algo == 'All':
        RunFDRCRISPhieRmix(name, isSim, rep_num)

    if algo == "TTest" or algo == 'All':
        RunFDRTTest(name, isSim, rep_num)

    # plt.show()

def AnalyzeSim(rep_num = 1):
    RunAlgorithms('Simulation', 'All', rep_num)

def main():
    # Reading Harts Essential Genes
    if not os.path.isfile('./pickled/essential_nonessential_genes.p'):
        rd.ReadEssential_NonEssential_Genes()
        print('Finished Reading Harts Essentiality Genes')
    essential_genes, nonessential_genes = pickle.load(open( "./pickled/essential_nonessential_genes.p", "rb" ))

    if not os.path.isfile('./input/GBM-lib1_l2fc.csv'):
        rd.ProcessTKODataForDESEQ()

    plt.show()	
    RunAlgorithms('Simulation', 'All') #TKT

    # rd.ProcessTKODataForDESEQ()
    #with open(ff_name, "a") as ff:
    #    print('#################################', file = ff)

    #filenames = ['HCT116_1-lib1', 'HCT116_2-lib1', 'DLD1-lib1', 'HeLa-lib1', 'RPE1-lib1','GBM-lib1']
    # RunAlgorithms(filenames[1], 'All')


    #for i in range(len(filenames)):
    #    with open(ff_name, "a") as ff:
    #        print(filenames[i], file = ff)
    #    RunAlgorithms(filenames[i], 'All')
    #return()

main()
# plt.show()
# AnalyzeSim()
# RunAlgorithms('Simulation', 'CRISPhieRmix', 1)

