import numpy as np
import pandas as pd
import pickle
import os.path
from subprocess import call
# import seaborn as sns
# import matplotlib.pyplot as plt
from sklearn import metrics
import sys


ff_name = "TKO_out.txt"
# variable = str(sys.argv[1])
# replication_number = str(sys.argv[2])
variable = 'TKO'
with open(ff_name, "w") as ff:
  print("opening TKO_out.txt", file = ff)

controlFile = 'NonHuman_TKO_Control_Guides.csv'

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
    call("python ./Algorithms/RSA/RSA.py -o ./output/%s_output.csv ./input/%s_Input.csv" %(filename,filename), shell=True)
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
    call('rscript ./R_Scripts/RunCRISPhieRmix.R %s %s %s' %(name, dataset, controlGuides), shell = True)
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

    with open("./results/%s_ROC_AUC_%s.txt"%(variable + '_'+ name+'_results', str(rep_num)), "a") as myfile:
        myfile.write(str(areaROC))
        myfile.write('\n')

    with open("./results/%s_PRC_AUC_%s.txt"%(variable + '_'+ name+'_results', str(rep_num)), "a") as myfile:
        myfile.write(str(pr_AUC))
        myfile.write('\n')


    fpr = fpr[:, np.newaxis]
    tpr = tpr[:, np.newaxis]
    roc_data = np.concatenate((fpr, tpr), axis = 1)

    with open("./results/%s_ROC_%s.npy"%(variable + '_'+ name+'_results', str(rep_num)), "wb") as myfile:
        np.save(file=myfile, arr=roc_data)

    precision = precision[:, np.newaxis]
    recall = recall[:, np.newaxis]
    pr_data = np.concatenate((recall, precision), axis = 1)

    with open("./results/%s_PR_%s.npy"%(variable + '_'+ name+'_results', str(rep_num)), "wb") as myfile:
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



def CheckFDR(output, fdrs, essential, totalGenes, name, rep_num = 1):
    totalEssentialGenes = len(essential)
    totalRows = output.shape[0]
    essentialCounted = 0
    genesSeen = set()
    essentialCountedSet = set()
    print('total essential gene count = %d' %(totalEssentialGenes))

    FDR_actual = np.zeros(totalGenes)
    FDR_predicted = np.zeros(totalGenes)


    for i in range(totalRows):
        gene = output.iloc[i,:].loc['Gene_ID']
        if gene not in genesSeen:
            genesSeen.add(gene)
            if gene in essential:
                essentialCounted = essentialCounted + 1
                essentialCountedSet.add(gene)
            FDR_actual[len(genesSeen)-1] = 1-float(essentialCounted)/len(genesSeen)
            FDR_predicted[len(genesSeen)-1] = fdrs[len(genesSeen)-1]

    FDR_predicted = FDR_predicted[FDR_predicted<.3]
    FDR_actual = FDR_actual[0:len(FDR_predicted)]

    sns.lineplot(FDR_predicted, FDR_actual, label='%s_FDR'%(name))
    plt.ylim((0,1))
    plt.xlim((0,1))

    FDR_actual = FDR_actual[:,np.newaxis]
    FDR_predicted = FDR_predicted[:,np.newaxis]
    areaFDR = np.concatenate((FDR_actual,FDR_predicted), axis = 1)
    areaFDR = areaFDR[areaFDR[:,0].argsort()]
    areaFDR = round(metrics.auc(areaFDR[:,0], areaFDR[:,1]),5)
    with open("./results/%s_FDR_AUC_%s.txt"%(name+'_results', str(rep_num)), "a") as myfile:
        myfile.write(str(areaFDR))
        myfile.write('\n')

    FDR_actual = FDR_actual[:,np.newaxis]
    FDR_predicted = FDR_predicted[:, np.newaxis]
    fdr_data = np.concatenate((FDR_predicted, FDR_actual), axis = 1)
    with open("./results/%s_FDR_%s.npy"%(name+'_results', str(rep_num)), "wb") as myfile:
        np.save(file=myfile, arr=fdr_data)

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
    # if algo == 'RSA' or algo == 'All':
    # Run RSA Analysis
    #     RunRSA('RSA_%s'%name, isSim, rep_num)
    #
    #
    # if algo == 'RRA' or algo == 'All':
        # Run MAGeCK RRA Analysis
        # RunMAGeCKRRA('MAGeCK_%s'%name, isSim, rep_num)

    if algo == 'MLE' or algo == 'All':
        # Run MAGeCK MLE Analysis
        RunMAGeCKMLE('MAGeCK_%s'%name, isSim, rep_num)
    #
    # if algo == "CRISPhieRmix" or algo == 'All':
    #     RunCRISPhieRmix(name, isSim, rep_num)


    # Running FDR Analysis here:
    # fig = plt.figure()
    # fig.suptitle('FDR: predicted vs actual')
    # if algo == 'RRA' or algo == 'All':
        # Run MAGeCK RRA Analysis
        # RunFDRARRA('MAGeCK_%s'%name, isSim, rep_num)

    if algo == 'MLE' or algo == 'All':
        # Run MAGeCK MLE Analysis
        RunFDRMLE('MAGeCK_%s'%name, isSim, rep_num)

    # if algo == "CRISPhieRmix" or algo == 'All':
    #     RunFDRCRISPhieRmix(name, isSim, rep_num)

    # plt.show()

def AnalyzeSim(rep_num = 1):
    RunAlgorithms('Simulation', 'All', rep_num)

def main():


    if not os.path.isfile('./pickled/essential_nonessential_genes.p'):
        rd.ReadEssential_NonEssential_Genes()
        print('Finished Reading Harts Essentiality Genes')
    essential_genes, nonessential_genes = pickle.load(open( "./pickled/essential_nonessential_genes.p", "rb" ))

    if not os.path.isfile('./input/GBM-lib1_l2fc.csv'):
        rd.ProcessTKODataForDESEQ()



    with open(ff_name, "a") as ff:
        print('#################################', file = ff)

    filenames = ['HCT116_1-lib1', 'HCT116_2-lib1', 'DLD1-lib1', 'HeLa-lib1', 'RPE1-lib1','GBM-lib1']


    for i in range(len(filenames)):
        with open(ff_name, "a") as ff:
            print(filenames[i], file = ff)
        RunAlgorithms(filenames[i], 'MLE')
    return()


main()
# plt.show()
# AnalyzeSim()

