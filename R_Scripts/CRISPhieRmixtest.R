f_name = 'Simulation'
dataset = 'Simulation'
controlGuideName = 'Simulation_l2fc_control.csv'

l2f_filename = paste("/Users/sbodapati/Documents/Qi Lab/CRISPR_Benchmarking/input/CRISPhieRmix_",f_name,"_l2fc",".csv", sep='')
l2fc = read.table(file = l2f_filename, header = TRUE, sep=',')
log2fc = l2fc$x
print(length(log2fc))

filename = paste("/Users/sbodapati/Documents/Qi Lab/CRISPR_Benchmarking/input/CRISPhieRmix_",f_name,".csv", sep='')
testCounts = read.table(file = filename, header = TRUE, sep=',')
geneID = testCounts$gene_id
geneID = factor(geneID,levels=unique(geneID))
print(length(geneID))




controlGuides_filename = paste("/Users/sbodapati/Documents/Qi Lab/CRISPR_Benchmarking/controls/",controlGuideName, sep='')
controlGuides = read.table(file = controlGuides_filename, header = TRUE, sep=',')
control2fc = controlGuides$x
print(length(control2fc))
print("________________here_New______________")
#install.packages('mixtools')
#install.packages('devtools')
#devtools::install_github("timydaley/CRISPhieRmix")
print("________________here_New1______________")
print(length(log2fc))
print(length(control2fc))
print(length(geneID))

#if (f_name == "Simulation"){
if (length(control2fc)>0){
  CRISPhieRmixResults = CRISPhieRmix::CRISPhieRmix(x = log2fc, geneIds = geneID, negCtrl = control2fc, BIMODAL=TRUE) # 
} else {
  CRISPhieRmixResults = CRISPhieRmix::CRISPhieRmix(x = log2fc, geneIds = geneID, BIMODAL=TRUE) # 
}
#}

#if (dataset == "TKO"){
# CRISPhieRmixResults = CRISPhieRmix::CRISPhieRmix(x = log2fc, geneIds = geneID, negCtrl = control2fc)
#  print("________________TKO______________")

#}

geneRanks = CRISPhieRmixResults$FDR 
genes = CRISPhieRmixResults$genes 

print(length(genes))
print(length(geneRanks))

final_data <- data.frame(genes = genes , ranks = geneRanks)
outputFilename = paste('/Users/sbodapati/Documents/Qi Lab/CRISPR_Benchmarking/output/','CRISPhieRmix_',f_name,'_Output_Test.csv', sep='')
write.csv(final_data, outputFilename)

