args <- commandArgs()

f_name = toString('Simulation')
dataset = toString('Simulation')
controlGuideName = toString('Simulation_l2fc_control.csv')

l2f_filename = paste("/Users/sbodapati/Documents/Qi Lab/CRISPR_Benchmarking/input/CRISPhieRmix_",f_name,"_l2fc",".csv", sep='')
l2fc = read.table(file = l2f_filename, header = TRUE, sep=',')
log2fc = l2fc$x
print(length(log2fc))

filename = paste("/Users/sbodapati/Documents/Qi Lab/CRISPR_Benchmarking/input/CRISPhieRmix_",f_name,".csv", sep='')
testCounts = read.table(file = filename, header = TRUE, sep=',')
geneID = testCounts$gene_id[1:length(log2fc)]
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
print(head(log2fc))
print(head(control2fc))
print(head(geneID))

#if (f_name == "Simulation"){
if (length(control2fc)>0){
  CRISPhieRmixResults = CRISPhieRmix::CRISPhieRmix(x = log2fc, geneIds = geneID, negCtrl = control2fc, BIMODAL=TRUE) # 
  geneRanks = CRISPhieRmixResults$negFDR
} else {
  CRISPhieRmixResults = CRISPhieRmix::CRISPhieRmix(x = log2fc, geneIds = geneID, BIMODAL=TRUE) # 
  geneRanks = CRISPhieRmixResults$negFDR
}

if(!exists("CRISPhieRmixResults")){
  # failure CRISPhieRmix, run in unimodal mode
  print('we get into this section!!!!!!!!!!!!!!!!!!!!')
  if (length(control2fc)>0){
    CRISPhieRmixResults = CRISPhieRmix::CRISPhieRmix(x = log2fc, geneIds = geneID, negCtrl = control2fc, mu = -4) # 
    geneRanks = CRISPhieRmixResults$FDR
  } else {
    CRISPhieRmixResults = CRISPhieRmix::CRISPhieRmix(x = log2fc, geneIds = geneID, mu=-4) # 
    geneRanks = CRISPhieRmixResults$FDR
  }
}
#}

#if (dataset == "TKO"){
# CRISPhieRmixResults = CRISPhieRmix::CRISPhieRmix(x = log2fc, geneIds = geneID, negCtrl = control2fc)
#  print("________________TKO______________")

#}

#geneRanks = CRISPhieRmixResults$negFDR
genes = CRISPhieRmixResults$genes 
final_data <- data.frame(genes = genes , ranks = geneRanks)
outputFilename = paste('/Users/sbodapati/Documents/Qi Lab/CRISPR_Benchmarking/output/','CRISPhieRmix_',f_name,'_Output.csv', sep='')
write.csv(final_data, outputFilename)

