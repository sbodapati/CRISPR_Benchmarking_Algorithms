

f_name= 'CRISPhieRmix_Simulation'

filename = paste("/Users/sbodapati/Documents/Qi Lab/Benchmarking/input/",f_name,".csv", sep='')
testCounts = read.table(file = filename, header = TRUE, sep=',')

l2f_filename = paste("/Users/sbodapati/Documents/Qi Lab/Benchmarking/input/",f_name,"_l2fc",".csv", sep='')
l2fc = read.table(file = l2f_filename, header = TRUE, sep=',')

controlGuideName = 'Simulation_l2fc_control.csv'
controlGuides_filename = paste("/Users/sbodapati/Documents/Qi Lab/Benchmarking/input/",controlGuideName, sep='')
controlGuides = read.table(file = controlGuides_filename, header = TRUE, sep=',')

log2fc = l2fc$x
control2fc = controlGuides$x

print("________________here_New______________")
#install.packages('mixtools')
#install.packages('devtools')
#devtools::install_github("timydaley/CRISPhieRmix")
print(controlGuides_filename)
print(controlGuides)
print(control2fc)
print(tail(log2fc))
print("________________here_New1______________")
CRISPhieRmixResults = CRISPhieRmix::CRISPhieRmix(x = log2fc, geneIds = testCounts$gene_id, negCtrl = control2fc, PLOT = TRUE, VERBOSE = TRUE) # 
print("________________SIM______________")


geneRanks = CRISPhieRmixResults$FDR
genes = CRISPhieRmixResults$genes 
final_data <- data.frame(genes = CRISPhieRmixResults$genes , ranks = CRISPhieRmixResults$FDR)
outputFilename = paste('/Users/sbodapati/Documents/Qi Lab/Benchmarking/output/','CRISPhieRmix_',f_name,'_Output.csv', sep='')
write.csv(final_data, outputFilename)

