args <- commandArgs()
print(args[6])
filename = paste("/Users/sbodapati/Documents/Qi Lab/CRISPR_Benchmarking/input/",toString(args[6]),".csv", sep='')
testCounts = read.table(file = filename, header = TRUE, sep=',')
counts = testCounts[ ,3:6]
colData = data.frame(condition = factor(c(0, 1, 1, 1))) # 1 is condition, 0 is baseline
rownames(colData) = colnames(counts)
# install DESeq2
#install.packages("BiocManager", repos = "http://cran.us.r-project.org")
#BiocManager::install("DESeq2")
# set up DESeq
testCountsDESeq2 = DESeq2::DESeqDataSetFromMatrix(countData = counts, 
                                                  colData = colData, 
                                                  design = ~ condition)
# compute
testCountsDESeq2 = DESeq2::DESeq(testCountsDESeq2)
# get results
testCountsDESeq2 = DESeq2::results(testCountsDESeq2)
log2fc = testCountsDESeq2$log2FoldChange
filename = paste("/Users/sbodapati/Documents/Qi Lab/CRISPR_Benchmarking/input/",toString(args[6]),"_l2fc.csv", sep='')
write.csv(log2fc,filename)
