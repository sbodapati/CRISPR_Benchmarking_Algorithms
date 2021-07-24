testCounts = read.table(file = "./input/Simulation.csv", header = TRUE, sep=',')
head(testCounts)
counts = testCounts[ ,3:6]
colData = data.frame(condition = factor(c(0, 0, 1, 1))) # 1 is condition, 0 is baseline
rownames(colData) = colnames(counts)
#install DESeq2
#install.packages('stringi', configure.args='--disable-cxx11', repos = "http://cran.us.r-project.org")
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
write.csv(log2fc,'./input/Simulation_l2fc.csv')
