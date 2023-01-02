# obtain the dataset from GEO library 
library(GEOquery)
my_id <- "GDS3309"
dat <- getGEO(my_id)

# extract the table containing gene expression data 
K <- dat@dataTable@table

# dimensions of the table 
dim(K)
[1] 22277    17

# Remove the first two columns which contains probe id and gene name for each row 
data <- K[,-c(1,2)]

# plot a boxplot to check if the data is normalized. If data is normalized than all columns would be on same level in the boxplot.  
boxplot(data)

# There are many large and small values in the box plot thereby we can’t have a clear picture. 
# In order to visualize the data better transform the data in to log2 
data2 <- log2(data)

# boxplot again 
boxplot(data2, main=expression(paste("Boxplot of the ", log2, " data")))

# the graph occurs to be normalized as all the columns are in same level 
# Remove low expression genes in the dataset by calculating the mean of each gene (row) 
test <- rowMeans(data)

# make it a quantile to get the genes below 25% quantile 
quantile(test, probs = c(.25, .5, .75))

# 25% quartile value is 5.733355, so remove all genes having a value below that value. 
test.rows <- which(test > 5.733355)

# make a new variable containing remaining genes
data3 <- data2[test.rows, ]

# make a label for samples with their class. First 8 samples are from normal population 
# while next 7 are from smoker group 
labels <- c(rep("control", 8), rep("smokers", 7))


# conduct a t-test as it contains only 2 factors and on only the first genes of both classes 
result <- t.test( data3[1, labels == 'control'], data3[1, labels == 'smokers'], var.equal=F )
result 

# the test statistic is -0.7 which means the controls have less average than the smokers
# it has 95% confidence interval 
# the difference between genes is not equal to 0 which means the samples are from 
# different populations. However, this assessed only from a small sample from the 
# population. We need to conduct this for all genes. 


# conduct t test for all genes 
result <- t.test( data3[c(1:7), labels == 'control'], data3[c(8:15), labels == 'smokers'], var.equal=F )
result 

# the test statistic is -2.2 which means the controls have less average than the smokers
# it has 95% confidence interval 
# the difference between genes is not equal to 0 which means the samples are from 
# different populations. 


# find p-values for all genes 
my.ttest <- function(v, labels)
{
  levels <- unique(labels)
  v1 <- v[ labels == levels[1]]
  v2 <- v[ labels == levels[2]]
  pval <- t.test(v1, v2, var.equal = F)$p.value
  pval
}

allpvalues <- apply(data2, 1, my.ttest, labels)

# adjust for multiplicity using BH method
padj <- p.adjust(allpvalues,method="BH")

# extract genes with p-values less than the threshold value 0.05 
sig.padj <- padj[padj < 0.05]
sig.padj.inds <- which(padj < 0.05)

# plot histogram to view the number of genes having p-value under threshold 
hist(sig.padj,col="darkmagenta",xlab="p-values",main="Genes with p-values under threshold",cex.main=0.9)

# change the column names corresponding to their class for all samples 
colnames(data3)[1:8] <- "c"
colnames(data3)[9:15] <- "s"

# grab all the rows of genes with under threshold p-value from data3
select.genes<- data3[c(sig.padj.inds), ]

# subset data with previously determined genes 
data3 <- data3[c(select.genes),]


# perform hierarchical cluster analysis on the data
dat.hca <- hclust(dist(t(data3),"man"),method="median")
plot(dat.hca, main = "HCA of control and smoker data with selected genes")


# the classic classifier algorithms - linear discriminant analysis (LDA)
# change colnames to class names in original data 
colnames(data2)[1:8] <- "c"
colnames(data2)[9:15] <- "s"
clas <- names(data2)


# make a data frame using the class labels and the data again 
data4 <- data.frame(clas, t(data2))

# make a training set using s
train <- rbind(data4[1:4,] , data4[9:11,])
test<- rbind(data4[5:8,], data4[11:15,])

# class names of test matrix
column1 <- test[,1]

# remove first column in test matrix
test1<- test[,-1]

# train the model using training set 
lda.dat <- lda(train$clas~.,train)

# predict the test sample 
pred.dat <- predict(lda.dat,test1)
	
# view confusion matrix 
table(pred.dat$class, column1)

# The total samples that are misclassified are 3. 
# In the first row, 2 samples from smoker got misclassified as normal and in last row 1 sample from control got misclassified as smoker. 

# plot the discriminant functions of all genes versus each other in an xy plot
plot(pred.dat$x, col = pred.dat$class, bg= as.numeric(factor(pred.dat$class)), pch=21,main='Scores on the 2 discriminant Variables')

# add legend 
legend(x = "topright", box.col = "black", box.lwd = 2 , legend=c("c", "s"),fill = c( "red", "black"))

# get the top 5 expressing genes 
tail(sort(p),5)
[1] 35915.66 36405.26 36988.29 38211.15 41425.52

which(p > 35915.65)
[1]  4419  6085  7304 12254 22241

# get least 5 expressing genes 
head(sort(p), 5)
[1] 1.722228 2.357062 2.628619 2.668868 2.728279

which(p < 2.728280)
[1]  7128 13683 15795 20264 21081

