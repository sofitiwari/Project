#title: "Code for the Machine Learning project"
#author: "Sofi Tiwari"
#date: "31/08/2019"
#R version 3.5.3 (2019-03-11)

#Loading the libraries. Before loading the libraries please install all packages
  #by typing install.packages("library name") in the R studio console
library(data.table) #for transposing data
library(ggplot2) #general plotting package in R
library(gridExtra) #plotting multiple graphs together
library(factoextra) #for plotting clusters
library(dbscan) #for density based clustering
library(FactoMineR) # for PCA
#Synapser package can be downloaded by typing in below command in the R Studio console:
  #install.packages("synapser", repos=c("http://ran.synapse.org", "http://cran.fhcrc.org"))
library(synapser) #for downloading data

#Downloaded the data from Synapse:Note, its required to create a synapse id 
#https://www.synapse.org/#!Synapse:syn4303551
synLogin('sofitiwari','mlproject')
# Obtain a pointer and download the data
syn4303551 = synGet(entity='syn4303551', downloadLocation=getwd())

#Reading the data
df = read.table(file = 'unc.edu_PANCAN_IlluminaGA_RNASeqV2.geneExp.tsv',
                sep = '\t', header = TRUE)

#Dropping the first 29 rows because the gene id's are not legible
df_mod = df[-c(1:29),]
rownames(df_mod)= c(1:nrow(df_mod))

#These ~20k genes are currently in the row. Transposing them to get as columns(features)
#and the rows will be the experiments
df_fin = transpose(df_mod)
rownames(df_fin) = colnames(df_mod)
colnames(df_fin) = df_fin[1,]
df_fin = df_fin[-1,]
df_fin = transpose(df_fin)
#Converting into a matrix for further use in plotting
df_matrix = as.matrix(sapply(df_fin,as.numeric))

#Checking for NA values. The function below counts the total NA values in each column
na_df = sapply(df_fin, function(x) sum(is.na(x)))

#no NA values
na_max = max(na_df)
summary(df_fin$V1)

#ELbow Plot
#Using elbow plot to check how many clusters would I need.
#It looks like a smaller number of cluster is not enough for these ~20k genes, 
#because the within-cluster sum of square(wss) is very high.
wss = (nrow(df_fin)-1)*sum(apply(df_fin,2,var))
for (i in 2:15) wss[i] = sum(kmeans(df_fin,
                                     centers=i,iter.max=30)$withinss)

options(scipen=3) 
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")
abline(v=10, col="red")

#K-means clustering on the data when clusters=100
set.seed(100)
k_fit1 = kmeans(df_fin, centers = 100,iter.max = 30)
#print(k_fit1)
x=as.data.frame(k_fit1$cluster)

#~50% of the data is clustered together. 100 clusters are not sufficient to visualise
  #any sort of clustering in the data
ggplot(data = x) +
  aes(x = `k_fit1$cluster`) +
  geom_histogram(bins = 30, fill = '#0c4c8a') +
  xlab("cluster") +
  theme_minimal()

#K-means clustering after reducing the number of clusters to 10. Its easier to visualise
set.seed(100)
k_fit2 = kmeans(df_fin, centers = 10,iter.max = 30)
#print(k_fit)
y=as.data.frame(k_fit2$cluster)
ggplot(y) +
  aes(x = `k_fit2$cluster`) +
  geom_histogram(bins = 30L, fill = "#0c4c8a") +
  xlab("cluster") +
  theme_minimal()

# Visualizing the 10 clusters
fviz_cluster(k_fit2, df_matrix, frame = FALSE, geom = "point")

#K-means clustering after reducing the number of clusters to 5
set.seed(100)
k_fit = kmeans(df_fin, centers = 5,iter.max = 30)
#print(k_fit)
z=as.data.frame(k_fit$cluster)
ggplot(z) +
  aes(x = `k_fit$cluster`) +
  geom_histogram(bins = 30L, fill = "#0c4c8a") +
  xlab("cluster") +
  theme_minimal()

# Visualizing 5 clusters
fviz_cluster(k_fit2, df_matrix, frame = FALSE, geom = "point",ellipse = TRUE)

#Hierarichal clustering(not feasible for huge dataset)
#does not seem very feasible to me, because it creates a distance matrix
#so computer will run out of space
#dist_mat = dist(df_fin, method = 'euclidean')

#Density based clustering
#https://en.proft.me/2017/02/3/density-based-clustering-r/
df_matrix = as.matrix(sapply(df_fin,as.numeric))
kNNdistplot(df_matrix, k=4)
abline(h=0.4, col="red")
db = dbscan(df_matrix, 0.4, 50)
hullplot(df_matrix, db$cluster)


#PCA
res.pca = PCA(df_matrix, graph = FALSE)
print(res.pca)

#Eigen values can be used to determine the number of components to keep after PCA
eig.val =  get_eigenvalue(res.pca)
#View(eig.val)
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 100))

#Contribution of the variable
fviz_pca_var(res.pca, col.var = "black")

# Contributions of variables to PC1
fviz_contrib(res.pca, choice = "var", axes = 1, top = 40)
# Contributions of variables to PC2
fviz_contrib(res.pca, choice = "var", axes = 2, top = 40)

#The real process starts now...

#Iteration 1
#Adding the k-means cluster to the df
df_fin1 = cbind(df_fin,y$`k_fit2$cluster`)
colnames(df_fin1)[820] = "k_cluster"
table(k_fit2$cluster)#9th cluster has highest population
#PCA on the most dense cluster
df_subset1 = df_fin1[which(df_fin1$k_cluster==9),]
res.pca = PCA(as.matrix(sapply(df_subset1,as.numeric)), graph = FALSE)
p1 = fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50), main = "Iteration 1")

#Defining some functions

#The cluster_subset function performs k-means clustering(cluster size=10)
#and returns the cluster number that has most number of genes
cluster_subset = function(df_subset)
{
  set.seed(100)
  iter = kmeans(df_subset[,-820], centers = 10,iter.max = 30) #k-means clustering
  df_fin = cbind(df_subset[,-820],iter$cluster)
  colnames(df_fin)[820] = "k_cluster"
  tt = table(iter$cluster) #frequency table
  max_cluster = as.numeric(names(tt[which.max(tt)])) #cluster with highest population
  df_subset = df_fin[which(df_fin$k_cluster==max_cluster),] #subset data
  return(df_subset) 
}

#The pca_plot function performs PCA on the clustered dataset and returns the plot
# explaining the percentage of variance by the first 10 principal components
pca_plot = function(df_subset,iteration)
{
  res.pca = PCA(as.matrix(sapply(df_subset,as.numeric)), graph = FALSE) #PCA
  plt = fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50), main = iteration) #plotting
  return(plt)
}

#Iterating the Clustering-->PCA-->Clustering-->PCA....process 11 times
subset2 = cluster_subset(df_subset1)
p2 = pca_plot(subset2, "Iteration 2")

subset3 = cluster_subset(subset2)
p3 = pca_plot(subset3, "Iteration 3")

subset4 = cluster_subset(subset3)
p4 = pca_plot(subset4, "Iteration 4")

subset5 = cluster_subset(subset4)
p5 = pca_plot(subset5, "Iteration 5")

subset6 = cluster_subset(subset5)
p6 = pca_plot(subset6, "Iteration 6")

subset7 = cluster_subset(subset6)
p7 = pca_plot(subset7, "Iteration 7")

subset8 = cluster_subset(subset7)
p8 = pca_plot(subset8, "Iteration 8")

subset9 = cluster_subset(subset8)
p9 = pca_plot(subset9, "Iteration 9")

subset10 = cluster_subset(subset9)
p10 = pca_plot(subset10, "Iteration 10")

subset11 = cluster_subset(subset10)
p11 = pca_plot(subset11, "Iteration 11")

grid.arrange(p1, p2, p3, nrow=3)
grid.arrange(p4, p5, p6, p7, nrow=4)
grid.arrange(p8, p9, p10, nrow=3)

#below code was just some tests to get the PPI graphs from StringDB,but did not go deeper
#code below is not a part of project

#library(STRINGdb) #connection to stringDB database
#library(igraph) #for manipulating the network
# #String DB connection (9606 is the id for humans)
# string_db = STRINGdb$new(version="10", species=9606,
#                           score_threshold=400, input_directory="" )
# full.graph = string_db$get_graph()
# nodes = vertex_attr(full.graph)
# 
# # see how many proteins do you have
# vcount(full.graph)
# 
# # find top 200 proteins with the highest degree
# top.degree.verticies = names(tail(sort(degree(full.graph)), 200))
# 
# # extract the relevant subgraph
# top.subgraph = induced_subgraph(full.graph, top.degree.verticies)
# 
# # count the number of proteins in it
# vcount(top.subgraph)
# 
# chk=as.list(df_mod$gene_id)
# class(chk)
# a=substring(chk, regexpr("|", chk) + 1)
# chk[1]
# substring("A1BG|1",1, regexpr("|", "A1BG|1"))
# 
# grep("B", "A1BG|1")
# 
# library('biomaRt')
# mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
# genes <- df$genes
# df<-df[,-4]
# G_list <- getBM(filters= "ensembl_peptide_id", attributes= c("ensembl_peptide_id","hgnc_symbol"),values=genes,mart= mart)
# merge(df,G_list,by.x="gene",by.y="ensembl_peptide_id")


