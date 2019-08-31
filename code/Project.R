
#Loading the libraries
library(data.table) #for transposing data
library(STRINGdb) #connection to stringDB database
library(igraph) #for manipulating the network

#Downloaded the data from Synapse:Note this might require a synapse id to register
#https://www.synapse.org/#!Synapse:syn4301332

#Reading the data
df = read.table(file = 'unc.edu_PANCAN_IlluminaGA_RNASeqV2.geneExp.tsv', sep = '\t', header = TRUE)

#Dropping the first 29 rows because the gene id's are not legible
df_mod = df[-c(1:29),]
rownames(df_mod)= c(1:nrow(df_mod))

#These 20k genes are currently in the row. Transposing them to get as columns(features)
#and the rows will be the experiments
df_fin = transpose(df_mod)
rownames(df_fin) = colnames(df_mod)
colnames(df_fin) = df_fin[1,]
df_fin = df_fin[-1,]
df_fin = transpose(df_fin)

#Checking for NA values. The function below counts the total NA values in each column
na_df = sapply(df_fin, function(x) sum(is.na(x)))

#no NA values
na_max = max(na_df)

#Was trying to summarize the info but maybe later
# options(scipen = 999)
# summary_df = do.call(cbind, lapply(df_fin[, 2:ncol(df_fin)], summary))
# summary_df_t = as.data.frame(round(t(summary_df),0))
# names(summary_df_t)[7] = paste("Missing_values")
# summary_df_t_2 = summary_df_t %>% 
#   mutate(obs = nrow(df_fin),Missing_prop = Missing_values / obs)
# print(summary_df_t_2)


#K-means- creating elbow curve to find the number of clusters that I will need
wss = (nrow(df_fin)-1)*sum(apply(df_fin,2,var))

for (i in 2:15) wss[i] <- sum(kmeans(df_fin,
                                     centers=i,iter.max=30)$withinss)

plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")


for (i in 2:50) wss[i] <- sum(kmeans(df_fin,
                                     centers=i,iter.max=30)$withinss)

plot(1:50, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

k_fit = kmeans(t(df_fin), centers = 100)
print(k_fit)

x=as.data.frame(k_fit$cluster)

#Hierarichal clustering

hc = hclust(df_fin)


#Connection to StringDB

string_db = STRINGdb$new(version="10", species=9606,
                          score_threshold=400, input_directory="" )
full.graph = string_db$get_graph()

# see how many proteins do you have    
vcount(full.graph)

# find top 200 proteins with the highest degree
top.degree.verticies = names(tail(sort(degree(full.graph)), 200))

# extract the relevant subgraph
top.subgraph = induced_subgraph(full.graph, top.degree.verticies)

# count the number of proteins in it
vcount(top.subgraph)
