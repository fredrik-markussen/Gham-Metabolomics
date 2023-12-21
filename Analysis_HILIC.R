#install.packages(c('ggalt','tidyverse', 'reshape2', 'circlize', 'ca', 'gridExtra' ))

library(dplyr)
#library(arrow)
library(gridExtra)
library(reshape2)
library(ca)
library(vegan)
library(circlize)
library(RColorBrewer)
library(ggalt)


#install complex heatmap if required . follow instuctions in console (remove "#")
#    if (!require("BiocManager", quietly = TRUE))
#       install.packages("BiocManager")

#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap) 

#BiocManager::install('PCAtools')
library(PCAtools)

#Import data
df.t <- read.csv("./data/long_format_HILIC-neg_FM_20231218.csv", header = TRUE)


ncol(df.t)

#Normalization porcess:
#normalization against the QC samples using the raw values. 
#aims to correct for batch effects, instrumental drift, systematic variations -> systemic variations are minimized

#log transformation to QC-normalized data.
#Log transformation stabilizes the variance across the range of data, making the data more suitable for subsequent analysis.
#It also helps to make the data more normally distributed, which is an assumption of many statistical techniques.
#Scale and center the data:
#Helps issues where some variables dominate solely due to their scale.
#Centering (subtracting the mean) can be particularly useful for techniques like PCA, where you're interested in the variation around the mean.


qc <- df.t[df.t$group == "qc", 4:ncol(df.t)] #extracts qc grpup

norm <- apply(qc, 2, median, na.rm = TRUE) #calculates the median for each compound

# Apply normalization factor to each compound
met.sc <- df.t
met.sc[,4:ncol(met.sc)] <- sweep(df.t[, 4:ncol(met.sc)], 2, norm, FUN = "/") #sweep across all columns, 2 indicates columns




met.sc <- scale(sqrt(met.sc[,4:ncol(met.sc)]), center = T, scale = T)
#view(met.sc)

#Making 4 plots showing effects of normalization
m.sc.plot <- cbind(df.t[2], met.sc)%>%
  melt(., id= "sample")%>%
  mutate(value =as.numeric(unlist(value)))%>%
  ggplot(aes(value, variable))+
  geom_boxplot(fill='gray20', color= 'gray50') +
  labs(y='Compound',
       x='Normalized Peak Intensity') + 
  theme_classic()+
  theme(
    axis.text.y = element_blank(),   
    axis.ticks.y = element_blank(),
    #plot.margin = margin(1, 1, 1, 1, "cm")
  )


m.sc.dens <- cbind(df.t[2], met.sc)%>%
  melt(., id= "sample")%>%
  mutate(value =as.numeric(unlist(value)))%>%
  ggplot(aes(value))+
  geom_density(linewidth=1, color ='grey20')+
  labs(title = 'After Normalization', 
       y='Density',
       x='Normalized Peak Intensity') + 
  theme_classic()+
  theme(
    axis.ticks.y = element_blank(),
    plot.title = element_text(hjust = 0.5),
    #plot.margin = margin(1, 1, 1, 1, "cm")
  )


df.t.unnorm.box <-
  cbind(df.t[2], df.t[4:ncol(df.t)])%>%
  melt(., id= "sample")%>%
  mutate(value =as.numeric(unlist(value)))%>%
  ggplot(aes(value, variable))+
  geom_boxplot(fill='gray20', color= 'gray50') +
  labs(y='Compound',
       x='Normalized Peak Intensity') + 
  theme_classic()+
  theme(
    axis.text.y = element_blank(),   
    axis.ticks.y = element_blank(),
    #plot.margin = margin(1, 1, 1, 1, "cm")
  )

df.t.unnorm.dens <-cbind(df.t[2], df.t[4:ncol(df.t)])%>%
  melt(., id= "sample")%>%
  mutate(value =as.numeric(unlist(value)))%>%
  ggplot(aes(value))+
  geom_density(linewidth=1, color ='grey20')+
  labs(title = 'Before Normalization', 
       y='Density',
       x='Normalized Peak Intensity') + 
  theme_classic()+
  theme(
    axis.ticks.y = element_blank(),
    plot.title = element_text(hjust = 0.5),
    #plot.margin = margin(1, 1, 1, 1, "cm")
  )


#Drawing the 4 plots together
grid.arrange(df.t.unnorm.dens,m.sc.dens, df.t.unnorm.box, m.sc.plot, ncol=2)


#Inspection of plot the normalization seems OK.

##############################################
##############################################


#Cluster analysis with heat map
h.df<-as.matrix(met.sc)
rownames(h.df) <- df.t[,2] #Making row names the samples
h.df<-as.matrix(h.df[-1,]) #removed blank 



##Heatmap preamble
c <-rev(brewer.pal(6,"RdBu")) #take reverse order of six colors from Brewer palette RdBu
col_fun = colorRamp2(c(-2,-1,0,1,2,4),c) #make gradient and supply color scale range

###define color for the groups with hexadecimal keys
an.col <- list(group = c("IBE1"= "#E41A1C", "IBE2"="#377EB8", "A1"="#4DAF4A", 
                         "A2"="#984EA3","ENT"= "#FF7F00",
                         "Test" ="#FFFF33","qc"="#FDBF6F", "blank"="#FDBF6A" ))

###Define vector that is the groups (-blank) 
ha <- HeatmapAnnotation(group = df.t[-1,3], col = an.col,
                        annotation_legend_param = list(
                          group = list(title = ""))
)


##Make heatmap
heat <- Heatmap(t(h.df),
                show_row_names=FALSE,
                clustering_distance_rows ="euclidean",
                clustering_distance_columns ="euclidean",
                clustering_method_rows= "ward.D",
                clustering_method_columns= "ward.D",
                column_dend_height = unit(2, "cm"), 
                row_dend_width = unit(2, "cm"),
                col = col_fun,
                top_annotation = ha,
                heatmap_legend_param = list(
                  title = "Relative abundance",
                  legend_height = unit(3, "cm"),
                  title_position = "leftcenter-rot"
                  
                ))

##Draw heatmap
heat

#k-means clustering

##Load required packages
library(FactoMineR)
library(factoextra)
library(cluster)


#inspect head of matrix to cluster. We want to use the scaled data for clustering.
head(h.df[,1:4])


#Don't know how many clusters is best so we make and elbow plot. 
fviz_nbclust(h.df, kmeans, method = "wss")


#Not totally clear what is the best so we may look at total intra-cluster variation
##for different values of k with their expected values for a distribution with no clustering

fviz_gap_stat(clusGap(h.df,
                      FUN = kmeans,
                      nstart = 25,
                      K.max = 10,
                      B = 50)
)

#looks like 9 clusters is optimal Which is a bit much even when including qc and tes animal (7)
#Perform k-means clustering with 7 clusters

km.h.df <- kmeans(h.df, centers = 7)

#Prints alot but maybe this is important
#Within cluster sum of squares by cluster:
#(between_SS / total_SS =  75.0 %)

#so we plot it
fviz_cluster(km.h.df, data = h.df, max.overlaps=70,
             ellipse = FALSE,
             star.plot=TRUE,#Addsegmentsfromcentroidstoitems 
             repel=TRUE,#Avoidlabeloverplotting(slow) 
             ggtheme=theme_minimal() )





#very clear separation of the early arousal groups but dont know how to compare against
#our groups to the clustering. 
#Pooled QC samples looks good (perfect overlap). 
#P1 (test animal, non hibernator) driving second component of k-means clustering. 
#Removing these two from analysis. 

##Heatmap preamble
h.df.2 <- h.df[2:28,] #removing P1 and QCs
#h.df.2 <-h.df.2[-11,]


c <-rev(brewer.pal(6,"RdBu")) #take reverse order of six colors from Brewer palette RdBu
col_fun = colorRamp2(c(-2,-1,0,1,2,4),c) #make gradient and supply color scale range

###define color for the groups with hexadecimal keys
an.col <- list(group = c("IBE1"= "#E41A1C", "IBE2"="#377EB8", "A1"="#4DAF4A", 
                         "A2"="#984EA3","ENT"= "#FF7F00"
))

###Define vector that is the groups (-blank , -P1, -QC) 
ha.2 <- HeatmapAnnotation(group = df.t[3:29,3], col = an.col,
                          annotation_legend_param = list(
                            group = list(title = ""))
)


##Make heatmap
heat.2 <- Heatmap(t(h.df.2),
                  show_row_names=FALSE,
                  clustering_distance_rows ="euclidean",
                  clustering_distance_columns ="euclidean",
                  clustering_method_rows= "ward.D",
                  clustering_method_columns= "ward.D",
                  column_dend_height = unit(2, "cm"), 
                  row_dend_width = unit(2, "cm"),
                  col = col_fun,
                  top_annotation = ha.2,
                  heatmap_legend_param = list(
                    title = "Relative abundance",
                    legend_height = unit(3, "cm"),
                    title_position = "leftcenter-rot"
                    
                  ))

##Draw heatmap
heat.2

#k-means clustering withouth QC and test
fviz_nbclust(h.df.2, kmeans, method = "wss")


#looks like 3 clusters is optimal which makes sense since it matches our numbers of "overall" states: arousal, early interbout, late interbout + entry.
#Perform k-means clustering with 7 clusters

km.h.df.2 <- kmeans(h.df.2, centers = 3, iter.max = 100000, nstart = 3 )

#Prints alot but maybe this is important
##Within cluster sum of squares by cluster:
## [1]    0.0000 1766.5678 1433.7513 1493.9281  255.9023    0.0000  828.1447
###(between_SS / total_SS =  60.2 %)


#so we plot it
fviz_cluster(km.h.df.2, data = h.df.2, max.overlaps=70,
             ellipse = T,
             star.plot=F,#Addsegmentsfromcentroidstoitems 
             repel=F,#Avoidlabeloverplotting(slow) 
             ggtheme=theme_minimal() )

#interesting separation pattern. Two components seem to explain 50% of the variation within the data. 
#Need to plot the groups on top of k-means plot. 
#P27 sample had thawed once so included in analysis similar thawed sample (P28) of P5 (not thawed). 
#Since P28 and P5 looks nearly identical in heatmap and k-means cluster, 
#P27 is deemed valid sample to include further.


#since we have log10 transformed, scaled and centered our data (thus interval data)
#PCA analysis is appropriate. 

metadata<-data.frame(row.names = rownames(h.df.2))

metadata$group <- df.t[3:29,3] #!test & QC and P2


pca.h.df.2 <- pca(t(h.df.2), metadata=metadata ,removeVar = 0.1)

screeplot(pca.h.df.2, axisLabSize = 18, titleLabSize = 22)

pairsplot(pca.h.df.2,
          components = getComponents(pca.h.df.2, c(1:5)),
          triangle = TRUE, trianglelabSize = 12,
          hline = 0, vline = 0,
          pointSize = 0.75,
          gridlines.major = FALSE, gridlines.minor = FALSE,
          title = 'Pairs plot', plotaxes = FALSE,
          margingaps = unit(c(-0.01, -0.01, -0.01, -0.01), 'cm')
)

#3 components explain 67% of the variation in the data. but mostly 1 and 2. No me than 3 in loading plot is probable sensible
#

biplot(pca.h.df.2,
       colby = "group",
       colkey = c('A1' = '#1f77b4',
                  'A2' = '#ff7f0e', 
                  'IBE1' = '#2ca02c', 
                  'IBE2' = '#d62728', 
                  'ENT' = '#9467bd'),
       lab= NULL,
       pointSize = 5,
       showLoadings = TRUE,
       ntopLoadings = 5,
       lengthLoadingsArrowsFactor = 1.5,
       sizeLoadingsNames = 3,
       showLoadingsNames = T,
       hline = 0, vline = 0,
       legendPosition = 'right',
       encircle = TRUE
)



pca.loadings <-plotloadings(pca.h.df.2,
                            components = getComponents(pca.h.df.2, c(1)),
                            rangeRetain = 0.01,
                            labSize = 4.0,
                            absolute =T,
                            title = 'Loadings plot',
                            subtitle = 'Top PC',
                            caption = 'Top 1% variables',
                            shapeSizeRange = c(0.1, 16),
                            col = c('white', 'firebrick'),
                            drawConnectors = T,
) + coord_flip()


pca.loadings


################################################################################
#Sorry from here things are messy :) 

a <- as.data.frame(pca.h.df.2$loadings) %>% 
  select(PC1) %>% 
  tibble::rownames_to_column()


plot_loading <- a%>%
  filter( PC1< (min(a$PC1) + ((max(a$PC1) - min(a$PC1)) * 0.01)) | PC1>(max(a$PC1) - ((max(a$PC1) - min(a$PC1)) * 0.01))) %>%
  ggplot(aes(reorder(rowname, PC1, mean),y=PC1))+geom_point()+coord_flip()

plot_loading

b <- as.data.frame(pca.h.df.2$loadings) %>% 
  select(PC2) %>% 
  rownames_to_column()


plot_loading2 <- b%>%
  filter( PC2< (min(b$PC2) + ((max(b$PC2) - min(b$PC2)) * 0.1)) | PC2>(max(b$PC2) - ((max(b$PC2) - min(b$PC2)) * 0.1))) %>%
  ggplot(aes(reorder(rowname, PC2, mean),y=PC2))+geom_point()+scale_x_discrete(position = "top") +coord_flip()

plot_loading2



offset <- (max(a$PC1) - min(a$PC1)) * 0.1

(max(a$PC1) - ((max(a$PC1) - min(a$PC1)) * 0.1))
(min(a$PC1) + ((max(a$PC1) - min(a$PC1)) * 0.1))


((max(a$PC1)-min(a$PC1))*0.1)
a

range(a$PC1)

getLoadings(pca.h.df.2)

pca.h.df.2[["xvars"]]




pca.h.df.2 <- PCA(h.df.2,  graph = FALSE, )
get_eig(pca.h.df.2)

var <- get_pca_var(pca.h.df.2)

mx.cos2 <- as.matrix(var$cos2)



fviz_pca_var(pca.h.df.2, col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping
)

corrplot(var$contrib, is.corr=FALSE) 



# Contributions of variables to PC2
fviz_contrib(pca.h.df.2, choice = "var", axes = 1:2, top = 20)


fviz_pca_ind(pca.h.df.2,
             geom = "point",
             col.ind = df.t[2:28,3],
             addEllipses = TRUE, ellipse.level = 0.7,# color by groups
             repel = T # Avoid text overlapping (slow if many points)
)
head(var$cos2, 4)



library(corrplot)


corrplot(var$cos2, type="upper", order="hclust", 
         p.mat = p.mat, sig.level = 0.01, insig = "blank")

