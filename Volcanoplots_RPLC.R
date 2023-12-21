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
df.t <- read.csv("./data/long_format_RPCL_FM_20221006.csv", header = TRUE)



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




met.sc <- scale(sqrt(met.sc[,4:ncol(met.sc)]), center = F, scale = T)



#Making 4 plots showing effects of normalization
m.sc.plot <- cbind(df.t[2], met.sc)%>%
  melt(., id= "sample")%>%
  mutate(value =as.numeric(unlist(value)))%>%
  ggplot(aes(value, sample))+
  geom_boxplot()

m.sc.dens <- cbind(df.t[2], met.sc)%>%
  melt(., id= "sample")%>%
  mutate(value =as.numeric(unlist(value)))%>%
  ggplot(aes(value))+
  geom_density()+geom_vline(xintercept = 1, linetype= "dashed", color="gray50")

df.t.unnorm.box <-cbind(df.t[2], df.t[4:ncol(df.t)])%>%
  melt(., id= "sample")%>%
  mutate(value =as.numeric(unlist(value)))%>%
  ggplot(aes(value, sample))+
  geom_boxplot()

df.t.unnorm.dens <-cbind(df.t[2], df.t[4:ncol(df.t)])%>%
  melt(., id= "sample")%>%
  mutate(value =as.numeric(unlist(value)))%>%
  ggplot(aes(value))+
  geom_density()

#Drawing the 4 plots together
grid.arrange(df.t.unnorm.dens,m.sc.dens, df.t.unnorm.box, m.sc.plot, ncol=2)


df.2 <- cbind(df.t[,2:3], met.sc)




######ENT VS A1#####

#Vulcanoplot between grups: 
df.ent<-df.2 %>% filter(group == 'ENT')
df.a1 <- df.2 %>% filter(group == 'A1')

#species means
df.ent_means <- as.data.frame(colMeans(df.ent[, -c(1:2)], na.rm = TRUE))
names(df.ent_means)[1] <- "mean"
df.a1_means <- as.data.frame(colMeans(df.a1[, -c(1:2)], na.rm = TRUE))
names(df.a1_means)[1] <- "mean"

summary(df.ent_means)
sum(df.ent_means$mean==0)
summary(df.a1_means)
sum(df.ent_means$mean==0)


#log2 fold change
log2fc <- log2(df.a1_means$mean / df.ent_means$mean)

#calc the p vals
p.val <- sapply(colnames(df.2)[-c(1:2)], function(compound) {
  t.test(df.ent[[compound]], df.a1[[compound]], alternative = "two.sided", var.equal = TRUE)$p.value
})

p.val.adj <- p.adjust(p.val, method = "BH")


volcano<- data.frame(
  Compound = colnames(df.2)[-c(1:2)],
  log2FC = log2fc,
  negLog10P = -log10(p.val.adj)
)
library(ggrepel)


# Create the volcano plot
ggplot(volcano, aes(x = log2FC, y = negLog10P)) +
  geom_point(aes(color = negLog10P)) +
  scale_y_continuous(limits= c(0,  max(volcano$negLog10P*1.2, na.rm = TRUE))) +
  scale_x_continuous(limits= c(-4,  4)) +
  labs(x = "Log2FC", y = "-Log10 P-value (Adjusted)",title = expression(atop("RPLC-Pos Log2FC", "ENT vs A1")))  +
  scale_color_gradient(low = "steelblue4", high = "firebrick1") +
  #geom_hline(yintercept = 5, linetype = "dashed", color = "darkred") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkred") +
  geom_vline(xintercept = 0, linetype= "dashed", color="gray50")+
  #geom_vline(xintercept = -2, linetype= "dashed", color="gray50")+
  #geom_vline(xintercept = 2, linetype= "dashed", color="gray50")+
  theme_classic() +
  theme(
    text = element_text(family = "arial",size = 16, color = 'gray20'),
    legend.position = "none",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, vjust = 2.5),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line = element_line(color = "black"))


volcano %>%
  filter(negLog10P > 5, (log2FC < -2 | log2FC > 2)) %>%
  ggplot(aes(x = 1, y = negLog10P, label = Compound)) +
  geom_point() +
  geom_text_repel(
    force        = 0.5,
    nudge_x      = 0.01,
    direction    = "y",
    hjust        = 0.5,
    segment.size = 0.2
  )  +
  theme_classic() +
  labs(x = "", y = "-Log10 P-value (Adjusted)") +
  theme(
    plot.margin = margin(1, 1, 1, 1, "cm"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )



#####IBE1 VS ENT#####

#Vulcanoplot between grups: 
df.ent<-df.2 %>% filter(group == 'ENT')
df.IBE1 <- df.2 %>% filter(group == 'IBE1')

#species means
df.ent_means <- as.data.frame(colMeans(df.ent[, -c(1:2)], na.rm = TRUE))
names(df.ent_means)[1] <- "mean"

df.IBE1_means <- as.data.frame(colMeans(df.IBE1[, -c(1:2)], na.rm = TRUE))
names(df.IBE1_means)[1] <- "mean"



summary(df.ent_means)
sum(df.ent_means$mean==0)
summary(df.IBE1_means)
sum(df.ent_means$mean==0)


#log2 fold change
log2fc <- log2(df.IBE1_means$mean / df.ent_means$mean)



#calc the p vals
p.val <- sapply(colnames(df.2)[-c(1:2)], function(compound) {
  t.test(df.ent[[compound]], df.a1[[compound]], alternative = "two.sided", var.equal = TRUE)$p.value
})

p.val.adj <- p.adjust(p.val, method = "BH")


volcano<- data.frame(
  Compound = colnames(df.2)[-c(1:2)],
  log2FC = log2fc,
  negLog10P = -log10(p.val.adj)
)



# Create the volcano plot
ggplot(volcano, aes(x = log2FC, y = negLog10P)) +
  geom_point(aes(color = negLog10P)) +
  scale_y_continuous(limits= c(0,  max(volcano$negLog10P*1.2, na.rm = TRUE))) +
  scale_x_continuous(limits= c(-3,  3)) +
  labs(x = "Log2FC", y = "-Log10 P-value (Adjusted)",title = expression(atop("HILIC-neg Log2FC", "ENT vs IBE1")))  +
  scale_color_gradient(low = "steelblue4", high = "firebrick1") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkred") +
  geom_vline(xintercept = 0, linetype= "dashed", color="gray50")+
  theme_classic() +
  theme(
    text = element_text(family = "arial",size = 16, color = 'gray20'),
    legend.position = "none",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, vjust = 2.5),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line = element_line(color = "black"))



library(ggrepel)
volcano %>%
  filter(negLog10P > 5, (log2FC < -1 | log2FC > 1)) %>%
  ggplot(aes(x = 1, y = negLog10P, label = Compound)) +
  geom_point() +
  geom_text_repel(
    force        = 0.5,
    nudge_x      = 0.01,
    direction    = "y",
    hjust        = 0.5,
    segment.size = 0.2
  )  +
  theme_classic() +
  labs(x = "", y = "-Log10 P-value (Adjusted)") +
  theme(
    plot.margin = margin(1, 1, 1, 1, "cm"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )


#####IBE1 VS A1#####

#Vulcanoplot between grups: 

#Vulcanoplot between grups: 
df.IBE1<-df.2 %>% filter(group == 'IBE1')
df.a1 <- df.2 %>% filter(group == 'A1')

#species means
df.IBE1_means <- as.data.frame(colMeans(df.IBE1[, -c(1:2)], na.rm = TRUE))
names(df.IBE1_means)[1] <- "mean"
df.a1_means <- as.data.frame(colMeans(df.a1[, -c(1:2)], na.rm = TRUE))
names(df.a1_means)[1] <- "mean"

summary(df.IBE1_means)
sum(df.IBE1_means$mean==0)
summary(df.a1_means)
sum(df.IBE1_means$mean==0)


#log2 fold change
log2fc <- log2(df.a1_means$mean / df.IBE1_means$mean)

#calc the p vals
p.val <- sapply(colnames(df.2)[-c(1:2)], function(compound) {
  t.test(df.IBE1[[compound]], df.a1[[compound]], alternative = "two.sided", var.equal = TRUE)$p.value
})

p.val.adj <- p.adjust(p.val, method = "BH")


volcano<- data.frame(
  Compound = colnames(df.2)[-c(1:2)],
  log2FC = log2fc,
  negLog10P = -log10(p.val.adj)
)
library(ggrepel)


# Create the volcano plot
ggplot(volcano, aes(x = log2FC, y = negLog10P)) +
  geom_point(aes(color = negLog10P)) +
  scale_y_continuous(limits= c(0,  max(volcano$negLog10P*1.2, na.rm = TRUE))) +
  scale_x_continuous(limits= c(-4,  4)) +
  labs(x = "Log2FC", y = "-Log10 P-value (Adjusted)",title = expression(atop("HILIC-neg Log2FC", "IBE1 vs A1")))  +
  scale_color_gradient(low = "steelblue4", high = "firebrick1") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkred") +
  geom_vline(xintercept = 0, linetype= "dashed", color="gray50")+
  theme_classic() +
  theme(
    text = element_text(family = "arial",size = 16, color = 'gray20'),
    legend.position = "none",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, vjust = 2.5),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line = element_line(color = "black"))

volcano %>%
  filter(negLog10P > 2, (log2FC < -1 | log2FC > 1)) %>%
  ggplot(aes(x = 1, y = negLog10P, label = Compound)) +
  geom_point() +
  geom_text_repel(
    force        = 0.5,
    nudge_x      = 0.01,
    direction    = "y",
    hjust        = 0.5,
    segment.size = 0.2
  )  +
  theme_classic() +
  labs(x = "", y = "-Log10 P-value (Adjusted)") +
  theme(
    plot.margin = margin(1, 1, 1, 1, "cm"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )


#####IBE2 VS ENT#####

#Vulcanoplot between grups: 
df.ent<-df.2 %>% filter(group == 'ENT')
df.IBE2 <- df.2 %>% filter(group == 'IBE2')

#species means
df.ent_means <- as.data.frame(colMeans(df.ent[, -c(1:2)], na.rm = TRUE))
names(df.ent_means)[1] <- "mean"

df.IBE2_means <- as.data.frame(colMeans(df.IBE2[, -c(1:2)], na.rm = TRUE))
names(df.IBE2_means)[1] <- "mean"



summary(df.ent_means)
sum(df.ent_means$mean==0)
summary(df.IBE2_means)
sum(df.ent_means$mean==0)


#log2 fold change
log2fc <- log2(df.IBE2_means$mean / df.ent_means$mean)



#calc the p vals
p.val <- sapply(colnames(df.2)[-c(1:2)], function(compound) {
  t.test(df.ent[[compound]], df.a1[[compound]], alternative = "two.sided", var.equal = TRUE)$p.value
})

p.val.adj <- p.adjust(p.val, method = "BH")


volcano<- data.frame(
  Compound = colnames(df.2)[-c(1:2)],
  log2FC = log2fc,
  negLog10P = -log10(p.val.adj)
)



# Create the volcano plot
ggplot(volcano, aes(x = log2FC, y = negLog10P)) +
  geom_point(aes(color = negLog10P)) +
  scale_y_continuous(limits= c(0,  max(volcano$negLog10P*1.2, na.rm = TRUE))) +
  scale_x_continuous(limits= c(-3,  3)) +
  labs(x = "Log2FC", y = "-Log10 P-value (Adjusted)",title = expression(atop("HILIC-neg Log2FC", "ENT vs IBE2")))  +
  scale_color_gradient(low = "steelblue4", high = "firebrick1") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkred") +
  geom_vline(xintercept = 0, linetype= "dashed", color="gray50")+
  theme_classic() +
  theme(
    text = element_text(family = "arial",size = 16, color = 'gray20'),
    legend.position = "none",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, vjust = 2.5),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line = element_line(color = "black"))



library(ggrepel)
volcano %>%
  filter(negLog10P > 5, (log2FC < -1 | log2FC > 1)) %>%
  ggplot(aes(x = 1, y = negLog10P, label = Compound)) +
  geom_point() +
  geom_text_repel(
    force        = 0.5,
    nudge_x      = 0.01,
    direction    = "y",
    hjust        = 0.5,
    segment.size = 0.2
  )  +
  theme_classic() +
  labs(x = "", y = "-Log10 P-value (Adjusted)") +
  theme(
    plot.margin = margin(1, 1, 1, 1, "cm"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

#####IBE2 VS IBE1#####

#Vulcanoplot between grups: 
df.IBE1<-df.2 %>% filter(group == 'IBE1')
df.IBE2 <- df.2 %>% filter(group == 'IBE2')

#species means
df.IBE1_means <- as.data.frame(colMeans(df.IBE1[, -c(1:2)], na.rm = TRUE))
names(df.IBE1_means)[1] <- "mean"

df.IBE2_means <- as.data.frame(colMeans(df.IBE2[, -c(1:2)], na.rm = TRUE))
names(df.IBE2_means)[1] <- "mean"



summary(df.IBE1_means)
sum(df.IBE1_means$mean==0)
summary(df.IBE2_means)
sum(df.IBE1_means$mean==0)


#log2 fold change
log2fc <- log2(df.IBE1_means$mean / df.IBE2_means$mean)



#calc the p vals
p.val <- sapply(colnames(df.2)[-c(1:2)], function(compound) {
  t.test(df.IBE1[[compound]], df.a1[[compound]], alternative = "two.sided", var.equal = TRUE)$p.value
})

p.val.adj <- p.adjust(p.val, method = "BH")


volcano<- data.frame(
  Compound = colnames(df.2)[-c(1:2)],
  log2FC = log2fc,
  negLog10P = -log10(p.val.adj)
)



# Create the volcano plot
ggplot(volcano, aes(x = log2FC, y = negLog10P)) +
  geom_point(aes(color = negLog10P)) +
  scale_y_continuous(limits= c(0,  max(volcano$negLog10P*1.2, na.rm = TRUE))) +
  scale_x_continuous(limits= c(-3,  3)) +
  labs(x = "Log2FC", y = "-Log10 P-value (Adjusted)",title = expression(atop("HILIC-neg Log2FC", "IBE2 vs IBE1")))  +
  scale_color_gradient(low = "steelblue4", high = "firebrick1") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkred") +
  geom_vline(xintercept = 0, linetype= "dashed", color="gray50")+
  theme_classic() +
  theme(
    text = element_text(family = "arial",size = 16, color = 'gray20'),
    legend.position = "none",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, vjust = 2.5),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line = element_line(color = "black"))



library(ggrepel)
volcano %>%
  filter(negLog10P > 5, (log2FC < -1 | log2FC > 1)) %>%
  ggplot(aes(x = 1, y = negLog10P, label = Compound)) +
  geom_point() +
  geom_text_repel(
    force        = 0.5,
    nudge_x      = 0.01,
    direction    = "y",
    hjust        = 0.5,
    segment.size = 0.2
  )  +
  theme_classic() +
  labs(x = "", y = "-Log10 P-value (Adjusted)") +
  theme(
    plot.margin = margin(1, 1, 1, 1, "cm"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )




#####A2 VS A1#####

#Vulcanoplot between grups: 

#Vulcanoplot between grups: 
df.A2<-df.2 %>% filter(group == 'A2')
df.a1 <- df.2 %>% filter(group == 'A1')

#species means
df.A2_means <- as.data.frame(colMeans(df.A2[, -c(1:2)], na.rm = TRUE))
names(df.A2_means)[1] <- "mean"
df.a1_means <- as.data.frame(colMeans(df.a1[, -c(1:2)], na.rm = TRUE))
names(df.a1_means)[1] <- "mean"

summary(df.A2_means)
sum(df.A2_means$mean==0)
summary(df.a1_means)
sum(df.A2_means$mean==0)


#log2 fold change
log2fc <- log2(df.A2_means$mean / df.a1_means$mean)

#calc the p vals
p.val <- sapply(colnames(df.2)[-c(1:2)], function(compound) {
  t.test(df.A2[[compound]], df.a1[[compound]], alternative = "two.sided", var.equal = TRUE)$p.value
})

p.val.adj <- p.adjust(p.val, method = "BH")


volcano<- data.frame(
  Compound = colnames(df.2)[-c(1:2)],
  log2FC = log2fc,
  negLog10P = -log10(p.val.adj)
)



# Create the volcano plot
ggplot(volcano, aes(x = log2FC, y = negLog10P)) +
  geom_point(aes(color = negLog10P)) +
  scale_y_continuous(limits= c(0,  max(volcano$negLog10P*1.2, na.rm = TRUE))) +
  scale_x_continuous(limits= c(-4,  4)) +
  labs(x = "Log2FC", y = "-Log10 P-value (Adjusted)",title = expression(atop("HILIC-neg Log2FC", "A1 vs A2")))  +
  scale_color_gradient(low = "steelblue4", high = "firebrick1") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkred") +
  geom_vline(xintercept = 0, linetype= "dashed", color="gray50")+
  theme_classic() +
  theme(
    text = element_text(family = "arial",size = 16, color = 'gray20'),
    legend.position = "none",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, vjust = 2.5),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line = element_line(color = "black"))

volcano %>%
  filter(negLog10P > 2, (log2FC < -0.5 | log2FC > 0-5)) %>%
  ggplot(aes(x = 1, y = negLog10P, label = Compound)) +
  geom_point() +
  geom_text_repel(
    force        = 0.5,
    nudge_x      = 0.01,
    direction    = "y",
    hjust        = 0.5,
    segment.size = 0.2
  )  +
  theme_classic() +
  labs(x = "", y = "-Log10 P-value (Adjusted)") +
  theme(
    plot.margin = margin(1, 1, 1, 1, "cm"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
