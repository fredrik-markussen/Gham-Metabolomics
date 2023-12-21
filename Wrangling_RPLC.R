#Clean script FM
library(dplyr)
library(arrow, warn.conflicts = T)


######## Read dataframe: 

df<-read.csv("./data/230314.CD_compounds_RPLC_Fredrik.csv", sep="")
str(df)


## Select name, CalcMW, RT, norm_area and peak rating, Select out group area and group CV:
df_sel<- df %>% 
  select(species=Name, tags=Tags, CalcMW=Calc..MW, RTmin=RT..min., contains(c("Norm", "Peak")), -contains("Gr"))



##  add column that have CalcMW/RT, is unique 

df_sel$CalcMWRT<-paste(df_sel$CalcMW ,df_sel$RTmin ,sep="/")


## reorder so CalcMWRT is first:

df_sel<- df_sel %>% 
  relocate(tags, CalcMWRT)


##Fill in Species based on CalcMWRT
sum(df_sel$species == "") #how many unidentified metabolites are there

df_sel$species <- ifelse(df_sel$species == "",  df_sel$CalcMWRT, df_sel$species ) #if else statement filling in MW/RT if metabolite ID is missing

sum(df_sel$tags == "") #how many tag classes
#assing level 4 class (un-IDable) as "F"
df_sel$tags <- ifelse(df_sel$tags == "",  'F', df_sel$tags ) #if else statement filling in MW/RT if metabolite ID is missing

#check to see distribution of confidence categories are. To keep in mind for later filtering. 
barplot(table(df_sel$tags), xlab='conf level category', ylab='count')

str(df_sel)


## Rename variables to see which sample is which individual

library(stringr)
#Replacing odd character strings 
df_sel<- df_sel %>% 
  rename_all(~str_replace_all(.,"\\.", "_"))%>%
  rename_all(~str_replace_all(.,"raw__F\\d+\\_", "")) %>%
  rename_all(~str_replace_all(.,".Area__sample_", ""))


#AND renaming sample location in 96-well plate to match Animal Individual ID
df_sel<- df_sel %>% 
  rename(P15=Norm__Area__220923_SG_Fredrik_RP_biological_sample_1_,
         P19=Norm__Area__220923_SG_Fredrik_RP_biological_sample_2_,
         P26=Norm__Area__220923_SG_Fredrik_RP_biological_sample_3_,
         P28=Norm__Area__220923_SG_Fredrik_RP_biological_sample_4_,
         P6=Norm__Area__220923_SG_Fredrik_RP_biological_sample_5_,
         P5=Norm__Area__220923_SG_Fredrik_RP_biological_sample_6_,         
         P9=Norm__Area__220923_SG_Fredrik_RP_biological_sample_7_,
         P22=Norm__Area__220923_SG_Fredrik_RP_biological_sample_8_,
         P21=Norm__Area__220923_SG_Fredrik_RP_biological_sample_9_,  
         P23=Norm__Area__220923_SG_Fredrik_RP_biological_sample_10_, 
         P4=Norm__Area__220923_SG_Fredrik_RP_biological_sample_11_ ,
         P14=Norm__Area__220923_SG_Fredrik_RP_biological_sample_12_ ,
         P3=Norm__Area__220923_SG_Fredrik_RP_biological_sample_13_ ,
         P18=Norm__Area__220923_SG_Fredrik_RP_biological_sample_14_ ,
         P16=Norm__Area__220923_SG_Fredrik_RP_biological_sample_15_ ,
         P20=Norm__Area__220923_SG_Fredrik_RP_biological_sample_16_ ,
         P17=Norm__Area__220923_SG_Fredrik_RP_biological_sample_17_ ,
         P24=Norm__Area__220923_SG_Fredrik_RP_biological_sample_18_ ,
         P11=Norm__Area__220923_SG_Fredrik_RP_biological_sample_19_ ,
         P27=Norm__Area__220923_SG_Fredrik_RP_biological_sample_20_,
         P8=Norm__Area__220923_SG_Fredrik_RP_biological_sample_21_ ,
         P10=Norm__Area__220923_SG_Fredrik_RP_biological_sample_22_,
         P7=Norm__Area__220923_SG_Fredrik_RP_biological_sample_23_,
         P2=Norm__Area__220923_SG_Fredrik_RP_biological_sample_24_ ,
         P1=Norm__Area__220923_SG_Fredrik_RP_biological_sample_25_ ,
         P13=Norm__Area__220923_SG_Fredrik_RP_biological_sample_26_ ,
         P12=Norm__Area__220923_SG_Fredrik_RP_biological_sample_27_ ,
         P25=Norm__Area__220923_SG_Fredrik_RP_biological_sample_28_,
         blank=Norm__Area__220923_SG_Fredrik_RP_extractionblank_2_,
         QC1=Norm__Area__220923_SG_Fredrik_RP_QCPool_1_             ,
         QC2=Norm__Area__220923_SG_Fredrik_RP_QCPool_2_            ,
         QC3=Norm__Area__220923_SG_Fredrik_RP_QCPool_3_            ,
         QC4=Norm__Area__220923_SG_Fredrik_RP_QCPool_4_            ,
         QC5=Norm__Area__220923_SG_Fredrik_RP_QCPool_5_             ,
         QC6=Norm__Area__220923_SG_Fredrik_RP_QCPool_6_            ,
         QC7=Norm__Area__220923_SG_Fredrik_RP_QCPool_7_             ,
         PR_P15=Peak_Rating__220923_SG_Fredrik_RP_biological_sample_1_ ,
         PR_P19=Peak_Rating__220923_SG_Fredrik_RP_biological_sample_2_ ,
         PR_P26=Peak_Rating__220923_SG_Fredrik_RP_biological_sample_3_ ,
         PR_P28=Peak_Rating__220923_SG_Fredrik_RP_biological_sample_4_,
         PR_P6=Peak_Rating__220923_SG_Fredrik_RP_biological_sample_5_ ,
         PR_P5=Peak_Rating__220923_SG_Fredrik_RP_biological_sample_6_,
         PR_P9=Peak_Rating__220923_SG_Fredrik_RP_biological_sample_7_ ,
         PR_P22=Peak_Rating__220923_SG_Fredrik_RP_biological_sample_8_ ,
         PR_P21=Peak_Rating__220923_SG_Fredrik_RP_biological_sample_9_ ,
         PR_P23=Peak_Rating__220923_SG_Fredrik_RP_biological_sample_10_,
         PR_P4=Peak_Rating__220923_SG_Fredrik_RP_biological_sample_11_,
         PR_P14=Peak_Rating__220923_SG_Fredrik_RP_biological_sample_12_,
         PR_P3=Peak_Rating__220923_SG_Fredrik_RP_biological_sample_13_,
         PR_P18=Peak_Rating__220923_SG_Fredrik_RP_biological_sample_14_,
         PR_P16=Peak_Rating__220923_SG_Fredrik_RP_biological_sample_15_,
         PR_P20=Peak_Rating__220923_SG_Fredrik_RP_biological_sample_16_,
         PR_P17 =Peak_Rating__220923_SG_Fredrik_RP_biological_sample_17_,
         PR_P24 =Peak_Rating__220923_SG_Fredrik_RP_biological_sample_18_,
         PR_P11 =Peak_Rating__220923_SG_Fredrik_RP_biological_sample_19_,
         PR_P27=Peak_Rating__220923_SG_Fredrik_RP_biological_sample_20_,
         PR_P8=Peak_Rating__220923_SG_Fredrik_RP_biological_sample_21_,
         PR_P10=Peak_Rating__220923_SG_Fredrik_RP_biological_sample_22_,
         PR_P7=Peak_Rating__220923_SG_Fredrik_RP_biological_sample_23_,
         PR_P2=Peak_Rating__220923_SG_Fredrik_RP_biological_sample_24_,
         PR_P1=Peak_Rating__220923_SG_Fredrik_RP_biological_sample_25_,
         PR_P13=Peak_Rating__220923_SG_Fredrik_RP_biological_sample_26_,
         PR_P12=Peak_Rating__220923_SG_Fredrik_RP_biological_sample_27_,
         PR_P25=Peak_Rating__220923_SG_Fredrik_RP_biological_sample_28_,
         PR_blank=Peak_Rating__220923_SG_Fredrik_RP_extractionblank_2_   ,
         PR_QC1=Peak_Rating__220923_SG_Fredrik_RP_QCPool_1_           ,
         PR_QC2=Peak_Rating__220923_SG_Fredrik_RP_QCPool_2_            ,
         PR_QC3=Peak_Rating__220923_SG_Fredrik_RP_QCPool_3_            ,
         PR_QC4=Peak_Rating__220923_SG_Fredrik_RP_QCPool_4_           ,
         PR_QC5=Peak_Rating__220923_SG_Fredrik_RP_QCPool_5_            ,
         PR_QC6=Peak_Rating__220923_SG_Fredrik_RP_QCPool_6_           ,
         PR_QC7=Peak_Rating__220923_SG_Fredrik_RP_QCPool_7_
  )


str(df_sel)


#find index range to perform peak filtering on
str_which(names(df_sel), "PR")


# filter out only metabolites that have at least 4 samples with > 5 in peakrating (Already filtered by SG in CD, shuld return unchanged)
df_PR <- df_sel%>%
  filter(rowSums(.[,41:76] >=5, na.rm = T) >=4) #rowsums switches to counts if we put a condition behind the selection (>=5). Then we filter on the basis if that count is >=4


#inspect structure:
str(df_PR)

# Write file to have it. 
#write_delim(df_PR, "./data/df_RPLC_peakrated_20221006_FM.txt")

#Consider filtering on tags at this point. 
#Filering to keep MSI level 1-3 (A-E). 
df_PR <- df_PR %>%
  filter(tags==c('A','B','C','D'))

barplot(table(df_PR$tags), xlab='conf level category', ylab='count')


#Then need to select columns only required for Metaboanalyst, Species, Norm_area, group info:
df_MA <- df_PR %>% 
  select("CalcMWRT" , contains("P"), contains ("QC"),  -species, -contains("PR"), blank)


#inspect structure:
str(df_MA)

# assign groups to dataset: 
#Info about groups:
#./data/"Randomization ID HILIC_NEG RPLC pos.xlsx"

# metaboanalyst are not offering flexibility of changing order of apperance, so need to do this by editing group names... 
group<- factor(c("CalcMWRT", 
                 'A2',
                 'IBE2',
                 'IBE1',
                 'ENT',
                 'ENT',
                 'ENT',
                 'A1',
                 'IBE1',
                 'IBE2',
                 'IBE1',
                 'ENT',
                 'A2',
                 'ENT',
                 'IBE2',
                 'A2',
                 'IBE2',
                 'IBE2',
                 'IBE1',
                 'A1',
                 'ENT',
                 'A1',
                 'A1',
                 'A1',
                 'ENT',
                 'Test',
                 'A2',
                 'A2',
                 'IBE1',
                rep('QC',7),
                'blank'
))

# add a row below column names as group info for metaboanalyst to recognize 

df_MA[1,]=group
names(df_MA)[names(df_MA) == "CalcMWRT"] <- "Sample"

#check if saples have duplicates length and length unique need to be the same. if miss match => make.unique
summary(df_MA[,1])
summary(unique(df_MA[,1])) 


#making unique
df_MA[,1] <- make.unique(df_MA[,1])

#write file
readr::write_delim(df_MA[,-26], "data/df_RPLC-pos_metaboanalyst_20220311_FM.txt")



##Convert to long format 
df_long <- df_PR %>% 
  select(species , contains("P"),contains ("qc"), -contains("PR"), blank)

#store ID as vector
ids <- names(df_long[,-1])

#Melt into long format and integrety check
melt_df <- reshape2::melt(df_long, id.vars = c("species"))
melt_df$variable <- as.character(melt_df$variable)
melt_df$value <- as.numeric(melt_df$value)


str(melt_df)
summary(melt_df)
unique(melt_df$variable)

#give Sample group ID
melt_df <- melt_df%>%
  mutate(group = case_when(
    variable == "P15"  ~ 'A2',
    variable == "P19" ~ 'IBE2',
    variable == "P26"  ~'IBE1',
    variable == "P28"  ~'ENT',
    variable == "P6"  ~'ENT',
    variable == "P5"  ~'ENT',
    variable == "P9"  ~'A1',
    variable == "P22"  ~ 'IBE1',
    variable == "P21"  ~ 'IBE2',
    variable == "P23"  ~ 'IBE1',
    variable == "P4"  ~'ENT',
    variable == "P14"  ~'A2',
    variable == "P3"  ~'ENT',
    variable == "P18"  ~'IBE2',
    variable == "P16"  ~'A2',
    variable == "P20"  ~'IBE2',
    variable == "P17"  ~'IBE2',
    variable == "P24"  ~'IBE1',
    variable == "P11"  ~'A1',
    variable == "P27"  ~'ENT',
    variable == "P8"  ~'A1',
    variable == "P10"  ~'A1',
    variable == "P7"  ~'A1',
    variable == "P2"  ~'ENT',
    variable == "P1"  ~'Test',
    variable == "P13"  ~ 'A2',
    variable == "P12"  ~'A2',
    variable == "P25"  ~'IBE1',
    variable == "QC1"  ~ "qc",
    variable == "QC2"  ~ "qc",
    variable == "QC3"  ~ "qc",
    variable == "QC4"  ~ "qc",
    variable == "QC5"  ~ "qc",
    variable == "QC6"  ~ "qc",
    variable == "QC7"  ~ "qc",
    variable == "blank"  ~ "blank"
  ))

#check is struckure ok:
unique(melt_df$group)
sum(is.na(melt_df$group))

#renam cols
names(melt_df)<- c("species","sample", "norm_area","group")


#cast df to long wide format
t_df<- reshape2::dcast(melt_df, sample+group~ species, value.var = "norm_area", fun.aggregate = sum)

names(t_df)<-iconv(names(t_df), "UTF-8", "ASCII", sub = "", ) #removes all specal characters from df

#write to disk
write.csv(t_df, "./data/long_format_RPCL_FM_20221006.csv")
