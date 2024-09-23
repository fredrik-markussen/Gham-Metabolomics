library(tidyverse)
library(dplyr)
library(arrow, warn.conflicts = T)
library(stringr)


blood.hilic.df<-read.csv("./data/230314.CD_compounds_HILIC_Fredrik_csv.csv", sep="", check.names =TRUE)


name.blood.hilic<- as.data.frame(blood.hilic.df$Name)
name.blood.hilic$type <- paste("HILIC")
colnames(name.blood.hilic) <- c("Name", "plasma_type")


df<-read.csv("./data/230314.CD_compounds_RPLC_Fredrik.csv", sep="", check.names = T)
name.blood.RPLC <- as.data.frame(df$Name)
name.blood.RPLC$type <- paste("RPLC")
colnames(name.blood.RPLC) <- c("Name", "plasma_type")

plasma.names <- full_join(name.blood.hilic, name.blood.RPLC) 


#if name blank, remove row
plasma.names <- plasma.names %>% filter(!is.na(Name))


csf.RPLC <- read_delim("data/csf/20240214_maur_csf_RPLC_pos.csv", 
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE)

name.csf.RPLC <- as.data.frame(csf.RPLC$Name)
name.csf.RPLC$type <- paste("RPLC")
colnames(name.csf.RPLC) <- c("Name", "csf_type")

csf.hilic <- read_delim("data/csf/20240214_maur_csf_HILIC_neg.csv", 
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE)

name.csf.hilic <- as.data.frame(csf.hilic$Name)
name.csf.hilic$type <- paste("HILIC")
colnames(name.csf.hilic) <- c("Name", "csf_type")


csf.names <- full_join(name.csf.RPLC, name.csf.hilic)

csf.names <- csf.names %>% filter(!is.na(Name))





head(plasma.names)
head(csf.names)


#############################
#clean names function standardising naming conventions and removing unwanted characters
clean.names <- function(name) {
  if (is.na(name)) {
    return(NA)
  }
  if (name == "" || str_detect(name, "Similar to")) {
    return(NA)  # Return NA for empty rows or rows containing "similar to..."
  } else {
    name <- iconv(name, to = "UTF-8")  # Convert to UTF-8 using iconv
    name <- str_trim(name)  # Trim whitespace
    name <- str_replace_all(name, "[\\{\\}\\[\\]?�]", "")
    name <- str_replace_all(name, "-�", "")
    name <- str_to_title(name)  # Capitalize each word
    # Add more cleaning rules as needed
    return(name)
  }
}

# Escape special characters in both dataframes
plasma.names$Name <- sapply(plasma.names$Name, clean.names)
csf.names$Name <- sapply(csf.names$Name, clean.names)

# Filter out NA values that could result from escaping errors
plasma.names <- plasma.names %>% filter(!is.na(Name))
csf.names <- csf.names %>% filter(!is.na(Name))

# Exact matches
exact_matches <- inner_join(plasma.names, csf.names)

# Partial matches (case insensitive)
partial_match <- function(plasma, csf) {
  matches <- lapply(plasma$Name, function(name) {
    matched <- csf[grepl(name, csf$Name, ignore.case = TRUE), ]
    if (nrow(matched) > 0) {
      matched$plasma_name <- name
    }
    return(matched)
  })
  do.call(rbind, matches)
}

partial_matches <- partial_match(plasma.names, csf.names)


head(partial_matches)

# Combine exact and partial matches
combined_matches <- bind_rows(
  exact_matches,
  plasma.names %>%
    inner_join(partial_matches, by = c("Name" = "plasma_name")) %>%
    select(Name, plasma_type, csf_type)
) %>%
  distinct()

# View similar names
print(combined_matches)


write.csv(combined_matches, file = "data/CSF_PLASMA_combined_matches.csv", row.names = FALSE)


##########################################################
# Load necessary packages
library(dplyr)
library(stringdist)

# Define a function to compute column similarities
name_sim <- function(df1, df2, threshold = 0.70, method = "lv") {
  df1_names <- df1$Name
  df2_names <- df2$Name
  similar_names <- data.frame(Name1 = character(), Name2 = character(), Similarity = numeric(), stringsAsFactors = FALSE)
  
  for (name1 in df1_names) {
    for (name2 in df2_names) {
      similarity <- stringdist::stringdist(tolower(name1), tolower(name2), method = method)
      max_chars <- max(nchar(name1), nchar(name2))
      normalized_similarity <- 1 - (similarity / max_chars)
      if (normalized_similarity >= threshold) {
        similar_names <- rbind(similar_names, data.frame(Name1 = name1, Name2 = name2, Similarity = normalized_similarity))
      }
    }
  }
  return(similar_names)
}

# Compute similarities between plasma.names and csf.names
similar_names <- name_sim(plasma.names, csf.names, threshold = 0.70)

# Prepare data for joining
similar_names <- similar_names %>%
  rename(plasma_name = Name1, csf_name = Name2)




# Join the data frames based on similar names
partial_matches <- plasma.names %>%
  inner_join(similar_names, by = c("Name" = "plasma_name")) %>%
  inner_join(csf.names, by = c("csf_name" = "Name"))

# Combine exact and partial matches
combined_matches <- bind_rows(
  exact_matches,
  partial_matches %>%
    select(Name = csf_name, plasma_type, csf_type)
) %>%
  distinct()

# View similar names
print(combined_matches)

# Count the number of similar names
num_similar_names <- nrow(combined_matches)
print(paste("Number of similar names:", num_similar_names))




##############################################################
HS_csf <- read.csv("./data/csf/Human_csf_referenceDB.csv", header = TRUE, sep = ",")

summary(HS_csf)

hs.names <- as.data.frame(HS_csf$NAME)
hs.names$type <- paste("Ref")
colnames(hs.names) <- c("Name", "csf_type")

head(hs.names)

# Compute similarities between plasma.names and csf.names
similar_names <- name_sim(hs.names, csf.names, threshold = 0.70)

similar_names <- similar_names %>%
  rename(HS_csf_Name = Name1, MA_csf_Name = Name2)
