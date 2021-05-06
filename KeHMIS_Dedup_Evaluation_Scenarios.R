## This script runs classifiers to learn patterns associated with matched / unmatched records.
##
## Author: Yoni Friedman, Palladium
## Last Edited: April 11, 2021

# Load libraries and data ----------------------
# Set Working Directory
setwd("~/Kenya/Deduplication")

# Load Libraries
library(dplyr); library(stringr); library(lubridate); library(stringdist); library(ggplot2)

# Load Data
dat <- read.csv('./Data/7_Master_Patient_Indices.csv', stringsAsFactors = FALSE)

# Identify complete duplicates ------------------------
# True duplicates agree on all fields, including biometric key and demogrpahic fields
# These are not analytically useful, and we will drop
recs_unique <- unique(dat)

# Compute number of record pairs -------------------
num_records <- as.numeric(nrow(recs_unique))

num_record_pairs <- ((num_records*num_records)-num_records)/2

# Compute number of positives --------------------
dup_biometrics <- recs_unique %>%
  group_by(BiometricKey) %>%
  summarize(count = n()) %>%
  filter(count > 1) %>%
  ungroup() %>%
  .$BiometricKey

# First, identify all records with a match
match_keys <- recs_unique[recs_unique$BiometricKey %in% dup_biometrics, ]

# Create a counter for each biometric key so that we create exact duplicate record pairs
# ie we want the record pairs where the biometric keys match but differ on another variable
match_keys <- match_keys %>%
  group_by(BiometricKey) %>%
  mutate(rownum = as.numeric(row_number() + runif(1)))

# Join on itself by biometric keys, and drop key
matches <- merge(match_keys, match_keys, by = "BiometricKey") %>%
  filter(rownum.x != rownum.y) %>%
  mutate(rownumsum = rownum.x + rownum.y) %>%
  group_by(BiometricKey, rownumsum) %>%
  mutate(rownum = row_number()) %>%
  filter(rownum == 1)

num_positives <- nrow(matches)
num_negatives <- num_record_pairs - num_positives

print(num_record_pairs)
print(num_positives)
print(num_negatives)

# Prep dataset ----------------------------------
# Prep YOB, MOB and Gender
recs_unique$YOB <- substr(recs_unique$dmPKValueDoB,
                          nchar(recs_unique$dmPKValueDoB) - 7,
                          nchar(recs_unique$dmPKValueDoB) - 4)

recs_unique$MOB <- substr(recs_unique$dmPKValueDoB,
                          nchar(recs_unique$dmPKValueDoB) - 3,
                          nchar(recs_unique$dmPKValueDoB) - 2)

recs_unique$Gender <- ifelse(recs_unique$Gender == 1, "M", "F")

# For double metaphone, use last three characters
recs_unique$dmLastName_AfterSC <- gsub(".*;", "", recs_unique$dmLastName)

# Implement KeHMIS Algorithm  --------------------------
recs_split <- split(recs_unique, list(recs_unique$Gender, recs_unique$YOB), drop = TRUE)

# Calculate PPV and Sensitivity for current PKV, adding middle name, and dropping MOB, for different thresholds
# For each group, create full index
thresholds <- c(.94, .95, .96, .97, .98)
performance <- data.frame()
performance_mn <- data.frame()
performance_mob <- data.frame()
performance_sx <- data.frame()
performance_mnmob <- data.frame()

for(i in 1:length(recs_split)){
  print(i)
  dat <- recs_split[[i]]
  if(nrow(dat) == 1){next}
  dat$rownum <- as.numeric(row.names(dat)) + runif(nrow(dat))
  recs_match <- merge(dat, dat, by = c("Gender", "YOB"))
  
  # Drop duplicate 
  recs_match <- recs_match %>% filter(rownum.x != rownum.y)
  recs_match <- recs_match %>%
    mutate(rownumsum = rownum.x + rownum.y) %>%
    group_by(rownumsum) %>%
    mutate(rownum = row_number()) %>%
    filter(rownum == 1)
  
  # Prep label
  recs_match$label <- ifelse(recs_match$BiometricKey.x == recs_match$BiometricKey.y, 1, 0)
  
  # Generate PKVs
  recs_match$PKVX <- paste0(recs_match$Gender, recs_match$sxFirstName.x, recs_match$dmLastName_AfterSC.x,
                            recs_match$YOB, recs_match$MOB.x)
  recs_match$PKVX_MN <- paste0(recs_match$Gender, recs_match$sxFirstName.x, recs_match$sxMiddleName.x,
                               recs_match$dmLastName_AfterSC.x, recs_match$YOB, recs_match$MOB.x)
  recs_match$PKVX_MOB <- paste0(recs_match$Gender, recs_match$sxFirstName.x, recs_match$dmLastName_AfterSC.x, recs_match$YOB)
  recs_match$PKVX_SX <- paste0(recs_match$Gender, recs_match$sxFirstName.x, recs_match$sxLastName.x,
                               recs_match$YOB, recs_match$MOB.x)
  recs_match$PKVX_MNMOB <- paste0(recs_match$Gender, recs_match$sxFirstName.x, recs_match$sxMiddleName.x,
                            recs_match$dmLastName_AfterSC.x, recs_match$YOB)
  
  recs_match$PKVY <- paste0(recs_match$Gender, recs_match$sxFirstName.y, recs_match$dmLastName_AfterSC.y,
                            recs_match$YOB, recs_match$MOB.y)
  recs_match$PKVY_MN <- paste0(recs_match$Gender, recs_match$sxFirstName.y, recs_match$sxMiddleName.y,
                               recs_match$dmLastName_AfterSC.y, recs_match$YOB, recs_match$MOB.y)
  recs_match$PKVY_MOB <- paste0(recs_match$Gender, recs_match$sxFirstName.y, recs_match$dmLastName_AfterSC.y, recs_match$YOB)
  recs_match$PKVY_SX <- paste0(recs_match$Gender, recs_match$sxFirstName.y, recs_match$sxLastName.y,
                               recs_match$YOB, recs_match$MOB.y)
  recs_match$PKVY_MNMOB <- paste0(recs_match$Gender, recs_match$sxFirstName.y, recs_match$sxMiddleName.y,
                                 recs_match$dmLastName_AfterSC.y, recs_match$YOB)
  
  # Compare PKVs
  recs_match$jw <- stringsim(recs_match$PKVX, recs_match$PKVY, method = "jw")
  recs_match$jw_mn <- stringsim(recs_match$PKVX_MN, recs_match$PKVY_MN, method = "jw")
  recs_match$jw_mob <- stringsim(recs_match$PKVX_MOB, recs_match$PKVY_MOB, method = "jw")
  recs_match$jw_sx <- stringsim(recs_match$PKVX_SX, recs_match$PKVY_SX, method = "jw")
  recs_match$jw_mnmob <- stringsim(recs_match$PKVX_MNMOB, recs_match$PKVY_MNMOB, method = "jw")
  
  # Loop through different thresholds
  for(j in 1:length(thresholds)){
    
    recs_match$pred <- ifelse(recs_match$jw >= thresholds[j], 1, 0)
    
    TP <- nrow(filter(recs_match, pred == 1 & label == 1))
    FP <- nrow(filter(recs_match, pred == 1 & label == 0))
    TN <- nrow(filter(recs_match, pred == 0 & label == 0))
    FN <- nrow(filter(recs_match, pred == 0 & label == 1))
    Num_Preds <- nrow(recs_match)
    
    perf_tmp <- data.frame(TP = TP,
                           FP = FP,
                           TN = TN,
                           FN = FN,
                           Num_Preds = Num_Preds,
                           Threshold = thresholds[j],
                           Formula = "Current")
    performance <- rbind(performance, perf_tmp)
    
  }
  
  # Loop through different thresholds for formulation with middle name
  for(j in 1:length(thresholds)){
    
    recs_match$pred <- ifelse(recs_match$jw_mn >= thresholds[j], 1, 0)
    
    TP <- nrow(filter(recs_match, pred == 1 & label == 1))
    FP <- nrow(filter(recs_match, pred == 1 & label == 0))
    TN <- nrow(filter(recs_match, pred == 0 & label == 0))
    FN <- nrow(filter(recs_match, pred == 0 & label == 1))
    Num_Preds <- nrow(recs_match)
    
    perf_tmp <- data.frame(TP = TP,
                           FP = FP,
                           TN = TN,
                           FN = FN,
                           Num_Preds = Num_Preds,
                           Threshold = thresholds[j],
                           Formula = "MN")
    performance_mn <- rbind(performance_mn, perf_tmp)
    
  }
  
  # Loop through different thresholds for formulation without month of birth
  for(j in 1:length(thresholds)){
    
    recs_match$pred <- ifelse(recs_match$jw_mob >= thresholds[j], 1, 0)
    
    TP <- nrow(filter(recs_match, pred == 1 & label == 1))
    FP <- nrow(filter(recs_match, pred == 1 & label == 0))
    TN <- nrow(filter(recs_match, pred == 0 & label == 0))
    FN <- nrow(filter(recs_match, pred == 0 & label == 1))
    Num_Preds <- nrow(recs_match)
    
    perf_tmp <- data.frame(TP = TP,
                           FP = FP,
                           TN = TN,
                           FN = FN,
                           Num_Preds = Num_Preds,
                           Threshold = thresholds[j],
                           Formula = "MOB")
    performance_mob <- rbind(performance_mob, perf_tmp)
    
  }
  
  # Loop through different thresholds for formulation without month of birth
  for(j in 1:length(thresholds)){
    
    recs_match$pred <- ifelse(recs_match$jw_sx >= thresholds[j], 1, 0)
    
    TP <- nrow(filter(recs_match, pred == 1 & label == 1))
    FP <- nrow(filter(recs_match, pred == 1 & label == 0))
    TN <- nrow(filter(recs_match, pred == 0 & label == 0))
    FN <- nrow(filter(recs_match, pred == 0 & label == 1))
    Num_Preds <- nrow(recs_match)
    
    perf_tmp <- data.frame(TP = TP,
                           FP = FP,
                           TN = TN,
                           FN = FN,
                           Num_Preds = Num_Preds,
                           Threshold = thresholds[j],
                           Formula = "SX")
    performance_sx <- rbind(performance_sx, perf_tmp)
    
  }
  
  # Loop through different thresholds for formulation with middle name and without month of birth
  for(j in 1:length(thresholds)){
    
    recs_match$pred <- ifelse(recs_match$jw_mnmob >= thresholds[j], 1, 0)
    
    TP <- nrow(filter(recs_match, pred == 1 & label == 1))
    FP <- nrow(filter(recs_match, pred == 1 & label == 0))
    TN <- nrow(filter(recs_match, pred == 0 & label == 0))
    FN <- nrow(filter(recs_match, pred == 0 & label == 1))
    Num_Preds <- nrow(recs_match)
    
    perf_tmp <- data.frame(TP = TP,
                           FP = FP,
                           TN = TN,
                           FN = FN,
                           Num_Preds = Num_Preds,
                           Threshold = thresholds[j],
                           Formula = "MN_MOB")
    performance_sx <- rbind(performance_sx, perf_tmp)
    
  }
  
  
}

performance_summary <- performance %>%
  rbind(., performance_mn) %>%
  rbind(., performance_mob) %>%
  rbind(., performance_sx) %>%
  rbind(., performance_mnmob) %>%
  group_by(Formula, Threshold) %>%
  summarize(PPV = sum(TP, na.rm = TRUE) / (sum(TP, na.rm = TRUE) + sum(FP, na.rm = TRUE)),
            Sensitivity = sum(TP, na.rm = TRUE) / num_positives) %>%
  mutate(F1 = (2 * PPV * Sensitivity) / (PPV + Sensitivity))

print(performance_summary)

# Implement KeHMIS Algorithm by Facility --------------------------
recs_split <- split(recs_unique, list(recs_unique$Gender, recs_unique$YOB, recs_unique$FacilityCode), drop = TRUE)


# Calculate PPV and Sensitivity for current PKV, adding middle name, and dropping MOB, for different thresholds
# For each group, create full index
thresholds <- c(.94, .95, .96, .97, .98)
performance <- data.frame()
performance_mn <- data.frame()
performance_mob <- data.frame()
performance_sx <- data.frame()
performance_mnmob <- data.frame()

for(i in 1:length(recs_split)){
  print(i)
  dat <- recs_split[[i]]
  if(nrow(dat) == 1){next}
  dat$rownum <- as.numeric(row.names(dat)) + runif(nrow(dat))
  recs_match <- merge(dat, dat, by = c("Gender", "YOB"))
  
  # Drop duplicate 
  recs_match <- recs_match %>% filter(rownum.x != rownum.y)
  recs_match <- recs_match %>%
    mutate(rownumsum = rownum.x + rownum.y) %>%
    group_by(rownumsum) %>%
    mutate(rownum = row_number()) %>%
    filter(rownum == 1)
  
  # Prep label
  recs_match$label <- ifelse(recs_match$BiometricKey.x == recs_match$BiometricKey.y, 1, 0)
  
  # Generate PKVs
  recs_match$PKVX <- paste0(recs_match$Gender, recs_match$sxFirstName.x, recs_match$dmLastName_AfterSC.x,
                            recs_match$YOB, recs_match$MOB.x)
  recs_match$PKVX_MN <- paste0(recs_match$Gender, recs_match$sxFirstName.x, recs_match$sxMiddleName.x,
                               recs_match$dmLastName_AfterSC.x, recs_match$YOB, recs_match$MOB.x)
  recs_match$PKVX_MOB <- paste0(recs_match$Gender, recs_match$sxFirstName.x, recs_match$dmLastName_AfterSC.x, recs_match$YOB)
  recs_match$PKVX_SX <- paste0(recs_match$Gender, recs_match$sxFirstName.x, recs_match$sxLastName.x,
                               recs_match$YOB, recs_match$MOB.x)
  recs_match$PKVX_MNMOB <- paste0(recs_match$Gender, recs_match$sxFirstName.x, recs_match$sxMiddleName.x,
                                  recs_match$dmLastName_AfterSC.x, recs_match$YOB)
  
  recs_match$PKVY <- paste0(recs_match$Gender, recs_match$sxFirstName.y, recs_match$dmLastName_AfterSC.y,
                            recs_match$YOB, recs_match$MOB.y)
  recs_match$PKVY_MN <- paste0(recs_match$Gender, recs_match$sxFirstName.y, recs_match$sxMiddleName.y,
                               recs_match$dmLastName_AfterSC.y, recs_match$YOB, recs_match$MOB.y)
  recs_match$PKVY_MOB <- paste0(recs_match$Gender, recs_match$sxFirstName.y, recs_match$dmLastName_AfterSC.y, recs_match$YOB)
  recs_match$PKVY_SX <- paste0(recs_match$Gender, recs_match$sxFirstName.y, recs_match$sxLastName.y,
                               recs_match$YOB, recs_match$MOB.y)
  recs_match$PKVY_MNMOB <- paste0(recs_match$Gender, recs_match$sxFirstName.y, recs_match$sxMiddleName.y,
                                  recs_match$dmLastName_AfterSC.y, recs_match$YOB)
  
  # Compare PKVs
  recs_match$jw <- stringsim(recs_match$PKVX, recs_match$PKVY, method = "jw")
  recs_match$jw_mn <- stringsim(recs_match$PKVX_MN, recs_match$PKVY_MN, method = "jw")
  recs_match$jw_mob <- stringsim(recs_match$PKVX_MOB, recs_match$PKVY_MOB, method = "jw")
  recs_match$jw_sx <- stringsim(recs_match$PKVX_SX, recs_match$PKVY_SX, method = "jw")
  recs_match$jw_mnmob <- stringsim(recs_match$PKVX_MNMOB, recs_match$PKVY_MNMOB, method = "jw")
  
  # Loop through different thresholds
  for(j in 1:length(thresholds)){
    
    recs_match$pred <- ifelse(recs_match$jw >= thresholds[j], 1, 0)
    
    TP <- nrow(filter(recs_match, pred == 1 & label == 1))
    FP <- nrow(filter(recs_match, pred == 1 & label == 0))
    TN <- nrow(filter(recs_match, pred == 0 & label == 0))
    FN <- nrow(filter(recs_match, pred == 0 & label == 1))
    Num_Preds <- nrow(recs_match)
    
    perf_tmp <- data.frame(TP = TP,
                           FP = FP,
                           TN = TN,
                           FN = FN,
                           Num_Preds = Num_Preds,
                           Threshold = thresholds[j],
                           Formula = "Current")
    performance <- rbind(performance, perf_tmp)
    
  }
  
  # Loop through different thresholds for formulation with middle name
  for(j in 1:length(thresholds)){
    
    recs_match$pred <- ifelse(recs_match$jw_mn >= thresholds[j], 1, 0)
    
    TP <- nrow(filter(recs_match, pred == 1 & label == 1))
    FP <- nrow(filter(recs_match, pred == 1 & label == 0))
    TN <- nrow(filter(recs_match, pred == 0 & label == 0))
    FN <- nrow(filter(recs_match, pred == 0 & label == 1))
    Num_Preds <- nrow(recs_match)
    
    perf_tmp <- data.frame(TP = TP,
                           FP = FP,
                           TN = TN,
                           FN = FN,
                           Num_Preds = Num_Preds,
                           Threshold = thresholds[j],
                           Formula = "MN")
    performance_mn <- rbind(performance_mn, perf_tmp)
    
  }
  
  # Loop through different thresholds for formulation without month of birth
  for(j in 1:length(thresholds)){
    
    recs_match$pred <- ifelse(recs_match$jw_mob >= thresholds[j], 1, 0)
    
    TP <- nrow(filter(recs_match, pred == 1 & label == 1))
    FP <- nrow(filter(recs_match, pred == 1 & label == 0))
    TN <- nrow(filter(recs_match, pred == 0 & label == 0))
    FN <- nrow(filter(recs_match, pred == 0 & label == 1))
    Num_Preds <- nrow(recs_match)
    
    perf_tmp <- data.frame(TP = TP,
                           FP = FP,
                           TN = TN,
                           FN = FN,
                           Num_Preds = Num_Preds,
                           Threshold = thresholds[j],
                           Formula = "MOB")
    performance_mob <- rbind(performance_mob, perf_tmp)
    
  }
  
  # Loop through different thresholds for formulation without month of birth
  for(j in 1:length(thresholds)){
    
    recs_match$pred <- ifelse(recs_match$jw_sx >= thresholds[j], 1, 0)
    
    TP <- nrow(filter(recs_match, pred == 1 & label == 1))
    FP <- nrow(filter(recs_match, pred == 1 & label == 0))
    TN <- nrow(filter(recs_match, pred == 0 & label == 0))
    FN <- nrow(filter(recs_match, pred == 0 & label == 1))
    Num_Preds <- nrow(recs_match)
    
    perf_tmp <- data.frame(TP = TP,
                           FP = FP,
                           TN = TN,
                           FN = FN,
                           Num_Preds = Num_Preds,
                           Threshold = thresholds[j],
                           Formula = "SX")
    performance_sx <- rbind(performance_sx, perf_tmp)
    
  }
  
  # Loop through different thresholds for formulation with middle name and without month of birth
  for(j in 1:length(thresholds)){
    
    recs_match$pred <- ifelse(recs_match$jw_mnmob >= thresholds[j], 1, 0)
    
    TP <- nrow(filter(recs_match, pred == 1 & label == 1))
    FP <- nrow(filter(recs_match, pred == 1 & label == 0))
    TN <- nrow(filter(recs_match, pred == 0 & label == 0))
    FN <- nrow(filter(recs_match, pred == 0 & label == 1))
    Num_Preds <- nrow(recs_match)
    
    perf_tmp <- data.frame(TP = TP,
                           FP = FP,
                           TN = TN,
                           FN = FN,
                           Num_Preds = Num_Preds,
                           Threshold = thresholds[j],
                           Formula = "MN_MOB")
    performance_sx <- rbind(performance_sx, perf_tmp)
    
  }
  
  
}

performance_summary <- performance %>%
  rbind(., performance_mn) %>%
  rbind(., performance_mob) %>%
  rbind(., performance_sx) %>%
  rbind(., performance_mnmob) %>%
  group_by(Formula, Threshold) %>%
  summarize(PPV = sum(TP, na.rm = TRUE) / (sum(TP, na.rm = TRUE) + sum(FP, na.rm = TRUE)),
            Sensitivity = sum(TP, na.rm = TRUE) / num_positives) %>%
  mutate(F1 = (2 * PPV * Sensitivity) / (PPV + Sensitivity))

print(performance_summary)







