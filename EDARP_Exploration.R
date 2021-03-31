## This script explores records validated by biometric keys produced by the EDARP Program. 
## It generates features using demographic data, creates a paired-record dataset, and runs 
## classifiers to learn patterns associated with matched / unmatched records.
##
## Author: Yoni Friedman, Palladium
## Last Edited: March 31, 2021


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
print(sprintf("Of %f total records in the HTS dataset, there are %f unique records",
              nrow(dat), nrow(recs_unique)))


# Explore matches -------------------------------

# Idetify biometrics that appear more than once (these will disagree on at least one other field)
dup_biometrics <- recs_unique %>%
  group_by(BiometricKey) %>%
  summarize(count = n()) %>%
  filter(count > 1) %>%
  ungroup() %>%
  .$BiometricKey

# Create shell to hold biometric keys
shell <- data.frame(BiometricKey = dup_biometrics)

# Loop through each biometric key and identify the number of unique values for each column
for(i in 1:length(dup_biometrics)){
  df_temp <- recs_unique %>% filter(BiometricKey == dup_biometrics[i])
  # for each columns, get number of unique values
  shell[i, 2:15] <- apply(df_temp[, 2:15], 2, n_distinct)
}

# Assign original column names to shell
names(shell)[2:15] <- names(dat)[2:15]

# When biometric keys appear more than once, on what fields do they disagree?
apply(shell, 2, function(x){sum(x>1)})


# Create paired dataset ----------------------
# If we paired all records with each other, we'd produce almost 2 billion record pairs
# Instead, let's take a subset of all record pairs that includes all matches and a random
# subset of records that do not match

# First, identify all records with a match
match_keys <- recs_unique[recs_unique$BiometricKey %in% dup_biometrics, ]

# Create a counter for each biometric key so that we create exact duplicate record pairs
# ie we want the record pairs where the biometric keys match but differ on another variable
match_keys <- match_keys %>%
  group_by(BiometricKey) %>%
  mutate(rownum = row_number())

# Join on itself by biometric keys, and drop key
matches <- merge(match_keys, match_keys, by = "BiometricKey") %>%
  filter(rownum.x != rownum.y) %>%
  select(-BiometricKey, -rownum.x, -rownum.y)

# Name records as A / B
names(matches) <- gsub("\\.x", "A", names(matches))
names(matches) <- gsub("\\.y", "B", names(matches))

# Create Label
matches$target <- 1

## Randomly select negative pairs

# First, get keys that do not match and select 30,000 at random
unmatch_keys <- recs_unique[!recs_unique$BiometricKey %in% dup_biometrics, ]

# Set seed for reproducability
set.seed(2231)
# Sample 30,000 rows
samp_a <- sample(1:nrow(unmatch_keys), 30000, replace = TRUE)
unmatch_keys <- unmatch_keys[samp_a, ]
# Drop biometric key
unmatch_keys <- unmatch_keys %>% select(-BiometricKey)
# Rename variables for A record
names(unmatch_keys) <- paste0(names(unmatch_keys), "A")

# Next, select any 30,000 records at random (even those that might match since they will not match here)
samp_b <- sample(1:nrow(recs_unique), 30000, replace = TRUE)
unmatch_other <- recs_unique[samp_b, ]
# Drop biometric key
unmatch_other <- unmatch_other %>% select(-BiometricKey)
# Rename variables for A record
names(unmatch_other) <- paste0(names(unmatch_other), "B")

# Horizontally stack to create pairs
unmatches <- cbind(unmatch_keys, unmatch_other)

# Create target
unmatches$target <- 0

# Now, vertically stack matches and unmatches
matches_labeled <- rbind(matches, unmatches)

# Generate Facility Features ------------------
# Create flag for whether facility codes match
matches_labeled$facilityMatch <- ifelse(matches_labeled$FacilityCodeA == matches_labeled$FacilityCodeB, 1, 0)
table(matches_labeled$facilityMatch, matches_labeled$target, dnn = c("Facility Match", "Target"))

## EXPLORE Geocode facilities and calculate proximity

# Generate Gender features -----------------

# Generate feature for whether records match on gender
matches_labeled$genderMatch <- ifelse(matches_labeled$GenderA == matches_labeled$GenderB, 1, 0)
table(matches_labeled$genderMatch, matches_labeled$target, dnn = c("Gender Match", "Target"))

# Generate name-based features for Soundex encodings -------------
# Create indicators for whether names match on first letter of first/last name
matches_labeled$flFNSXMatch <- ifelse(substr(matches_labeled$sxFirstNameA, 1, 1) ==
                                        substr(matches_labeled$sxFirstNameB, 1, 1), 1, 0)
matches_labeled$flMNSXMatch <- ifelse(substr(matches_labeled$sxMiddleNameA, 1, 1) ==
                                        substr(matches_labeled$sxMiddleNameB, 1, 1), 1, 0)
matches_labeled$flLNSXMatch <- ifelse(substr(matches_labeled$sxLastNameA, 1, 1) ==
                                        substr(matches_labeled$sxLastNameB, 1, 1), 1, 0)
table(matches_labeled$flFNSXMatch, matches_labeled$target, dnn = c("First Name FL Match", "Target"))
table(matches_labeled$flLNSXMatch, matches_labeled$target, dnn = c("Last Name FL Match", "Target"))

# Calculate string distance metrics
## Use of these on encodings is controversial - would be best to apply directly on names
# Will try Hamming, Levenshtein, Jaro-Winkler, and Longest Common String
matches_labeled$fnSXHamming <- stringsim(matches_labeled$sxFirstNameA, matches_labeled$sxFirstNameB, method = "hamming")
matches_labeled$fnSXLV <- stringsim(matches_labeled$sxFirstNameA, matches_labeled$sxFirstNameB, method = "lv")
matches_labeled$fnSXJW <- stringsim(matches_labeled$sxFirstNameA, matches_labeled$sxFirstNameB, method = "jw")
matches_labeled$fnSXLCS <- stringsim(matches_labeled$sxFirstNameA, matches_labeled$sxFirstNameB, method = "lcs")
matches_labeled$mnSXHamming <- stringsim(matches_labeled$sxMiddleNameA, matches_labeled$sxMiddleNameB, method = "hamming")
matches_labeled$mnSXLV <- stringsim(matches_labeled$sxMiddleNameA, matches_labeled$sxMiddleNameB, method = "lv")
matches_labeled$mnSXJW <- stringsim(matches_labeled$sxMiddleNameA, matches_labeled$sxMiddleNameB, method = "jw")
matches_labeled$mnSXLCS <- stringsim(matches_labeled$sxMiddleNameA, matches_labeled$sxMiddleNameB, method = "lcs")
matches_labeled$lnSXHamming <- stringsim(matches_labeled$sxLastNameA, matches_labeled$sxLastNameB, method = "hamming")
matches_labeled$lnSXLV <- stringsim(matches_labeled$sxLastNameA, matches_labeled$sxLastNameB, method = "lv")
matches_labeled$lnSXJW <- stringsim(matches_labeled$sxLastNameA, matches_labeled$sxLastNameB, method = "jw")
matches_labeled$lnSXLCS <- stringsim(matches_labeled$sxLastNameA, matches_labeled$sxLastNameB, method = "lcs")

# Plot histograms of similarity metrics by match/unmatch
matches_labeled %>%
  ggplot(aes(x=fnSXHamming)) +
  geom_histogram() +
  facet_wrap(~target)

matches_labeled %>%
  ggplot(aes(x=mnSXHamming)) +
  geom_histogram() +
  facet_wrap(~target)

matches_labeled %>%
  ggplot(aes(x=lnSXHamming)) +
  geom_histogram() +
  facet_wrap(~target)

matches_labeled %>%
  ggplot(aes(x=fnSXJW)) +
  geom_histogram() +
  facet_wrap(~target)

matches_labeled %>%
  ggplot(aes(x=mnSXJW)) +
  geom_histogram() +
  facet_wrap(~target)

matches_labeled %>%
  ggplot(aes(x=lnSXJW)) +
  geom_histogram() +
  facet_wrap(~target)

# Generate name-based features for Double Metaphone encodings -------------

# Some names have double encodings. For these let's compare all encodings across records and take the 
# max value

# loop through record pairs
for(i in 1:nrow(matches_labeled)){
  
  # Get record pair
  dat <- matches_labeled[i, ]
  
  # Calculate similarity metrics for first name -----------------
  # Get Record A first name and Split
  rec_a_fn <- unlist(str_split(dat$dmFirstNameA, ";"))
  
  # Get Record B first name and Split
  rec_b_fn <- unlist(str_split(dat$dmFirstNameB, ";"))
  
  # Create grid
  grid <- expand.grid(rec_a_fn = rec_a_fn, rec_b_fn = rec_b_fn, stringsAsFactors = FALSE)
  
  # Calculate Similarity Metrics
  grid <- grid %>% mutate("Hamming" = stringsim(rec_a_fn, rec_b_fn, method = "hamming"),
                          "LV" = stringsim(rec_a_fn, rec_b_fn, method = "lv"),
                          "JW" = stringsim(rec_a_fn, rec_b_fn, method = "jw"),
                          "LCS" = stringsim(rec_a_fn, rec_b_fn, method = "lcs"))
  # Get highest value
  matches_labeled[i, "fnDMHamming"] <- max(grid$Hamming)
  matches_labeled[i, "fnDMLV"] <- max(grid$LV)
  matches_labeled[i, "fnDMJW"] <- max(grid$JW)
  matches_labeled[i, "fnDMLCS"] <- max(grid$LCS)
  
  # Repeate for Middle Name -----------------
  # Get Record A middle name and Split
  rec_a_mn <- unlist(str_split(dat$dmMiddleNameA, ";"))
  
  # Get Record B middle name and Split
  rec_b_mn <- unlist(str_split(dat$dmMiddleNameB, ";"))
  
  # Create grid
  grid <- expand.grid(rec_a_mn = rec_a_mn, rec_b_mn = rec_b_mn, stringsAsFactors = FALSE)
  
  # Calculate Similarity Metrics
  grid <- grid %>% mutate("Hamming" = stringsim(rec_a_mn, rec_b_mn, method = "hamming"),
                          "LV" = stringsim(rec_a_mn, rec_b_mn, method = "lv"),
                          "JW" = stringsim(rec_a_mn, rec_b_mn, method = "jw"),
                          "LCS" = stringsim(rec_a_mn, rec_b_mn, method = "lcs"))
  # Get highest value
  matches_labeled[i, "mnDMHamming"] <- max(grid$Hamming)
  matches_labeled[i, "mnDMLV"] <- max(grid$LV)
  matches_labeled[i, "mnDMJW"] <- max(grid$JW)
  matches_labeled[i, "mnDMLCS"] <- max(grid$LCS)
  
  # Repeat for Last Name --------------------
  # Get Record A first name and Split
  rec_a_ln <- unlist(str_split(dat$dmLastNameA, ";"))
  
  # Get Record B first name and Split
  rec_b_ln <- unlist(str_split(dat$dmLastNameB, ";"))
  
  # Create grid
  grid <- expand.grid(rec_a_fn = rec_a_ln, rec_b_ln = rec_b_ln, stringsAsFactors = FALSE)
  
  # Calculate Similarity Metrics
  grid <- grid %>% mutate("Hamming" = stringsim(rec_a_ln, rec_b_ln, method = "hamming"),
                          "LV" = stringsim(rec_a_ln, rec_b_ln, method = "lv"),
                          "JW" = stringsim(rec_a_ln, rec_b_ln, method = "jw"),
                          "LCS" = stringsim(rec_a_ln, rec_b_ln, method = "lcs"))
  # Get highest value
  matches_labeled[i, "lnDMHamming"] <- max(grid$Hamming)
  matches_labeled[i, "lnDMLV"] <- max(grid$LV)
  matches_labeled[i, "lnDMJW"] <- max(grid$JW)
  matches_labeled[i, "lnDMLCS"] <- max(grid$LCS)
}

# Plot histograms of similarity metrics by match/unmatch
matches_labeled %>%
  ggplot(aes(x=fnDMHamming)) +
  geom_histogram() +
  facet_wrap(~target)

matches_labeled %>%
  ggplot(aes(x=mnDMHamming)) +
  geom_histogram() +
  facet_wrap(~target)

matches_labeled %>%
  ggplot(aes(x=lnDMHamming)) +
  geom_histogram() +
  facet_wrap(~target)

matches_labeled %>%
  ggplot(aes(x=fnDMJW)) +
  geom_histogram() +
  facet_wrap(~target)

matches_labeled %>%
  ggplot(aes(x=mnDMJW)) +
  geom_histogram() +
  facet_wrap(~target)

matches_labeled %>%
  ggplot(aes(x=lnDMJW)) +
  geom_histogram() +
  facet_wrap(~target)

# Create DOB-based features ------------------

# Get year, month, day of birth
matches_labeled$YOBA <- substr(matches_labeled$dmPKValueDoBA,
                               nchar(matches_labeled$dmPKValueDoBA) - 7,
                               nchar(matches_labeled$dmPKValueDoBA) - 4)
matches_labeled$MOBA <- substr(matches_labeled$dmPKValueDoBA,
                               nchar(matches_labeled$dmPKValueDoBA) - 3,
                               nchar(matches_labeled$dmPKValueDoBA) - 2)
matches_labeled$DOBA <- substr(matches_labeled$dmPKValueDoBA,
                               nchar(matches_labeled$dmPKValueDoBA) - 1,
                               nchar(matches_labeled$dmPKValueDoBA))
matches_labeled$YOBB <- substr(matches_labeled$dmPKValueDoBB,
                               nchar(matches_labeled$dmPKValueDoBB) - 7,
                               nchar(matches_labeled$dmPKValueDoBB) - 4)
matches_labeled$MOBB <- substr(matches_labeled$dmPKValueDoBB,
                               nchar(matches_labeled$dmPKValueDoBB) - 3,
                               nchar(matches_labeled$dmPKValueDoBB) - 2)
matches_labeled$DOBB <- substr(matches_labeled$dmPKValueDoBB,
                               nchar(matches_labeled$dmPKValueDoBB) - 1,
                               nchar(matches_labeled$dmPKValueDoBB))

# For each date element, get string similarities (will be 1 if match)
matches_labeled$YOBHamming <- stringsim(matches_labeled$YOBA, matches_labeled$YOBB, method = "hamming")
matches_labeled$YOBLV <- stringsim(matches_labeled$YOBA, matches_labeled$YOBB, method = "lv")
matches_labeled$YOBJW <- stringsim(matches_labeled$YOBA, matches_labeled$YOBB, method = "jw")
matches_labeled$YOBLCS <- stringsim(matches_labeled$YOBA, matches_labeled$YOBB, method = "lcs")
matches_labeled$MOBHamming <- stringsim(matches_labeled$MOBA, matches_labeled$MOBB, method = "hamming")
matches_labeled$MOBLV <- stringsim(matches_labeled$MOBA, matches_labeled$MOBB, method = "lv")
matches_labeled$MOBJW <- stringsim(matches_labeled$MOBA, matches_labeled$MOBB, method = "jw")
matches_labeled$MOBLCS <- stringsim(matches_labeled$MOBA, matches_labeled$MOBB, method = "lcs")
matches_labeled$DOBHamming <- stringsim(matches_labeled$DOBA, matches_labeled$DOBB, method = "hamming")
matches_labeled$DOBLV <- stringsim(matches_labeled$DOBA, matches_labeled$DOBB, method = "lv")
matches_labeled$DOBJW <- stringsim(matches_labeled$DOBA, matches_labeled$DOBB, method = "jw")
matches_labeled$DOBLCS <- stringsim(matches_labeled$DOBA, matches_labeled$DOBB, method = "lcs")

## EXPLORE other date comparisons, month/day switches

# Save Out ---------------------

matches_out <- matches_labeled %>%
  select(match("target", names(matches_labeled)):match("lnDMLCS", names(matches_labeled)), 
         match("YOBHamming", names(matches_labeled)):ncol(.))

saveRDS(matches_out, './matches_prep.rds')
