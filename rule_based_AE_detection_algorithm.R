################################################################
# Title: Rule-based Adverse Event Detection Algorithm        
#
# Inputs:
# 1. Ontology file (txt)  -- list of AE terms/phrases
# 2. Notes file (txt)     -- clinical notes
#
# Outputs:
#   - note_level.csv  : numbers of AEs per note
#   - patient_level.csv : AE present per patient
#
# Usage
# Rscript detect_AEs.R --ontology ontology.txt --notes notes.txt --output results.csv
#
# Dependencies 
# R >= 4.2.0
# Packages:optparse, lubridate, tidytext, word2vec, stringr, quanteda, tokenizers, data.table, dplyr
################################################################

setwd("R:/WorkArea-jsg9018")
suppressPackageStartupMessages({
  library(optparse)
  library(readxl)
  library(lubridate)
  library(tidytext)
  library(word2vec)
  library(stringr)
  library(ggplot2)
  library(quanteda)
  library(quanteda.textstats)
  library(quanteda.textplots)
  library(tokenizers)
  library(data.table)
  library(dplyr)
})
options(width=500)

option_list <- list(
  make_option("--ontology", type = "character", help = "Ontology txt file"),
  make_option("--notes", type = "character", help = "Notes txt file"),
  make_option("--output_prefix", type = "character", help = "Prefix for output files")
)

opt <- parge_args(OptionParser(option_list = option_list))

if(any(sapply(opt, is.null))) {
  stop("Error: Missing arguments.")
}

# Load data
ont <- read_lines(opt$ontology)
df <- read_lines(opt$notes)

# Basic checks
if (length(ontology) == 0) stop("Ontology file is empty")
if (length(notes) == 0) stop("Notes file is empty")
if (!all(c("pat_id", "note_id", "note_text") %in% colnames(notes))) {
  stop("Notes file must have columns: pat_id, note_id, note_text")
} 

# Read in ADE ontology
ont <- read.csv("ADE_ontology.txt", sep="\t", header = T, na="")

# Read in notes and select just pat_id and note_text
#df <- read_excel("ADE_NewNotes_18Mar24_jsgannotated.xlsx")
#df <- df %>% dplyr::select(pat_id, note_text, note_id) %>% filter(!is.na(pat_id))

df <- df %>% mutate(note_id = paste0(pat_id, "_", note_id))

# Clean ADE list to create ADE dictionary
ontdf <- as.data.frame(ont)
ontlist <- as.list(ont)
ontlist <- lapply(ontlist, tolower)
ontlist <- lapply(ontlist, gsub, pattern = " ", replacement= "_")
ontlist <- lapply(ontlist, gsub, pattern = "NA", replacement= NA)
ontlist <- lapply(ontlist, function(x) x[!is.na(x)])
AE_dict <- dictionary(ontlist)

# Tokenize text and identify AEs
corp <- corpus(df, text_field = "note_text", docid_field = "note_id")
tok <- tokens_tolower(tokens_ngrams(tokens(corp, remove_punct = TRUE, remove_numbers = TRUE), n =1:5))

# Establish negation phrases and negation conjugates
negation_keywords <- c("not", "no", "without", "denies", "none",  "negative",
                       "ppx", "prophylaxis", "recommend", "allergies", "indication", "as needed", "reviewed", "recommend", "vaccination", "vaccine")
negative_conjugates <- c("but", "however", "though", "although", "except","aside", "yet", "does have")

# Function for detecting negated AEs
is_negated <- function(tokens, negation_keywords, negative_conjugates, window=7){
  negated <- rep(FALSE, length(tokens))
  for( i in seq_along(tokens)) {
    if(tokens[i] %in% unlist(AE_dict)){
      start <- max(1, i - window)
      end <- min(length(tokens), i + window)
      context <- tokens[start:end]
      negation_indices <- which(context %in% negation_keywords)
      if (length(negation_indices) > 0) {
        for (neg_idx in negation_indices) {
          if (!any(context[neg_idx:i] %in% negative_conjugates) && neg_idx < i) {
            negated[i] <- TRUE
          }
        }
      }
    }
  }
  return(negated)
}

# Generate results by detecting AEs and negated AEs in each note
AE_detection <- lapply(seq_along(tok), function(doc_idx){
  token_list <- tok[[doc_idx]]
  document_result <- list()
  
  for(ae in names(AE_dict)) {
    ae_terms <- AE_dict[[ae]]
    ae_indices <- which(token_list %in% ae_terms)
    negated_flags <- is_negated(token_list, negation_keywords, negative_conjugates)
    
    total_events <- length(ae_indices)
    negated_events <- sum(negated_flags[ae_indices])
    non_negated_events <- total_events - negated_events
    
    document_result[[ae]] <- list(
      total_events = total_events,
      negated_events = negated_events,
      non_negated_events = non_negated_events
    )
  }
  return(document_result)
})

# Turn results into dataframe
doc_ids <- docnames(tok)
df_results <- do.call(rbind, lapply(seq_along(AE_detection), function(i) {
  doc_result <- results[[i]]
  doc_data <- c(doc_id = doc_ids[i])
  for (ae in names(doc_result)) {
    doc_data[paste0(ae)] <- doc_result[[ae]]$non_negated_events
    doc_data[paste0(ae, "_negated")] <- doc_result[[ae]]$negated_events
  }
  pb$tick()
  return(doc_data)
}))

df_results <- as.data.frame(df_results, stringsAsFactors = FALSE)
for (col in names(df_results)[-1]) {
  df_results[[col]] <- as.numeric(df_results[[col]])
}
note_level_results <- df_results %>% tidyr::separate(doc_id, into=c("pat_id", "note_id"), sep = "_", remove = TRUE)


note_level_results <- note_level_results %>% mutate(across(ends_with("_negated"), as.integer)) %>% 
  mutate(across(.cols = ends_with("_negated"),
                .fns = ~ {ae_name <- str_remove(cur_column(), "_negated")
                get(ae_name, note_level_results)  - .x
                },
                .names = "{str_remove(.col, '_negated')}")) %>%
  select(-ends_with("_negated")) %>% mutate(across(where(is.numeric), ~ pmax(., 0)))

write_csv(note_level_results, paste0(opt$output_prefix, "_note_level.csv"))

## Calculate AE detection  at a patient level
AE_cols <- names(AE_dict)
pat_level_results <- note_level_results %>% group_by(pat_id) %>% summarize(across(all_of(AE_cols), ~ as.integer(any(. == 1)), .names = "{.col}")) %>% ungroup()

write_csv(pat_level_results, file = paste0(opt$output_prefix, "patient_level"))

