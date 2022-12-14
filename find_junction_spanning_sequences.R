#!/usr/bin/env Rscript

library(optparse)

option_list <- list(
  make_option(c("--id"), dest = "id",
              help = "Identifier for the current FASTQ file (default: FASTQ file name"),

  make_option(c("--fastq"), dest = "fastq_file",
              help = "FASTQ file containing sequences to be searched for breakpoint junctions"),

  make_option(c("--flanking-sequences"), dest = "flanking_sequence_file", default = "flanking_sequences.csv",
              help = "CSV file containing flanking sequences on either side of the breakpoint junction; this is expected to contain ID, LeftFlankingSequence, RightFlankingSequence and JunctionSequence columns (default: flanking_sequences.csv)"),

  make_option(c("--output"), dest = "output_file", default = "junction_sequence_matches.txt",
              help = "Tab-delimited output file containing sequences that match breakpoint junctions (default: junction_sequence_matches.txt)"),

  make_option(c("--max-distance"), dest = "max_distance", type = "integer", default = 2,
              help = "Maximum edit distance for matches to each of the flanking sequences (default: 2)"),

  make_option(c("--umi-length"), dest = "umi_length", type = "integer", default = 0,
              help = "Number of UMI bases at beginning of read to omit from search (default: 0)")
)

option_parser <- OptionParser(usage = "usage: %prog [options]", option_list = option_list, add_help_option = TRUE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) args <- "--help"

opt <- parse_args(option_parser, args)

id <- opt$id
fastq_file <- opt$fastq_file
flanking_sequence_file <- opt$flanking_sequence_file
output_file <- opt$output_file
max_distance <- opt$max_distance
umi_length <- opt$umi_length

if (is.null(fastq_file)) stop("FASTQ file must be specified")
if (is.null(id)) id <- fastq_file

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(stringr)
  library(stringi)
})


# reverse complement helper function
# note stri_reverse function from stringi package, a stringr dependency
reverse_complement <- function(sequence) {
  stri_reverse(chartr("ACGTacgt", "TGCAtgca", sequence))
}


# function for finding a given query sequence with a set of sequences within a given edit distance
find_matches <- function(query_id, query_sequence, ids, sequences, max_distance = 0) {
  matches <- aregexec(query_sequence, sequences, max_distance, fixed = FALSE)
  tibble(
    query_id,
    query_sequence = query_sequence,
    id = ids,
    sequence = sequences,
    position = unlist(matches)
  ) %>%
    mutate(match = as.character(regmatches(sequences, matches))) %>%
    filter(position != -1)
}


# read the flanking sequences file
flanking_sequences <- read_csv(
  flanking_sequence_file,
  col_types = cols(
    ID = col_character(),
    LeftFlankingSequence = col_character(),
    RightFlankingSequence = col_character(),
    .default = col_character()
  )
)

# check the expected columns exist
expected_columns <- c("ID", "LeftFlankingSequence", "RightFlankingSequence")
missing_columns <- setdiff(expected_columns, colnames(flanking_sequences))
if (length(missing_columns) > 0) {
  stop("missing columns found in ", flanking_sequence_file, ": '", str_c(missing_columns, collapse = "', '"), "'")
}

# check there are no missing values for the expected columns
if (nrow(filter(flanking_sequences, is.na(ID))) > 0) {
  stop("missing IDs found in ", flanking_sequence_file)
}
if (nrow(filter(flanking_sequences, is.na(LeftFlankingSequence))) > 0 ||
    nrow(filter(flanking_sequences, is.na(RightFlankingSequence))) > 0) {
  stop("missing flanking sequences found in ", flanking_sequence_file)
}

# check for uniqueness of the ID column
duplicates <- flanking_sequences %>%
  count(ID) %>%
  filter(n > 1) %>%
  pull(ID)
if (length(duplicates) > 0) {
  stop("duplicate IDs found in ", flanking_sequence_file, ": '", str_c(duplicates, collapse = "', '"), "'")
}

# add optional JunctionSequence column if it doesn't exist
if (!"JunctionSequence" %in% colnames(flanking_sequences)) {
  flanking_sequences <- mutate(flanking_sequences, JunctionSequence = as.character(NA))
}

flanking_sequences <- select(flanking_sequences, ID, LeftFlankingSequence, RightFlankingSequence, JunctionSequence)

# get reverse complement sequences
flanking_sequences <- flanking_sequences %>%
  mutate(LeftFlankingSequenceReverseComplement = reverse_complement(LeftFlankingSequence)) %>%
  mutate(RightFlankingSequenceReverseComplement = reverse_complement(RightFlankingSequence)) %>%
  mutate(JunctionSequenceReverseComplement = reverse_complement(JunctionSequence))

query_sequences <- bind_rows(
  transmute(flanking_sequences, ID, Direction = "Forward",  LeftFlankingSequence, RightFlankingSequence, LeftQuerySequence = LeftFlankingSequence, RightQuerySequence = RightFlankingSequence, JunctionSequence),
  transmute(flanking_sequences, ID, Direction = "Reverse",  LeftFlankingSequence, RightFlankingSequence, LeftQuerySequence = RightFlankingSequenceReverseComplement, RightQuerySequence = LeftFlankingSequenceReverseComplement, JunctionSequence = JunctionSequenceReverseComplement)
) %>%
  arrange(ID, Direction) %>%
  mutate(QueryID = str_c(ID, "_", Direction)) %>%
  mutate(QuerySequence = str_c(LeftQuerySequence, ".*", RightQuerySequence))


# read fastq file in chunks

fastq <- gzfile(fastq_file, open = "r")

record_count <- 0
message(Sys.time(), "  ", record_count)

append <- FALSE

repeat {

  fastq_lines <- readLines(fastq, 4000)

  if (length(fastq_lines) == 0) break

  indexes <- seq(1, length(fastq_lines), 4)

  sequences <- tibble(SequenceID = fastq_lines[indexes], Sequence = fastq_lines[indexes + 1]) %>%
    mutate(ID = row_number())

  # extract ID from ID lines
  if (!all(str_detect(sequences$SequenceID, "^@"))) stop("Unexpected format in FASTQ file: ID lines should begin with @")
  sequences <- sequences %>%
    mutate(SequenceID = str_remove(SequenceID, "^@")) %>%
    mutate(SequenceID = str_remove(SequenceID, "[ #].*"))

  # trim UMI tag
  sequences <- sequences %>%
    mutate(UMI = str_sub(Sequence, 1, umi_length)) %>%
    mutate(Sequence = str_sub(Sequence, umi_length + 1))

  # find matches with up to twice the maximum edit distance specified as
  # searching for both left and right flanking sequence in single query
  matches <- query_sequences %>%
    select(query_id = QueryID, query_sequence = QuerySequence) %>%
    pmap_dfr(find_matches, ids = sequences$ID, sequences = sequences$Sequence, max_distance = 2 * max_distance)

  # join matches to left and right flanking sequences and compute distance to
  # the matching sequence for both
  matches <- matches %>%
    select(ID = id, QueryID = query_id, MatchingPosition = position, MatchingSequence = match) %>%
    left_join(sequences, by = "ID") %>%
    select(SequenceID, UMI, Sequence, QueryID, MatchingPosition, MatchingSequence) %>%
    left_join(select(query_sequences, SVID = ID, Direction, QueryID, LeftFlankingSequence, RightFlankingSequence, LeftQuerySequence, RightQuerySequence, JunctionSequence), by = "QueryID") %>%
    rowwise() %>%
    mutate(LeftDistance = as.vector(adist(LeftQuerySequence, MatchingSequence, partial = TRUE))) %>%
    mutate(RightDistance = as.vector(adist(RightQuerySequence, MatchingSequence, partial = TRUE))) %>%
    ungroup()

  # filter distances for each of the left and right flanking sequences that are
  # below the maximum distance specified
  matches <- matches %>%
    filter(LeftDistance <= max_distance, RightDistance <= max_distance)

  # assess match across junction sequence if available
  if (nrow(matches) == 0) {
    matches <- matches %>%
      mutate(JunctionPosition = NA, JunctionMatchingSequence = NA, JunctionDistance = NA, JunctionOverlap = NA)
  } else {
    matches <- matches %>%
      rowwise() %>%
      mutate(matches = aregexec(Sequence, JunctionSequence, max.distance = 0.25, fixed = TRUE)) %>%
      mutate(JunctionPosition = unlist(matches)) %>%
      mutate(JunctionPosition = na_if(JunctionPosition, -1)) %>%
      mutate(JunctionMatchingSequence = ifelse(is.na(JunctionPosition), NA, unlist(regmatches(JunctionSequence, matches)))) %>%
      mutate(JunctionDistance = ifelse(is.na(JunctionPosition), NA, as.vector(adist(Sequence, JunctionMatchingSequence, partial = TRUE)))) %>%
      select(-matches) %>%
      ungroup()

    matches <- matches %>%
      mutate(JunctionMidPoint = floor(nchar(JunctionSequence) / 2)) %>%
      mutate(JunctionOverlap = pmin(JunctionMidPoint - JunctionPosition + 1, JunctionPosition + nchar(JunctionMatchingSequence) - 1 - JunctionMidPoint))
  }

  matches <- matches %>%
    arrange(SequenceID, SVID, Direction) %>%
    mutate(ID = id) %>%
    select(ID, SequenceID, UMI, Sequence, SVID, LeftFlankingSequence, RightFlankingSequence, Direction, MatchingPosition, MatchingSequence, LeftDistance, RightDistance, JunctionPosition, JunctionMatchingSequence, JunctionDistance, JunctionOverlap)

  write_tsv(matches, output_file, append = append)

  record_count <- record_count + nrow(sequences)
  message(Sys.time(), "  ", record_count)

  append <- TRUE
}

close(fastq)
