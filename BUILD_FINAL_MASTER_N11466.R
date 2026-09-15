
# ======================================================================
# BUILD ENRICHED RECONSTRUCTED AGDS MASTER COHORT
# ======================================================================
#
# Starting reconstructed cohort:
#   N = 11,466
#   Cases = 6,570
#   Controls = 4,896
#
# Add from original AGDS .sav:
#   TYPE1A = osteoarthritis
#   20 psychiatric diagnosis variables
#
# Add from PRS profile:
#   PRS_total
#
# IMPORTANT:
#   - Does NOT alter the original reconstructed dataset.
#   - Does NOT exclude participants because psychiatric/OA variables
#     are missing.
#   - Does NOT exclude participants because PRS_total is missing.
#   - Hard QC checks prevent accidental row duplication.
#
# ======================================================================

library(haven)
library(readr)
library(dplyr)

# ======================================================================
# PATHS
# ======================================================================

project <- "/working/lab_miguelr/badeU/04_PRS_AGDS_Painsite"

cohort_file <- file.path(
  project,
  "10_AUDIT/PARTICIPANT_FLOW_RECONSTRUCTION/FINAL_N11466",
  "FINAL_AGDS_CASE_CONTROL_N11466_FULL_PHENOTYPE.tsv"
)

agds_file <- paste0(
  "/working/lab_miguelr/badeU/Reference/AGDS/",
  "Combined_AGDS_2020_Freeze_20251205.sav"
)

prs_file <- file.path(
  project,
  "04_prs_models/profiles/PRS_profile",
  "AGDS_only_PRS_total.txt"
)

outdir <- file.path(
  project,
  "11_FINAL_RECONSTRUCTED_MASTER_N11466"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)


# ======================================================================
# VARIABLES TO ADD FROM ORIGINAL AGDS
# ======================================================================

psych_vars <- c(
  "DXMDD",
  "DXBPD",
  "DXPDMD",
  "DXSCZ",
  "DXANOR",
  "DXBUL",
  "DXADHD",
  "DXASD",
  "DXTOUR",
  "DXANX",
  "DXPANIC",
  "DXOCD",
  "DXHOARD",
  "DXPTSD",
  "DXPHOB",
  "DXSAD",
  "DXSANX",
  "DXAGORA",
  "DXPERSD",
  "DXSUD"
)

clinical_vars <- c(
  "TYPE1A",
  psych_vars
)


# ======================================================================
# STEP 1 — LOAD RECONSTRUCTED N=11,466 COHORT
# ======================================================================

cat("\n============================================================\n")
cat("STEP 1 — LOAD RECONSTRUCTED COHORT\n")
cat("============================================================\n")

cohort <- read_tsv(
  cohort_file,
  show_col_types = FALSE
)

cat("Rows:", nrow(cohort), "\n")
cat("Columns:", ncol(cohort), "\n")
cat("Cases:", sum(cohort$analysis_status == "Case"), "\n")
cat("Controls:", sum(cohort$analysis_status == "Control"), "\n")

stopifnot(nrow(cohort) == 11466)
stopifnot(sum(cohort$analysis_status == "Case") == 6570)
stopifnot(sum(cohort$analysis_status == "Control") == 4896)

cohort <- cohort %>%
  mutate(
    STUDYID_MATCH = trimws(as.character(STUDYID)),
    IID_MATCH = trimws(as.character(IID))
  )

stopifnot(n_distinct(cohort$STUDYID_MATCH) == 11466)
stopifnot(n_distinct(cohort$IID_MATCH) == 11466)


# ======================================================================
# STEP 2 — LOAD REQUIRED VARIABLES FROM ORIGINAL AGDS SAV
# ======================================================================

cat("\n============================================================\n")
cat("STEP 2 — LOAD CLINICAL VARIABLES FROM ORIGINAL AGDS\n")
cat("============================================================\n")

# Read source
agds <- read_sav(agds_file)

cat("Original AGDS rows:", nrow(agds), "\n")
cat("Original AGDS columns:", ncol(agds), "\n")

stopifnot(nrow(agds) == 22374)


# ======================================================================
# VERIFY ALL REQUESTED VARIABLES EXIST
# ======================================================================

required_agds <- c(
  "STUDYID",
  "TYPE1A",
  psych_vars
)

missing_agds_vars <- setdiff(
  required_agds,
  names(agds)
)

if (length(missing_agds_vars) > 0) {

  cat("\nERROR — VARIABLES MISSING FROM AGDS SOURCE:\n")
  print(missing_agds_vars)

  stop(
    "One or more requested variables are absent from the original AGDS SAV."
  )
}

cat("\nAll requested clinical variables found.\n")


# ======================================================================
# STEP 3 — EXTRACT CLINICAL VARIABLES
# ======================================================================

clinical <- agds %>%
  select(
    STUDYID,
    all_of(clinical_vars)
  ) %>%
  mutate(
    STUDYID_MATCH =
      trimws(as.character(STUDYID))
  )

cat(
  "Unique AGDS STUDYIDs:",
  n_distinct(clinical$STUDYID_MATCH),
  "\n"
)


# ======================================================================
# CHECK DUPLICATE SOURCE STUDYIDS
# ======================================================================

clinical_duplicates <- clinical %>%
  count(STUDYID_MATCH) %>%
  filter(n > 1)

cat(
  "Duplicated STUDYIDs in original source:",
  nrow(clinical_duplicates),
  "\n"
)

if (nrow(clinical_duplicates) > 0) {

  write_tsv(
    clinical_duplicates,
    file.path(
      outdir,
      "QC_DUPLICATED_STUDYIDS_IN_ORIGINAL_AGDS.tsv"
    )
  )

  # Our reconstructed cohort already contains unique participant IDs.
  # Restrict source to IDs required by final cohort, then retain one
  # source record per STUDYID only where required for the join.

  clinical <- clinical %>%
    filter(
      STUDYID_MATCH %in% cohort$STUDYID_MATCH
    ) %>%
    distinct(
      STUDYID_MATCH,
      .keep_all = TRUE
    )
}


# ======================================================================
# STEP 4 — JOIN OA + PSYCHIATRIC VARIABLES
# ======================================================================

cat("\n============================================================\n")
cat("STEP 4 — ADD OA + PSYCHIATRIC DIAGNOSES\n")
cat("============================================================\n")

n_before <- nrow(cohort)

master <- cohort %>%
  left_join(
    clinical %>%
      select(
        -STUDYID
      ),
    by = "STUDYID_MATCH"
  )

stopifnot(nrow(master) == n_before)

cat("Rows after clinical join:", nrow(master), "\n")


# ======================================================================
# REPORT COMPLETENESS
# ======================================================================

clinical_qc <- tibble(
  variable = clinical_vars,
  N_nonmissing = sapply(
    master[clinical_vars],
    function(x) sum(!is.na(x))
  ),
  N_missing = sapply(
    master[clinical_vars],
    function(x) sum(is.na(x))
  )
)

clinical_qc <- clinical_qc %>%
  mutate(
    percent_nonmissing =
      100 * N_nonmissing / nrow(master)
  )

cat("\nClinical variable completeness:\n")

print(
  clinical_qc,
  n = Inf
)

write_tsv(
  clinical_qc,
  file.path(
    outdir,
    "QC_CLINICAL_VARIABLE_COMPLETENESS.tsv"
  )
)


# ======================================================================
# STEP 5 — LOAD PRS_total FILE
# ======================================================================

cat("\n============================================================\n")
cat("STEP 5 — LOAD PRS_total\n")
cat("============================================================\n")

# First inspect delimiter/header automatically
prs <- read.table(
  prs_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("PRS rows:", nrow(prs), "\n")
cat("PRS columns:", ncol(prs), "\n")

cat("\nPRS column names:\n")
print(names(prs))


# ======================================================================
# IDENTIFY PRS VALUE COLUMN
# ======================================================================

if (!"PRS_total" %in% names(prs)) {

  cat("\nERROR: PRS_total not found.\n")
  cat("Available columns:\n")

  print(names(prs))

  stop(
    "PRS_total column is not present in PRS profile file."
  )
}


# ======================================================================
# STEP 6 — DETERMINE SAFE PARTICIPANT IDENTIFIER FOR PRS
#
# Preferred:
#   IID
#
# Fallback:
#   STUDYID
#
# We DO NOT guess another identifier.
# ======================================================================

if ("IID" %in% names(prs)) {

  prs_id <- "IID"

  prs <- prs %>%
    mutate(
      PRS_MATCH_ID =
        trimws(as.character(IID))
    )

  master <- master %>%
    mutate(
      PRS_MATCH_ID = IID_MATCH
    )

} else if ("STUDYID" %in% names(prs)) {

  prs_id <- "STUDYID"

  prs <- prs %>%
    mutate(
      PRS_MATCH_ID =
        trimws(as.character(STUDYID))
    )

  master <- master %>%
    mutate(
      PRS_MATCH_ID = STUDYID_MATCH
    )

} else {

  cat("\nERROR: Cannot safely identify participant ID in PRS file.\n")
  cat("PRS columns are:\n")

  print(names(prs))

  stop(
    paste0(
      "PRS file contains neither IID nor STUDYID. ",
      "No PRS values were joined."
    )
  )
}

cat("\nPRS matched using:", prs_id, "\n")


# ======================================================================
# STEP 7 — PRS DUPLICATE QC
# ======================================================================

prs_duplicates <- prs %>%
  filter(
    !is.na(PRS_MATCH_ID),
    PRS_MATCH_ID != ""
  ) %>%
  count(PRS_MATCH_ID) %>%
  filter(n > 1)

cat(
  "Duplicated participant IDs in PRS file:",
  nrow(prs_duplicates),
  "\n"
)

if (nrow(prs_duplicates) > 0) {

  write_tsv(
    prs_duplicates,
    file.path(
      outdir,
      "QC_DUPLICATED_IDS_IN_PRS_FILE.tsv"
    )
  )

  stop(
    paste0(
      "Duplicate participant IDs detected in PRS file. ",
      "Stopped rather than arbitrarily selecting a PRS value."
    )
  )
}


# ======================================================================
# STEP 8 — ADD PRS_total
# ======================================================================

prs_join <- prs %>%
  select(
    PRS_MATCH_ID,
    PRS_total
  )

n_before_prs <- nrow(master)

master <- master %>%
  left_join(
    prs_join,
    by = "PRS_MATCH_ID"
  )

stopifnot(nrow(master) == n_before_prs)


# ======================================================================
# PRS COMPLETENESS
# ======================================================================

n_prs <- sum(!is.na(master$PRS_total))
n_no_prs <- sum(is.na(master$PRS_total))

cat("\n============================================================\n")
cat("PRS MATCHING RESULTS\n")
cat("============================================================\n")

cat("Final cohort N:", nrow(master), "\n")
cat("PRS_total available:", n_prs, "\n")
cat("PRS_total missing:", n_no_prs, "\n")

cat(
  "PRS coverage:",
  round(100 * n_prs / nrow(master), 2),
  "%\n"
)


# ======================================================================
# PRS COVERAGE BY CASE/CONTROL
# ======================================================================

prs_by_status <- master %>%
  group_by(analysis_status) %>%
  summarise(
    N = n(),
    PRS_available = sum(!is.na(PRS_total)),
    PRS_missing = sum(is.na(PRS_total)),
    PRS_coverage_percent =
      100 * PRS_available / N,
    .groups = "drop"
  )

cat("\nPRS coverage by case/control:\n")

print(prs_by_status)

write_tsv(
  prs_by_status,
  file.path(
    outdir,
    "QC_PRS_COVERAGE_BY_CASE_CONTROL.tsv"
  )
)


# ======================================================================
# STEP 9 — REMOVE TEMPORARY MATCHING VARIABLES
# ======================================================================

master <- master %>%
  select(
    -STUDYID_MATCH,
    -IID_MATCH,
    -PRS_MATCH_ID
  )


# ======================================================================
# STEP 10 — PUT IMPORTANT VARIABLES NEAR FRONT
# ======================================================================

master <- master %>%
  relocate(
    FID,
    IID,
    SUPPLIED_ID,
    STUDYID,
    analysis_status,
    AGE,
    SEX,
    PAINMAIN,
    PAINAVG,
    TYPE1A,
    all_of(psych_vars),
    PRS_total
  )


# ======================================================================
# STEP 11 — FINAL HARD QC
# ======================================================================

cat("\n============================================================\n")
cat("FINAL MASTER DATASET QC\n")
cat("============================================================\n")

cat("Rows:", nrow(master), "\n")
cat("Columns:", ncol(master), "\n")

cat(
  "Unique STUDYID:",
  n_distinct(master$STUDYID),
  "\n"
)

cat(
  "Unique IID:",
  n_distinct(master$IID),
  "\n"
)

cat(
  "Cases:",
  sum(master$analysis_status == "Case"),
  "\n"
)

cat(
  "Controls:",
  sum(master$analysis_status == "Control"),
  "\n"
)

cat(
  "PRS_total nonmissing:",
  sum(!is.na(master$PRS_total)),
  "\n"
)

stopifnot(nrow(master) == 11466)
stopifnot(n_distinct(master$STUDYID) == 11466)
stopifnot(n_distinct(master$IID) == 11466)

stopifnot(
  sum(master$analysis_status == "Case") ==
    6570
)

stopifnot(
  sum(master$analysis_status == "Control") ==
    4896
)


# ======================================================================
# STEP 12 — SAVE MASTER TSV
# ======================================================================

master_tsv <- file.path(
  outdir,
  "FINAL_RECONSTRUCTED_AGDS_MASTER_N11466_WITH_OA_PSYCH_PRS.tsv"
)

write_tsv(
  master,
  master_tsv,
  na = "NA"
)


# ======================================================================
# STEP 13 — SAVE SPSS VERSION
#
# SPSS does not permit MABDOML.
# ======================================================================

master_sav <- master

if ("MABDOML." %in% names(master_sav)) {

  names(master_sav)[
    names(master_sav) == "MABDOML."
  ] <- "MABDOML"

}

sav_file <- file.path(
  outdir,
  "FINAL_RECONSTRUCTED_AGDS_MASTER_N11466_WITH_OA_PSYCH_PRS.sav"
)

write_sav(
  master_sav,
  sav_file
)


# ======================================================================
# STEP 14 — SAVE PARTICIPANTS WITHOUT PRS
# ======================================================================

master %>%
  filter(is.na(PRS_total)) %>%
  select(
    FID,
    IID,
    SUPPLIED_ID,
    STUDYID,
    analysis_status
  ) %>%
  write_tsv(
    file.path(
      outdir,
      "PARTICIPANTS_WITHOUT_PRS_TOTAL.tsv"
    )
  )


# ======================================================================
# STEP 15 — SAVE MASTER SUMMARY
# ======================================================================

summary <- tibble(
  metric = c(
    "Total participants",
    "Cases",
    "Controls",
    "PRS_total available",
    "PRS_total missing",
    "TYPE1A nonmissing",
    "Missing PAINMAIN retained",
    "Missing PAINAVG"
  ),

  N = c(
    nrow(master),
    sum(master$analysis_status == "Case"),
    sum(master$analysis_status == "Control"),
    sum(!is.na(master$PRS_total)),
    sum(is.na(master$PRS_total)),
    sum(!is.na(master$TYPE1A)),
    sum(is.na(master$PAINMAIN)),
    sum(is.na(master$PAINAVG))
  )
)

write_tsv(
  summary,
  file.path(
    outdir,
    "FINAL_MASTER_DATASET_SUMMARY.tsv"
  )
)


# ======================================================================
# STEP 16 — READ BACK TSV AND SAV
# ======================================================================

check_tsv <- read_tsv(
  master_tsv,
  show_col_types = FALSE
)

check_sav <- read_sav(
  sav_file
)

stopifnot(nrow(check_tsv) == 11466)
stopifnot(nrow(check_sav) == 11466)


# ======================================================================
# COMPLETE
# ======================================================================

cat("\n\n")
cat("============================================================\n")
cat("MASTER DATASET CREATED SUCCESSFULLY\n")
cat("============================================================\n")

cat("\nFinal N:", nrow(master), "\n")

cat(
  "Cases:",
  sum(master$analysis_status == "Case"),
  "\n"
)

cat(
  "Controls:",
  sum(master$analysis_status == "Control"),
  "\n"
)

cat(
  "PRS_total available:",
  sum(!is.na(master$PRS_total)),
  "\n"
)

cat(
  "PRS_total missing:",
  sum(is.na(master$PRS_total)),
  "\n"
)

cat("\nMASTER TSV:\n")
cat(master_tsv, "\n")

cat("\nMASTER SAV:\n")
cat(sav_file, "\n")

cat("\n============================================================\n")

