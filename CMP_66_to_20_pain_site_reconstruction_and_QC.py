cd /working/lab_miguelr/badeU/04_PRS_AGDS_Painsite

mkdir -p 12_FINAL_LCA_N6216_CASES/{00_INPUT,01_SCRIPTS,02_REGION_COLLAPSE,03_MODEL_FITS,04_MODEL_SELECTION,05_FINAL_SOLUTION,06_FIGURES,07_TABLES,08_LOGS}

cat > 12_FINAL_LCA_N6216_CASES/01_SCRIPTS/01_COLLAPSE_66_TO_20_REGIONS.R <<'RS'

library(readr)
library(dplyr)

# ============================================================
# PATHS
# ============================================================

project <- "/working/lab_miguelr/badeU/04_PRS_AGDS_Painsite"

master_file <- file.path(
  project,
  "11_FINAL_RECONSTRUCTED_MASTER_N11466",
  "FINAL_ANALYSIS_READY_AGDS_MASTER_N10804.tsv"
)

outdir <- file.path(
  project,
  "12_FINAL_LCA_N6216_CASES"
)

# ============================================================
# 66 ORIGINAL BODY POINTS -> 20 ANATOMICAL REGIONS
# ============================================================

region_mapping <- list(

  head = c(
    "MHEADRF", "MHEADLF", "MHEADLB", "MHEADRB"
  ),

  neck = c(
    "MNECKRF", "MNECKLF", "MNECKLB", "MNECKRB"
  ),

  shoulder = c(
    "MSHOURF", "MSHOULF", "MSHOULB", "MSHOURB"
  ),

  upperarm = c(
    "MUARMRF", "MUARMLF", "MUARMLB", "MUARMRB"
  ),

  elbow = c(
    "MELBOWRF", "MELBOWLF", "MELBOWLB", "MELBOWRB"
  ),

  lowerarm = c(
    "MLARMRF", "MLARMLF", "MLARMLB", "MLARMRB"
  ),

  wrist = c(
    "MWRISTRF", "MWRISTLF", "MWRISTLB", "MWRISTRB"
  ),

  hand = c(
    "MHANDRF", "MHANDLF", "MHANDLB", "MHANDRB"
  ),

  chest = c(
    "MCHESTR", "MCHESTL"
  ),

  abdomen = c(
    "MABDOMR", "MABDOML."
  ),

  upperback = c(
    "MUBACKL", "MUBACKR"
  ),

  lowerback = c(
    "MLBACKL", "MLBACKR"
  ),

  groin = c(
    "MGROINR", "MGROINL"
  ),

  buttock = c(
    "MBUML", "MBUMR"
  ),

  hip = c(
    "MHIPR", "MHIPL"
  ),

  upperleg = c(
    "MUPLEGRF", "MUPLEGLF", "MUPLEGLB", "MUPLEGRB"
  ),

  knee = c(
    "MKNEERF", "MKNEELF", "MKNEELB", "MKNEERB"
  ),

  lowerleg = c(
    "MLLEGRF", "MLLEGLF", "MLLEGLB", "MLLEGRB"
  ),

  ankle = c(
    "MANKLERF", "MANKLELF", "MANKLELB", "MANKLERB"
  ),

  foot = c(
    "MFOOTRF", "MFOOTLF", "MFOOTLB", "MFOOTRB"
  )
)

region_names <- names(region_mapping)
body_vars <- unique(unlist(region_mapping))

# ============================================================
# LOAD MASTER
# ============================================================

dat <- read_tsv(
  master_file,
  show_col_types = FALSE
)

cat("\n============================================================\n")
cat("MASTER DATASET\n")
cat("============================================================\n")

cat("N:", nrow(dat), "\n")
cat("Cases:", sum(dat$analysis_status == "Case"), "\n")
cat("Controls:", sum(dat$analysis_status == "Control"), "\n")

stopifnot(nrow(dat) == 10804)
stopifnot(sum(dat$analysis_status == "Case") == 6216)
stopifnot(sum(dat$analysis_status == "Control") == 4588)

# ============================================================
# VERIFY ALL 66 VARIABLES EXIST
# ============================================================

missing_vars <- setdiff(body_vars, names(dat))

if (length(missing_vars) > 0) {

  cat("\nMISSING BODY VARIABLES:\n")
  print(missing_vars)

  stop("Not all 66 body-point variables are present.")
}

cat("\nNumber of original body-point variables:", length(body_vars), "\n")

stopifnot(length(body_vars) == 66)

# ============================================================
# CASES ONLY
# ============================================================

cases <- dat %>%
  filter(analysis_status == "Case")

stopifnot(nrow(cases) == 6216)

# ============================================================
# VERIFY RAW BODY-MAP VALUES
#
# Expected:
# Like
# Neutral
# Dislike
# blank / NA
# ============================================================

all_values <- sort(
  unique(
    unlist(
      cases[body_vars],
      use.names = FALSE
    )
  )
)

cat("\n============================================================\n")
cat("RAW BODY-MAP VALUES IN N=6,216 CASES\n")
cat("============================================================\n")

print(all_values)

# ============================================================
# RECOUNT ORIGINAL POSITIVE BODY POINTS
#
# ONLY "Like" = positive
# ============================================================

like_matrix <- sapply(
  cases[body_vars],
  function(x) {
    !is.na(x) & trimws(as.character(x)) == "Like"
  }
)

cases$n_like_recalculated <- rowSums(like_matrix)

cat("\n============================================================\n")
cat("ORIGINAL 66-POINT POSITIVE COUNT\n")
cat("============================================================\n")

print(table(cases$n_like_recalculated))

cat(
  "\nCases with <2 Like points:",
  sum(cases$n_like_recalculated < 2),
  "\n"
)

cat(
  "Minimum Like points:",
  min(cases$n_like_recalculated),
  "\n"
)

cat(
  "Maximum Like points:",
  max(cases$n_like_recalculated),
  "\n"
)

# Current case definition requires >=2 original Like points.
stopifnot(all(cases$n_like_recalculated >= 2))

# ============================================================
# COLLAPSE TO 20 BINARY REGIONS
#
# Region = 1 if ANY constituent body point == "Like"
# Region = 0 otherwise.
#
# Neutral and Dislike are NOT positive.
# ============================================================

for (region in region_names) {

  vars <- region_mapping[[region]]

  region_like <- sapply(
    cases[vars],
    function(x) {
      !is.na(x) & trimws(as.character(x)) == "Like"
    }
  )

  if (is.vector(region_like)) {
    region_like <- matrix(region_like, ncol = 1)
  }

  cases[[region]] <- as.integer(
    rowSums(region_like) > 0
  )
}

# ============================================================
# VERIFY 20 REGION VARIABLES
# ============================================================

cat("\n============================================================\n")
cat("20-REGION QC\n")
cat("============================================================\n")

for (v in region_names) {

  cat("\n", v, "\n", sep = "")
  print(table(cases[[v]], useNA = "ifany"))

  stopifnot(
    all(cases[[v]] %in% c(0L, 1L))
  )

  stopifnot(
    sum(is.na(cases[[v]])) == 0
  )
}

# ============================================================
# NUMBER OF POSITIVE COLLAPSED REGIONS
# ============================================================

cases$n_regions_positive <- rowSums(
  cases[region_names]
)

cat("\n============================================================\n")
cat("NUMBER OF POSITIVE 20-REGION INDICATORS\n")
cat("============================================================\n")

print(
  table(cases$n_regions_positive)
)

cat(
  "\nMinimum positive regions:",
  min(cases$n_regions_positive),
  "\n"
)

cat(
  "Maximum positive regions:",
  max(cases$n_regions_positive),
  "\n"
)

# IMPORTANT:
# Do NOT exclude participants with one collapsed region.
#
# A participant can have >=2 original Like body points
# within one anatomical region.
# ============================================================

cat(
  "\nParticipants with exactly 1 collapsed region:",
  sum(cases$n_regions_positive == 1),
  "\n"
)

# ============================================================
# REGION PREVALENCE
# ============================================================

region_prevalence <- tibble(
  region = region_names,
  N_positive = sapply(
    cases[region_names],
    sum
  ),
  N_total = nrow(cases)
) %>%
  mutate(
    prevalence_percent =
      100 * N_positive / N_total
  )

cat("\n============================================================\n")
cat("REGION PREVALENCE\n")
cat("============================================================\n")

print(region_prevalence, n = 20)

# ============================================================
# LCA INPUT
# ============================================================

lca_input <- cases %>%
  select(
    IID,
    STUDYID,
    all_of(region_names)
  )

stopifnot(nrow(lca_input) == 6216)
stopifnot(n_distinct(lca_input$IID) == 6216)
stopifnot(n_distinct(lca_input$STUDYID) == 6216)

# ============================================================
# SAVE
# ============================================================

write_tsv(
  cases,
  file.path(
    outdir,
    "02_REGION_COLLAPSE",
    "FINAL_N6216_CASES_WITH_20_PAIN_REGIONS.tsv"
  )
)

write_tsv(
  lca_input,
  file.path(
    outdir,
    "00_INPUT",
    "LCA_INPUT_N6216_20_BINARY_REGIONS.tsv"
  )
)

write_tsv(
  region_prevalence,
  file.path(
    outdir,
    "02_REGION_COLLAPSE",
    "N6216_20_REGION_PREVALENCE.tsv"
  )
)

mapping_table <- bind_rows(
  lapply(
    names(region_mapping),
    function(region) {
      tibble(
        region = region,
        original_body_point =
          region_mapping[[region]]
      )
    }
  )
)

write_tsv(
  mapping_table,
  file.path(
    outdir,
    "02_REGION_COLLAPSE",
    "BODYPOINT_66_TO_REGION_20_MAPPING.tsv"
  )
)

cat("\n============================================================\n")
cat("REGION COLLAPSE COMPLETE\n")
cat("============================================================\n")

cat("Cases:", nrow(cases), "\n")
cat("Original body points:", length(body_vars), "\n")
cat("Collapsed regions:", length(region_names), "\n")

cat("\nLCA input:\n")
cat(
  file.path(
    outdir,
    "00_INPUT",
    "LCA_INPUT_N6216_20_BINARY_REGIONS.tsv"
  ),
  "\n"
)

RS

Rscript 12_FINAL_LCA_N6216_CASES/01_SCRIPTS/01_COLLAPSE_66_TO_20_REGIONS.R \
  2>&1 | tee 12_FINAL_LCA_N6216_CASES/08_LOGS/01_COLLAPSE_66_TO_20_REGIONS.log
