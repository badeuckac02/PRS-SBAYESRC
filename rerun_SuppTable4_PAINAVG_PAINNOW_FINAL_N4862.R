cd /working/lab_miguelr/badeU/04_PRS_AGDS_Painsite

Rscript - <<'RS'

library(MASS)

infile <- "08_LCA_FINAL_N4862/FINAL_LCA10_CASES_PLUS_CONTROLS_N11557_WITH_STANDARDISED_PGS.tsv"
outdir <- "08_LCA_FINAL_N4862"

dat <- read.delim(
  infile,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Final CMP cases only
cases <- dat[dat$analysis_group == "LCA10 case", ]

cat("FINAL CASE N =", nrow(cases), "\n")
print(table(cases$SEX))

# Expected:
# Male   = 1151
# Female = 3711

run_polr <- function(d, outcome, interaction = FALSE) {

  pcs <- paste0("PC", 1:10)

  if (interaction) {

    f <- as.formula(
      paste0(
        "ordered(", outcome, ") ~ CMP_PGS_Z * factor(SEX) + ",
        "AGE + MDD + ",
        paste(pcs, collapse = " + ")
      )
    )

  } else {

    f <- as.formula(
      paste0(
        "ordered(", outcome, ") ~ CMP_PGS_Z + AGE + MDD + ",
        paste(pcs, collapse = " + ")
      )
    )
  }

  MASS::polr(
    f,
    data = d,
    Hess = TRUE,
    method = "logistic"
  )
}

extract_term <- function(model, term) {

  x <- coef(summary(model))

  beta <- x[term, "Value"]
  se   <- x[term, "Std. Error"]

  z <- beta / se
  p <- 2 * pnorm(abs(z), lower.tail = FALSE)

  data.frame(
    beta = beta,
    SE = se,
    OR = exp(beta),
    lower95 = exp(beta - 1.96 * se),
    upper95 = exp(beta + 1.96 * se),
    P = p
  )
}

make_S4 <- function(outcome) {

  female <- cases[cases$SEX == 2, ]
  male   <- cases[cases$SEX == 1, ]

  # Sex-stratified models
  fit_f <- run_polr(female, outcome)
  fit_m <- run_polr(male, outcome)

  # Combined interaction model
  fit_int <- run_polr(cases, outcome, interaction = TRUE)

  pgs_f <- extract_term(fit_f, "CMP_PGS_Z")
  age_f <- extract_term(fit_f, "AGE")

  pgs_m <- extract_term(fit_m, "CMP_PGS_Z")
  age_m <- extract_term(fit_m, "AGE")

  int_term <- grep(
    "CMP_PGS_Z.*factor\\(SEX\\)|factor\\(SEX\\).*CMP_PGS_Z",
    names(coef(fit_int)),
    value = TRUE
  )

  int_res <- extract_term(fit_int, int_term)

  vars <- c(
    outcome, "CMP_PGS_Z", "AGE", "MDD",
    paste0("PC", 1:10)
  )

  n_f <- sum(complete.cases(female[, vars]))
  n_m <- sum(complete.cases(male[, vars]))

  tab <- data.frame(
    Sex = c("Female", "Male"),
    N = c(n_f, n_m),

    PGS_OR = c(pgs_f$OR, pgs_m$OR),
    PGS_CI_L = c(pgs_f$lower95, pgs_m$lower95),
    PGS_CI_U = c(pgs_f$upper95, pgs_m$upper95),
    PGS_P = c(pgs_f$P, pgs_m$P),

    Age_OR = c(age_f$OR, age_m$OR),
    Age_CI_L = c(age_f$lower95, age_m$lower95),
    Age_CI_U = c(age_f$upper95, age_m$upper95),
    Age_P = c(age_f$P, age_m$P),

    PGS_sex_interaction_OR = c(int_res$OR, NA),
    PGS_sex_interaction_CI_L = c(int_res$lower95, NA),
    PGS_sex_interaction_CI_U = c(int_res$upper95, NA),
    PGS_sex_interaction_P = c(int_res$P, NA)
  )

  return(list(
    table = tab,
    female_model = fit_f,
    male_model = fit_m,
    interaction_model = fit_int
  ))
}

S4A <- make_S4("PAINAVG")
S4B <- make_S4("PAINNOW")

write.table(
  S4A$table,
  file.path(
    outdir,
    "Supplementary_Table_4A_AVERAGE_PAIN_FINAL_N4862.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  S4B$table,
  file.path(
    outdir,
    "Supplementary_Table_4B_CURRENT_PAIN_FINAL_N4862.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Complete model output for Miguel and Santiago
sink(
  file.path(
    outdir,
    "Supplementary_Table_4_COMPLETE_MODEL_OUTPUT_FINAL_N4862.txt"
  )
)

cat("============================================\n")
cat("SUPPLEMENTARY TABLE 4A — PAINAVG\n")
cat("============================================\n\n")

cat("Female model\n")
print(summary(S4A$female_model))

cat("\nMale model\n")
print(summary(S4A$male_model))

cat("\nCombined PGS x SEX interaction model\n")
print(summary(S4A$interaction_model))

cat("\n\n============================================\n")
cat("SUPPLEMENTARY TABLE 4B — PAINNOW\n")
cat("============================================\n\n")

cat("Female model\n")
print(summary(S4B$female_model))

cat("\nMale model\n")
print(summary(S4B$male_model))

cat("\nCombined PGS x SEX interaction model\n")
print(summary(S4B$interaction_model))

sink()

cat("\nS4A\n")
print(S4A$table, digits = 10)

cat("\nS4B\n")
print(S4B$table, digits = 10)

cat("\nPAINAVG raw available responses:\n")
print(table(cases$SEX[!is.na(cases$PAINAVG)]))

cat("\nPAINNOW raw available responses:\n")
print(table(cases$SEX[!is.na(cases$PAINNOW)]))

RS
