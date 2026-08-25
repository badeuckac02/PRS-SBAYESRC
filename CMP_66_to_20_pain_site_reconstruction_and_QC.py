import pandas as pd
import numpy as np

# ============================================================
# 66 ORIGINAL BODY POINTS -> 20 ANATOMICAL REGIONS
# ============================================================

region_mapping = {
    "head":      ["MHEADRF", "MHEADLF", "MHEADLB", "MHEADRB"],
    "neck":      ["MNECKRF", "MNECKLF", "MNECKLB", "MNECKRB"],
    "shoulder":  ["MSHOURF", "MSHOULF", "MSHOULB", "MSHOURB"],
    "upperarm":  ["MUARMRF", "MUARMLF", "MUARMLB", "MUARMRB"],
    "elbow":     ["MELBOWRF", "MELBOWLF", "MELBOWLB", "MELBOWRB"],
    "lowerarm":  ["MLARMRF", "MLARMLF", "MLARMLB", "MLARMRB"],
    "wrist":     ["MWRISTRF", "MWRISTLF", "MWRISTLB", "MWRISTRB"],
    "hand":      ["MHANDRF", "MHANDLF", "MHANDLB", "MHANDRB"],
    "chest":     ["MCHESTR", "MCHESTL"],
    "abdomen":   ["MABDOMR", "MABDOML."],
    "upperback": ["MUBACKL", "MUBACKR"],
    "lowerback": ["MLBACKL", "MLBACKR"],
    "groin":     ["MGROINR", "MGROINL"],
    "buttock":   ["MBUML", "MBUMR"],
    "hip":       ["MHIPR", "MHIPL"],
    "upperleg":  ["MUPLEGRF", "MUPLEGLF", "MUPLEGLB", "MUPLEGRB"],
    "knee":      ["MKNEERF", "MKNEELF", "MKNEELB", "MKNEERB"],
    "lowerleg":  ["MLLEGRF", "MLLEGLF", "MLLEGLB", "MLLEGRB"],
    "ankle":     ["MANKLERF", "MANKLELF", "MANKLELB", "MANKLERB"],
    "foot":      ["MFOOTRF", "MFOOTLF", "MFOOTLB", "MFOOTRB"]
}

raw66 = [
    variable
    for variables in region_mapping.values()
    for variable in variables
]

assert len(raw66) == 66

# ============================================================
# ORIGINAL QUESTIONNAIRE RESPONSE CODING
#
# Like    = pain present
# Neutral = pain absent
# Dislike = pain absent
# Missing = missing
# ============================================================

def pain_binary(value):

    if pd.isna(value):
        return np.nan

    value = str(value).strip()

    if value == "Like":
        return 1.0

    if value in ["Neutral", "Dislike"]:
        return 0.0

    return np.nan


for variable in raw66:
    df["BIN_" + variable] = df[variable].map(pain_binary)


# ============================================================
# COUNT POSITIVE ORIGINAL BODY POINTS
#
# This count is used for CMP multisite eligibility.
# ============================================================

binary66 = ["BIN_" + variable for variable in raw66]

df["raw66_positive_n"] = (
    df[binary66]
    .eq(1)
    .sum(axis=1)
)

df["raw66_missing_n"] = (
    df[binary66]
    .isna()
    .sum(axis=1)
)


# ============================================================
# CMP MULTISITE ELIGIBILITY
#
# Require >=2 positive ORIGINAL body points.
# ============================================================

df["CMP_multisite_eligible"] = (
    df["raw66_positive_n"] >= 2
)


# ============================================================
# COLLAPSE 66 BODY POINTS -> 20 ANATOMICAL REGIONS
#
# OR rule:
# region = 1 if ANY constituent point is "Like"
# region = 0 if no constituent point is "Like"
# region = NA only when ALL constituent responses are missing
# ============================================================

for region, variables in region_mapping.items():

    binary_variables = [
        "BIN_" + variable
        for variable in variables
    ]

    any_positive = (
        df[binary_variables]
        .eq(1)
        .any(axis=1)
    )

    all_missing = (
        df[binary_variables]
        .isna()
        .all(axis=1)
    )

    df[region] = np.where(
        all_missing,
        np.nan,
        any_positive.astype(int)
    )


# ============================================================
# NUMBER OF POSITIVE COLLAPSED ANATOMICAL REGIONS
# ============================================================

regions = list(region_mapping.keys())

df["pain_site_n"] = (
    df[regions]
    .sum(axis=1, min_count=1)
)
