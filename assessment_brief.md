# Assessment Brief: QHIA Practical Take-Home
## Part 2: Take-Home Practical Assessment

**Due Date:** 9th June 2026  
**Format:** Group work (pairs)  
**Assessment:** Same as Part 1

---

## Overview

Building on the in-class exercises, this assessment requires you to apply lessons from Lectures 1-6 to complete a quantitative health impact assessment (QHIA) in R. You will estimate health impacts of changes in physical activity (PA) for Bogotá's population.

---

## Task Summary

### 1. Derive Baseline and Scenario Exposures
- Use provided travel survey data to calculate baseline PA (MET-hours/week)
- Define scenario: shift from car to active transport (walking/cycling)
- Calculate exposure change by age group and sex

### 2. Calculate PIFs
- Use provided relative risks for health outcomes:
  - Physical activity: IHD, stroke, COPD, diabetes
  - Include citation/source for RRs
- Calculate age-specific PIFs for each outcome

### 3. Estimate Disease Burden Changes
- Apply PIFs to baseline burden from GBD data
- Calculate cases prevented, deaths prevented
- Use MSLT to compute life years gained
- Run separately for subgroups (age groups, sex)

### 4. Uncertainty Analysis
- Conduct sensitivity analysis on:
  - Relative risks (lower/upper CI)
  - Exposure levels
- Present results as ranges

---

## Data Files Provided

| File | Description |
|------|-------------|
| `travel_survey.csv` | Individual trip data with mode, duration, METs |
| `population_bogota.csv` | Population by age/sex |
| `mslt_colombia.csv` | GBD disease rates (incidence, prevalence, mortality) |
| `gbd_disease_burden.csv` | Baseline YLLs/YLDs by disease |

*Contact instructor if files are missing*

---

## Key Functions Available

```r
# From class exercises:
source("functions.R")

# Key functions:
RunLifeTable()       # Build cohort life table
RunDisease()         # Run disease model
calculate_pif()      # Calculate PIF from RR and exposure
```

---

## Deliverables

### 1. R Code (.txt)
- Commented code explaining each step
- Include all calculations from data loading to results

### 2. Report (max 1500 words)

| Section | Content | Approx. Words |
|---------|---------|---------------|
| Introduction | Public health problem, aims | 250 |
| Methods | Data, RR, PIF calculation, analysis | 400 |
| Results | 4 tables/figures max | 500 |
| Discussion | Interpretation, policy, limitations | 350 |

**Total: 1500 words max** (excluding references and code)

---

## Marking Criteria

| Criterion | Weight | Description |
|-----------|--------|-------------|
| Methods | 30% | Correct use of PIF formula, appropriate RR selection |
| Results | 30% | Accurate calculations, clear tables/figures |
| Interpretation | 20% | Meaningful discussion of findings |
| Code Quality | 10% | Well-commented, reproducible |
| Uncertainty | 10% | Appropriate sensitivity analysis |

---

## Questions?

Contact: [Instructor email]
Office hours: [Date/Time]
