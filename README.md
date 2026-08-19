# Final Silwood biodiversity-uncertainty project

## Project title

**Validation-calibrated propagation of automated species-identification uncertainty into Hill diversity and temporal partitioning**

This directory is the submission-ready computational archive for the MSc project. It is intended to let a supervisor, examiner, or future user trace the study from the original BirdNET export through data cleaning, manual-validation inputs, calibration modelling, full Monte Carlo uncertainty propagation, diversity calculations, robustness analyses, analytical benchmarking, and the exact figures and tables used in the thesis.

The archive contains the data, R code, computational outputs, and run logs needed to understand and reproduce the analysis. Literature PDFs, thesis drafts, supervisor feedback, and other writing materials are not part of the analytical archive.

---

# 1. What this project does

The project addresses a specific problem in automated biodiversity monitoring: a BirdNET confidence score is a classifier output, not a validated probability that a predicted species label is correct. The analysis therefore does not treat BirdNET confidence scores as probabilities and does not modify the mathematical definitions of Hill diversity.

Instead, the workflow:

1. cleans the BirdNET detection export;
2. defines an **event-species candidate** as one predicted-species label associated with one acoustic event;
3. uses a manually reviewed validation dataset to calibrate candidate correctness;
4. compares pooled and species-aware logistic calibration models;
5. assigns a **candidate-level correctness probability** to every cleaned event-species candidate;
6. propagates classification and calibration uncertainty through a full Monte Carlo framework;
7. calculates monthly Hill diversity at \(q=0\), \(q=1\), and \(q=2\), together with temporal alpha diversity, multiplicative beta diversity, and gamma diversity;
8. compares the probabilistic framework with naive retention, hard top-1 assignment, and three global BirdNET confidence-score thresholds;
9. evaluates sensitivity to represented sampling coverage; and
10. compares the full simulation framework with an analytical / plug-in benchmark.

Three interpretation boundaries are central throughout the project:

- **BirdNET confidence score is not a probability of correct identification.**
- **An event-species candidate is not an individual bird.**
- **The framework changes how uncertain automated identifications enter established biodiversity calculations; it does not redefine Hill diversity or multiplicative diversity partitioning.**

---

# 2. Reproducibility status

A complete final run of `run_all.R` was successfully executed on **19 August 2026**.

Reference run:

- start: **2026-08-19 04:53:54**
- completion: **2026-08-19 05:02:16**
- mode: **final**
- Monte Carlo simulations per probability scenario: **1,000**
- fixed-probability Monte Carlo simulations: **1,000**
- event-grouped cross-validation: **5 folds x 5 repeats**
- total recorded pipeline runtime: approximately **502 seconds (8 min 22 sec)** on the reference machine
- all default pipeline steps completed successfully

The most computationally expensive stages in the reference run were:

- full Monte Carlo uncertainty propagation: approximately **284 sec**
- analytical benchmark including fixed-probability Monte Carlo: approximately **92 sec**
- repeated event-grouped cross-validation: approximately **21 sec**

The exact execution record is stored in:

```text
log/run_20260819_045354/
```

Stable copies of the most recent successful run are stored as:

```text
log/latest_pipeline_log.txt
log/latest_pipeline_manifest.csv
log/latest_sessionInfo.txt
```

---

# 3. Quick start

## 3.1 Open the project

Open:

```text
silwood_project.Rproj.Rproj
```

in RStudio.

Alternatively, start R from the `final_project/` root directory.

`run_all.R` checks that it is being executed from the project root.

## 3.2 Install the required packages

Run once:

```r
source("install_packages.R")
```

or from a shell:

```bash
Rscript install_packages.R
```

The installer checks for and, if necessary, installs:

```text
tidyverse
here
lme4
Matrix
MASS
patchwork
lubridate
scales
ragg
digest
```

## 3.3 Rebuild the complete final analysis

From R:

```r
source("run_all.R")
```

or from a shell:

```bash
Rscript run_all.R
```

The default settings reproduce the final thesis analysis:

```text
SILWOOD_MODE=final
1,000 Monte Carlo simulations per probability scenario
1,000 fixed-probability Monte Carlo simulations
5-fold event-grouped cross-validation
5 cross-validation repeats
manual-validation decisions treated as frozen input
LOSO diagnostic disabled
final reporting outputs enabled
```

When the run finishes successfully, the exact manuscript outputs are available in:

```text
results/10_final_outputs/
```

---

# 4. Important reproducibility boundary: manual validation

Manual listening and ecological / biogeographic review are human judgements and cannot be recreated automatically from R code.

The final reviewed validation decisions are therefore treated as a **frozen analytical input**:

```text
data/02_manual_validation/final/
validation_sample_silwood_combined_final.csv
```

The project also preserves the completed Round-1 and Round-2 review files for audit.

The sampling-design scripts can reproduce how records were selected for review, but they cannot recreate the human decisions themselves. For this reason:

- scripts `02`, `10`, and `11` are audit / sampling-design scripts;
- they are not run in the default final pipeline;
- `run_all.R` does not overwrite the frozen final validation decisions.

This is intentional and is necessary to distinguish reproducible computational processing from irreducible human validation.

---

# 5. Directory structure

```text
final_project/
|
|-- README.md
|-- FILE_MANIFEST.csv
|-- silwood_project.Rproj.Rproj
|-- install_packages.R
|-- run_all.R
|
|-- data/
|   |-- 00_raw/
|   |-- 01_cleaned/
|   |-- 02_manual_validation/
|   |   |-- round1_design/
|   |   |-- round1_completed/
|   |   |-- round2_design/
|   |   |-- round2_completed/
|   |   `-- final/
|   |-- 03_analysis_ready/
|   |   |-- validation/
|   |   `-- calibration/
|   `-- 04_static_assets/
|
|-- scripts/
|   |-- 00_setup/
|   |-- 01_data_preparation/
|   |-- 02_manual_validation/
|   |-- 03_calibration/
|   |-- 04_uncertainty_propagation/
|   |-- 05_method_comparison/
|   |-- 06_sampling_coverage/
|   |-- 07_analytical_benchmark/
|   |-- 08_reporting/
|   `-- 99_optional_diagnostics/
|
|-- results/
|   |-- 00_data_preparation/
|   |-- 01_validation/
|   |-- 02_calibration/
|   |-- 03_model_validation/
|   |-- 04_probability_assignment/
|   |-- 05_monte_carlo/
|   |-- 06_method_comparison/
|   |-- 07_sampling_coverage/
|   |-- 08_coverage_robustness/
|   |-- 09_analytical_benchmark/
|   `-- 10_final_outputs/
|
`-- log/
    |-- reference_run/
    |-- run_YYYYMMDD_HHMMSS/
    |-- latest_pipeline_log.txt
    |-- latest_pipeline_manifest.csv
    `-- latest_sessionInfo.txt
```

The numbering of the `results/` folders follows the scientific workflow. The script numbers retain the original analysis-step numbers so that individual scripts can still be traced back to the development history of the project.

---

# 6. Data folders and what they contain

## 6.1 `data/00_raw/` — original detector export

Primary raw input:

```text
detections-export-proj_silwood_park(in).csv
```

This is the original Silwood Park BirdNET detection export supplied for the project.

Verified input size:

- raw export rows: **79,398**

The raw export contains both BirdNET species detections and anomaly-detection rows. The anomaly rows are retained in the raw source file but excluded from species-diversity analysis.

---

## 6.2 `data/01_cleaned/` — cleaned event-species candidate data

Generated by:

```text
scripts/01_data_preparation/01_clean_birdnet_data.R
```

Key files:

```text
birdnet_detections_clean.rds
birdnet_detections_clean.csv
acoustic_event_summary.rds
acoustic_event_summary.csv
monthly_summary.csv
multi_species_events.csv
species_month_counts_naive.csv
species_month_counts_hard_assignment.csv
species_month_counts_thresholds.csv
```

Verified cleaned dataset:

- valid `birdnet-lite` event-species candidates: **65,255**
- excluded anomaly-detection rows: **14,143**
- unique acoustic events: **63,985**
- predicted species labels: **229**
- represented calendar months: **23**
- cleaned source audio files: **35,067**
- date range: **21 December 2021 to 17 April 2024**
- site count: **1**
- recorder count: **1**

An acoustic event is reconstructed from source-audio identity and segment start / end times. Multiple predicted-species candidates can occur in the same acoustic event.

---

## 6.3 `data/02_manual_validation/` — manual-review design and decisions

### `round1_completed/`

Contains the completed Round-1 validation records.

Verified Round 1:

- **439** validation records
- **45** predicted species
- **437** unique acoustic events

Round 1 was designed to cover the BirdNET confidence-score range, common and less common predicted species, and both single-candidate and multi-candidate events.

### `round2_design/`

Contains the species sampling frame and Round-2 screening design.

Round-2 sampling targets were based on the total number of detections available for a predicted species:

```text
1-3 detections   -> all records
4-20             -> 4 records
21-100           -> 6 records
>100             -> 10 records
```

Low- and high-score strata were represented where possible.

### `round2_completed/`

Contains the completed Round-2 review workbook.

Verified Round 2:

- **708** validation records
- **184** species not represented in Round 1
- **705** unique acoustic events

### `final/`

Contains the frozen final combined manual-validation dataset and species summaries.

Final validation coverage:

- **1,147** validated event-species candidates
- **1,138** unique acoustic events
- all **229** predicted species labels represented
- **573** records classified as correct in the deliberately stratified validation sample
- unweighted confirmation proportion: **49.96%**

The 49.96% value describes the validation sample only. Because the two rounds had different sampling objectives and were not a simple random sample of the complete detection dataset, this proportion must **not** be interpreted as overall BirdNET accuracy for the full dataset.

Validation evidence basis:

- `audio_based`: **1,016** records
- `audio_plus_ecological_prior`: **131** records

The Primary analysis uses all 1,147 confident final judgements. The Audio-only sensitivity analysis uses only the 1,016 audio-based records.

---

## 6.4 `data/03_analysis_ready/` — model-ready generated inputs

### `validation/`

Generated by script 12:

```text
validation_analysis_ready.rds
validation_analysis_ready.csv
```

These files reconstruct model variables, including the final event-level multi-candidate indicator and the Primary / Audio-only analysis sets.

### `calibration/`

Generated by script 16:

```text
calibrated_candidates.rds
calibrated_candidates_minimal.csv
```

These files contain the candidate-level correctness probabilities used by the downstream diversity analyses.

---

## 6.5 `data/04_static_assets/`

Contains:

```text
Figure_01_validation_calibrated_framework.png
```

Figure 1 is a conceptual workflow diagram rather than a statistical output. It is therefore treated as a static manuscript asset. Script 27 copies it into the final manuscript-output directory.

---

# 7. Script-by-script analytical workflow

The default `run_all.R` pipeline executes the following computational chain.

## Step 1 — clean the raw BirdNET export

Script:

```text
scripts/01_data_preparation/01_clean_birdnet_data.R
```

Purpose:

- read the original 79,398-row export;
- separate `birdnet-lite` from anomaly-detection rows;
- retain valid predicted-species labels and BirdNET confidence scores;
- reconstruct event identifiers;
- create event-species candidate records;
- derive date / month and event-level summaries;
- generate deterministic candidate-count inputs used later for comparison;
- run fixed regression checks on row, event, species, and month counts.

Main outputs:

```text
data/01_cleaned/
results/00_data_preparation/
```

Expected checks:

```text
raw rows       79,398
clean rows     65,255
events         63,985
species        229
months         23
```

---

## Steps 2, 10, and 11 — reconstruct validation sampling designs

Scripts:

```text
scripts/02_manual_validation/02_generate_round1_validation_sample.R
scripts/02_manual_validation/10_build_species_sampling_frame.R
scripts/02_manual_validation/11_generate_round2_screening_sample.R
```

Purpose:

- reproduce the record-selection logic used for manual validation;
- document which species and score strata were targeted;
- preserve the sampling design for audit.

These steps are **not part of the default final run** because the completed human validation decisions are already frozen.

To rebuild the sampling-design outputs explicitly:

```bash
SILWOOD_REBUILD_VALIDATION_SAMPLING=true Rscript run_all.R
```

These scripts do not replace or overwrite the final manual decisions.

---

## Step 12 — prepare the frozen final validation dataset

Script:

```text
scripts/02_manual_validation/12_prepare_final_validation_data.R
```

Purpose:

- read the frozen combined manual-validation input;
- connect it to the cleaned event structure;
- check duplicates and event-field consistency;
- reconstruct analysis variables;
- define the Primary and Audio-only analysis sets;
- write model-ready validation data.

Verified quality checks:

```text
Total validation records             1,147
Primary analysis records             1,147
Audio-only sensitivity records       1,016
Ecological-prior records               131
Predicted species labels               229
Unique acoustic events               1,138
Round-1 records                        439
Round-2 records                        708
Duplicate detection IDs                  0
Duplicate validation IDs                 0
Event-field mismatches                   0
```

Outputs:

```text
data/03_analysis_ready/validation/
results/01_validation/
```

---

## Step 13 — fit calibration models

Script:

```text
scripts/03_calibration/13_fit_final_calibration_models.R
```

The binary response is whether a manually reviewed event-species candidate was correct.

BirdNET confidence scores are bounded to `[0.001, 0.999]` and logit transformed before model fitting.

Four nested Primary calibration models are fitted:

```text
M1  score-only model
M2  pooled event-aware model
M3  species-aware random-intercept model
M4  species-varying random-slope model
```

Conceptually:

```text
M1: correctness ~ logit(BirdNET confidence score)

M2: correctness ~ logit(BirdNET confidence score)
                  + multi-candidate event indicator

M3: M2
    + predicted-species random intercept

M4: M2
    + predicted-species random intercept
    + predicted-species random score slope
```

The Audio-only analysis refits M2 and M3 after excluding the 131 records whose validation basis included ecological / biogeographic prior information.

Verified Primary model fit:

```text
M1 AIC = 1342.03
M2 AIC = 1336.60
M3 AIC =  976.49
M4 AIC =  976.12
M3 species random-intercept SD = 2.327
M4 species score-slope SD      = 0.862
```

M4 has the lowest in-sample AIC, but final model selection is based on the prespecified event-grouped predictive-validation rule rather than AIC alone.

Outputs:

```text
results/02_calibration/models/
results/02_calibration/tables/
results/02_calibration/figures/
```

---

## Step 14 — repeated event-grouped predictive validation

Script:

```text
scripts/03_calibration/14_event_grouped_cross_validation.R
```

Purpose:

- evaluate probability prediction on held-out acoustic events;
- prevent candidates from the same acoustic event appearing in both training and test data;
- compare the four Primary calibration models on out-of-fold predictions.

Default final design:

```text
5 folds x 5 repeats = 25 fitted folds per model
```

Metrics:

- Brier score;
- binary log loss;
- calibration intercept;
- calibration slope;
- convergence / singularity diagnostics.

Verified event-grouped results:

| Model | Mean Brier | Mean log loss | Calibration intercept | Calibration slope |
|---|---:|---:|---:|---:|
| M1 score-only | 0.2010 | 0.5845 | 0.0004 | 0.9885 |
| M2 pooled event-aware | 0.1999 | 0.5821 | 0.0002 | 0.9828 |
| **M3 species-aware random intercept** | **0.1245** | **0.3954** | **0.0977** | **0.8398** |
| M4 species-varying random slope | 0.1258 | 0.4005 | 0.0929 | 0.8022 |

M3 had the lowest mean Brier score and log loss.

M4 did not satisfy the prespecified replacement rule because:

- it did not improve both predictive metrics relative to M3; and
- its maximum optimisation gradient across CV folds was approximately **113.73**, far above the stability threshold of 0.01.

Final downstream model:

```text
M3 species-aware random-intercept model
```

M4 is retained as a model-structure sensitivity analysis.

Outputs:

```text
results/03_model_validation/
```

---

## Optional Step 15 — leave-one-species-out diagnostic

Script:

```text
scripts/99_optional_diagnostics/15_leave_one_species_out_validation.R
```

This diagnostic is retained for transparency but is not part of the final thesis pipeline and is disabled by default.

To run it:

```bash
SILWOOD_RUN_LOSO=true Rscript run_all.R
```

It should be interpreted as an optional generalisation diagnostic rather than part of the final main model-selection workflow.

---

## Step 16 — assign calibrated probabilities to the full cleaned dataset

Script:

```text
scripts/03_calibration/16_apply_calibration_to_full_data.R
```

Purpose:

- apply the selected calibration models to all 65,255 cleaned event-species candidates;
- create point-estimate candidate-level correctness probabilities;
- apply pooled fallback probabilities where a species-specific effect is unavailable.

Verified probability assignment:

| Scenario | Candidates | Species | Mean probability | Expected correct candidate records | Pooled fallback |
|---|---:|---:|---:|---:|---|
| Primary species-aware hybrid | 65,255 | 229 | 0.7907 | 51,597.15 | 0 records / 0 species |
| Primary pooled | 65,255 | 229 | 0.5682 | 37,080.30 | all records pooled |
| Audio-only hybrid | 65,255 | 229 | 0.8005 | 52,237.99 | 1,919 records / 19 species |

The Primary M3 fit contains all 229 predicted species, so the Primary pooled fallback is not used in the actual Primary application.

The Audio-only species-aware fit contains 210 predicted species; the remaining 19 species use pooled Audio-only fallback probabilities.

These probabilities:

- are probabilities that individual event-species candidate labels are correct;
- do not sum to one within an acoustic event;
- are not occupancy probabilities;
- are not individual-presence probabilities;
- are not individual-bird abundance estimates.

Outputs:

```text
data/03_analysis_ready/calibration/
results/04_probability_assignment/
```

---

## Step 17 — Output 1: full Monte Carlo uncertainty propagation

Script:

```text
scripts/04_uncertainty_propagation/
17_monte_carlo_diversity_simulation_corrected.R
```

Default final run:

```text
1,000 simulations x 3 probability scenarios
```

Probability scenarios:

1. Primary species-aware;
2. Primary pooled;
3. Audio-only sensitivity.

For the Primary species-aware scenario, each simulation propagates three uncertainty layers:

1. **fixed-effect uncertainty**  
   fixed-effect coefficients are drawn jointly from their estimated multivariate normal distribution;

2. **conditional predicted-species random-intercept uncertainty**  
   species random intercepts are drawn approximately from their conditional estimate and conditional standard error;

3. **candidate-state uncertainty**  
   each event-species candidate is retained or removed with a Bernoulli draw based on its simulation-specific candidate-level correctness probability.

Candidates within an acoustic event are not forced to be mutually exclusive, because more than one species can occur in the same audio segment.

The simulation then calculates monthly Hill diversity:

```text
q = 0  species richness
q = 1  exponential Shannon diversity
q = 2  inverse Simpson diversity
```

and temporal:

```text
alpha diversity
multiplicative beta diversity
gamma diversity
```

under:

- equal-month weighting (primary); and
- abundance-weighted sensitivity analysis, where "abundance" means retained candidate-record totals, not individual birds.

Verified Monte Carlo quality checks:

```text
simulations per scenario                 1,000
monthly draw rows                       207,000
partition draw rows                      54,000
Primary count-array simulations           1,000
Primary count-array months                   23
Primary count-array species                 229
missing monthly diversity values              0
missing partition diversity values            0
```

Primary species-aware equal-month results (median [95% Monte Carlo simulation interval]):

| Diversity order | Alpha | Multiplicative beta | Gamma |
|---|---:|---:|---:|
| q = 0 richness | 49.130 [45.174, 53.565] | 3.219 [3.008, 3.421] | 158 [146, 171] |
| q = 1 exponential Shannon | 8.795 [8.089, 9.462] | 1.683 [1.666, 1.703] | 14.805 [13.551, 15.991] |
| q = 2 inverse Simpson | 4.711 [4.415, 5.005] | 1.733 [1.710, 1.761] | 8.174 [7.653, 8.710] |

Outputs:

```text
results/05_monte_carlo/
```

---

## Step 18 — compare alternative identification workflows

Script:

```text
scripts/05_method_comparison/18_compare_diversity_methods_final.R
```

Eight final analytical scenarios / workflows are compared:

```text
Primary species-aware
Primary pooled
Audio-only sensitivity
Naive: all candidates
Hard top-1
Threshold >= 0.50
Threshold >= 0.70
Threshold >= 0.90
```

Deterministic candidate retention:

| Deterministic workflow | Candidate records retained | Events represented | Species retained | Months |
|---|---:|---:|---:|---:|
| Naive | 65,255 | 63,985 | 229 | 23 |
| Hard top-1 | 63,985 | 63,985 | 226 | 23 |
| Threshold >= 0.50 | 55,285 | 54,416 | 217 | 23 |
| Threshold >= 0.70 | 29,142 | 29,022 | 148 | 23 |
| Threshold >= 0.90 | 11,935 | 11,934 | 65 | 23 |

The reference for percentage differences is the Primary species-aware scenario.

Selected verified comparisons:

- Primary pooled alpha diversity was **+19.73%**, **+27.20%**, and **+15.44%** relative to the Primary species-aware scenario for \(q=0,1,2\), respectively.
- Naive alpha diversity was **+60.97%**, **+61.11%**, and **+39.86%** for \(q=0,1,2\).
- Naive gamma diversity was **+44.94%**, **+65.30%**, and **+42.51%**.
- Naive richness beta was **-10.04%**, illustrating that multiplicative beta can move in a different direction from alpha and gamma.
- No single global BirdNET confidence-score threshold reproduced the Primary species-aware result across all diversity orders and partition components.

Complete values are stored in:

```text
results/06_method_comparison/tables/
```

Outputs:

```text
results/06_method_comparison/
```

---

## Step 19 — sampling-coverage proxies

Script:

```text
scripts/06_sampling_coverage/
19_sampling_effort_proxy_sensitivity_fixed5min.R
```

Coverage is reconstructed from the complete raw detection export.

Proxies:

- represented audio hours;
- active detection days;
- active detection hours.

Verified coverage reconstruction:

```text
raw detection rows                79,398
unique source audio IDs           36,759
months with represented audio         23
months with Primary diversity         23
low-coverage months                    9
```

Because one source recording represents five minutes, represented audio hours are reconstructed as:

```text
unique source audio IDs x 5/60
```

These hours refer only to recordings represented in the exported detections. Files with zero exported detections may be absent. Therefore represented audio hours are **not complete recorder uptime**.

Verified Spearman correlations with Primary species-aware monthly diversity:

| Diversity | Represented audio hours | Active detection hours | Active detection days |
|---|---:|---:|---:|
| q = 0 richness | **0.962** | **0.933** | **0.856** |
| q = 1 exponential Shannon | -0.275 | -0.291 | -0.219 |
| q = 2 inverse Simpson | -0.317 | -0.366 | -0.322 |

The correlations are diagnostic associations, not causal estimates of sampling effects.

Outputs:

```text
results/07_sampling_coverage/
```

---

## Step 20 — higher-coverage robustness analysis

Script:

```text
scripts/06_sampling_coverage/
20_low_coverage_month_sensitivity.R
```

Low-coverage months are defined from the lower quartile of represented audio hours or active detection days.

Verified:

```text
all represented months      23
low-coverage months           9
higher-coverage months       14
Primary simulations reused 1000
```

The same Primary Monte Carlo draws are repartitioned after excluding low-coverage months. This avoids adding a new source of Monte Carlo variation.

Across each deterministic workflow, the direction of the difference from the Primary species-aware scenario remained stable for all **9/9** combinations of:

```text
3 Hill orders x 3 partition components
```

Mean absolute change in the percentage difference after excluding low-coverage months:

```text
Hard top-1             2.42 percentage points
Naive                  2.72
Threshold >= 0.50      1.62
Threshold >= 0.70      1.23
Threshold >= 0.90      4.63
```

This supports robustness of the **relative method comparison**, while absolute temporal ecological interpretation remains limited by incomplete sampling coverage.

Outputs:

```text
results/08_coverage_robustness/
```

---

## Step 21 — Output 2: analytical benchmark

Script:

```text
scripts/07_analytical_benchmark/21_analytical_vs_monte_carlo.R
```

The analytical benchmark starts from the Primary species-aware point-estimate candidate-level correctness probabilities.

For species \(s\) in month \(j\), it calculates:

- expected correct candidate count;
- under conditional independence, probability that at least one candidate record for the species is correct;
- expected richness;
- expected relative candidate-record frequencies;
- plug-in \(q=1\) and \(q=2\) Hill numbers.

It compares:

1. analytical / plug-in estimates;
2. fixed-probability Monte Carlo, where the fitted probabilities are held constant and only candidate Bernoulli states are simulated;
3. full Monte Carlo, which also propagates calibration-parameter uncertainty.

Verified analytical-benchmark quality checks:

```text
species-month expectation rows       5,267
monthly analytical rows                 69
analytical partition rows               18
fixed-probability simulations        1,000
monthly comparison rows                 69
partition comparison rows               18
missing fixed monthly values              0
missing fixed partition values            0
```

Monthly analytical / plug-in agreement with fixed-probability Monte Carlo:

| Diversity order | Mean absolute percentage error | Maximum absolute percentage error | Proportion inside fixed-MC 95% interval |
|---|---:|---:|---:|
| q = 0 richness | **0.150%** | 0.505% | 1.00 |
| q = 1 exponential Shannon | **1.844%** | 5.737% | 1.00 |
| q = 2 inverse Simpson | **0.132%** | 0.980% | 1.00 |

Interpretation:

- analytical expected richness closely reproduces the fixed-probability Monte Carlo mean under the shared conditional-independence assumption;
- the inverse Simpson plug-in estimate is also very close in this dataset;
- exponential Shannon shows the largest nonlinear plug-in discrepancy;
- the analytical benchmark provides point estimates and does not replace the full Monte Carlo uncertainty intervals.

Outputs:

```text
results/09_analytical_benchmark/
```

---

# 8. Reporting scripts and exact thesis outputs

The core analytical scripts write detailed intermediate results. Scripts 22-27 then rebuild and collect the presentation-ready outputs used in the thesis.

## Step 22 — Figure 2

```text
scripts/08_reporting/22_rebuild_manuscript_figure2.R
```

Rebuilds the four-panel calibration / predictive-validation figure:

- multi-candidate calibration curves;
- single-candidate calibration curves;
- selected species random-intercept deviations;
- event-grouped out-of-fold Brier score and log loss.

## Step 23 — Figure 3

```text
scripts/08_reporting/23_rebuild_manuscript_figure3.R
```

Rebuilds monthly uncertainty-aware Hill diversity for \(q=0,1,2\).

## Step 24 — Figure 4

```text
scripts/08_reporting/24_rebuild_manuscript_figure4.R
```

Rebuilds the alternative-workflow comparison heatmap. Partition components are displayed consistently as:

```text
alpha -> multiplicative beta -> gamma
```

## Step 25 — Figure 5

```text
scripts/08_reporting/25_rebuild_manuscript_figure5.R
```

Rebuilds the analytical benchmark versus fixed-probability Monte Carlo figure.

## Step 26 — Supplementary outputs

```text
scripts/08_reporting/26_rebuild_supplementary_outputs.R
```

Rebuilds:

```text
Figure S1  all species-aware random intercepts
Figure S2  species-varying random-slope sensitivity
Figure S3  probability assignment
Figure S4  sampling-coverage correlations
Figure S5  higher-coverage robustness
Figure S6  analytical approximation error
```

and Supplementary Tables S1-S6.

## Step 27 — collect final outputs

```text
scripts/08_reporting/27_collect_final_outputs.R
```

Collects the exact thesis outputs into:

```text
results/10_final_outputs/
```

and writes:

```text
FINAL_OUTPUT_INDEX.csv
```

---

# 9. Where to find the exact thesis figures and tables

For a reader who does not want to inspect every intermediate result, start here:

```text
results/10_final_outputs/
```

## Main figures

```text
main_figures/
    Figure_01_validation_calibrated_framework.png
    Figure_02_calibration_and_validation.pdf
    Figure_02_calibration_and_validation.png
    Figure_03_output1_monthly_Hill_diversity.pdf
    Figure_03_output1_monthly_Hill_diversity.png
    Figure_04_workflow_difference_from_M3.pdf
    Figure_04_workflow_difference_from_M3.png
    Figure_05_output2_analytical_benchmark.pdf
    Figure_05_output2_analytical_benchmark.png
```

Figure 1 is the static conceptual workflow. Figures 2-5 are regenerated from saved analytical results by scripts 22-25.

## Main tables

```text
main_tables/
    Table_01_mathematical_notation.csv
    Table_02_calibration_model_selection.csv
    Table_03_primary_species_aware_diversity.csv
    Table_04_analytical_benchmark.csv
```

## Supplementary figures

```text
supplementary_figures/
    Figure_S01_M3_all_species_random_intercepts.pdf/.png
    Figure_S02_M4_random_slope_sensitivity.pdf/.png
    Figure_S03_probability_assignment.pdf/.png
    Figure_S04_sampling_coverage_correlations.pdf/.png
    Figure_S05_low_coverage_robustness.pdf/.png
    Figure_S06_analytical_approximation_error.pdf/.png
```

## Supplementary tables

```text
supplementary_tables/
    Table_S01_validation_by_species.csv
    Table_S02_calibration_fixed_effects.csv
    Table_S03_full_workflow_comparison.csv
    Table_S04_coverage_correlations.csv
    Table_S05_low_coverage_robustness.csv
    Table_S06_abundance_weighted_partitions.csv
```

Optional LaTeX table fragments are also written to:

```text
supplementary_latex/
```

The complete selected-output inventory is:

```text
results/10_final_outputs/FINAL_OUTPUT_INDEX.csv
```

---

# 10. Result folders in analytical order

Each `results/` folder contains the detailed, machine-readable intermediate outputs for one stage.

```text
results/00_data_preparation/
```

Cleaning summaries and quality checks.

```text
results/01_validation/
```

Validation summaries, species summaries, audit tables, and example records.

```text
results/02_calibration/
```

Fitted model objects, fixed effects, species effects, model-comparison tables, and calibration figures.

```text
results/03_model_validation/
```

Event-grouped CV predictions, repeat-level metrics, model summaries, fold diagnostics, selected-model decision, and CV figures.

```text
results/04_probability_assignment/
```

Full-data probability-application summaries and species-level probability summaries.

```text
results/05_monte_carlo/
```

Full Monte Carlo result object, monthly summaries, partition summaries, quality checks, and intermediate figures.

```text
results/06_method_comparison/
```

Monthly and partition-level comparisons among probabilistic and deterministic workflows.

```text
results/07_sampling_coverage/
```

Monthly coverage proxies, diversity-coverage joins, correlations, low-coverage-month classification, quality checks, and coverage figures.

```text
results/08_coverage_robustness/
```

Higher-coverage partitions, differences from the Primary species-aware scenario, robustness summaries, and figures.

```text
results/09_analytical_benchmark/
```

Species-month expectations, analytical monthly and partition estimates, fixed-probability Monte Carlo summaries, analytical-error summaries, and figures.

```text
results/10_final_outputs/
```

The exact subset selected for the thesis and Supplementary Material.

---

# 11. `run_all.R` controls

The default command:

```bash
Rscript run_all.R
```

runs the complete final pipeline from cleaning through final figure / table collection.

Optional environment variables:

```text
SILWOOD_MODE=quick|final
SILWOOD_N_SIMULATIONS=<integer>
SILWOOD_N_FOLDS=<integer>
SILWOOD_N_REPEATS=<integer>
SILWOOD_REBUILD_VALIDATION_SAMPLING=true|false
SILWOOD_RUN_LOSO=true|false
SILWOOD_BUILD_REPORTING=true|false
SILWOOD_START_AT=<step number>
SILWOOD_END_AT=<step number>
SILWOOD_CLEAN_RESULTS=true|false
```

Advanced path overrides supported by the shared configuration are:

```text
SILWOOD_OUTPUT_ROOT=<directory>
SILWOOD_PROCESSED_ROOT=<directory>
```

## Quick smoke test

A quick test uses fewer simulations and less cross-validation:

```bash
SILWOOD_MODE=quick Rscript run_all.R
```

Default quick-mode settings:

```text
20 Monte Carlo simulations
3 CV folds
1 CV repeat
```

Quick mode is intended only to test that the pipeline runs. It should **not** be used to reproduce the final thesis numerical results.

## Re-run only part of the pipeline

Example:

```bash
SILWOOD_START_AT=17 SILWOOD_END_AT=21 Rscript run_all.R
```

This runs steps 17-21 only.

Partial execution assumes that all prerequisite data and model objects from earlier steps already exist. `run_all.R` does not automatically infer and rebuild missing upstream dependencies.

## Delete generated results before a clean rebuild

```bash
SILWOOD_CLEAN_RESULTS=true Rscript run_all.R
```

This removes generated result folders before rebuilding when the run starts from step 1. Frozen raw and manual-validation inputs are not removed.

---

# 12. Logs and audit trail

Every `run_all.R` execution creates:

```text
log/run_YYYYMMDD_HHMMSS/
```

with:

```text
pipeline_log.txt
pipeline_manifest.csv
sessionInfo.txt
step_XX_*.log
```

`pipeline_log.txt` records the sequence and runtime of the steps.

`pipeline_manifest.csv` records:

- step number;
- step name;
- script path;
- per-step log file;
- start time;
- completion time;
- elapsed seconds;
- success / failure state;
- error message if a step fails.

Each `step_XX_*.log` contains the R messages generated by the corresponding script.

If a step fails, the pipeline stops immediately and reports the step-specific log file.

After a successful run, the master log, manifest, and `sessionInfo()` are copied to stable `latest_*` filenames.

---

# 13. Reference software environment

The verified full run used:

```text
R 4.3.3
Ubuntu 24.04.3 LTS
time zone: Europe/London
```

Principal attached package versions included:

```text
tidyverse 2.0.0
dplyr 1.1.4
tidyr 1.3.1
readr 2.1.5
ggplot2 4.0.1
lme4 2.0-6
Matrix 1.6-5
MASS 7.3-60.0.1
patchwork 1.3.2
lubridate 1.9.5
scales 1.4.0
```

The complete package and platform record is stored in:

```text
log/latest_sessionInfo.txt
```

Random-number generation is made reproducible through the central project seed:

```text
20260729
```

used by the shared pipeline configuration.

---

# 14. Quality-control checkpoints

The project contains explicit regression / quality checks so that a future run fails or becomes visibly inconsistent if the main analytical inputs change unexpectedly.

Core cleaning checks:

```text
79,398 raw export rows
65,255 cleaned event-species candidates
63,985 acoustic events
229 predicted species
23 represented months
```

Core validation checks:

```text
1,147 validation records
1,138 validation events
229 predicted species represented
1,016 audio-based records
131 audio-plus-ecological-prior records
0 duplicate detection IDs
0 duplicate validation IDs
0 event-field mismatches
```

Core Monte Carlo checks:

```text
1,000 simulations per scenario
207,000 monthly draw rows
54,000 partition draw rows
0 missing monthly diversity values
0 missing partition values
```

Core coverage checks:

```text
23 represented months
9 low-coverage months
14 higher-coverage months
1,000 Primary simulations reused
```

Core analytical-benchmark checks:

```text
5,267 species-month expectation rows
69 monthly analytical rows
18 analytical partition rows
1,000 fixed-probability simulations
0 missing fixed-probability monthly values
0 missing fixed-probability partition values
```

---

# 15. How the computational outputs map to the thesis questions

## Calibration question

Relevant code / outputs:

```text
scripts/03_calibration/
results/02_calibration/
results/03_model_validation/
Figure 2
Table 2
Figures S1-S3
Table S2
```

Main computational result: the species-aware random-intercept model substantially improved out-of-fold Brier score and log loss relative to pooled models.

## Output 1: propagate classification uncertainty into biodiversity estimates

Relevant code / outputs:

```text
scripts/04_uncertainty_propagation/
results/05_monte_carlo/
Figure 3
Table 3
Table S6
```

Main computational result: candidate-level classification and calibration uncertainty is propagated into monthly Hill diversity and temporal alpha, multiplicative beta, and gamma diversity with 95% Monte Carlo simulation intervals.

## Compare with conventional workflows

Relevant code / outputs:

```text
scripts/05_method_comparison/
results/06_method_comparison/
Figure 4
Table S3
```

Main computational result: naive retention, hard top-1 assignment, and global score thresholds systematically alter downstream diversity estimates; no single threshold reproduces the Primary species-aware result across all orders and partition components.

## Sampling-coverage robustness

Relevant code / outputs:

```text
scripts/06_sampling_coverage/
results/07_sampling_coverage/
results/08_coverage_robustness/
Figures S4-S5
Tables S4-S5
```

Main computational result: richness is strongly associated with represented sampling coverage, but exclusion of lower-coverage months does not reverse the direction of the principal deterministic-workflow comparisons.

## Output 2: analytical benchmark

Relevant code / outputs:

```text
scripts/07_analytical_benchmark/
results/09_analytical_benchmark/
Figure 5
Table 4
Figure S6
```

Main computational result: analytical expected richness and the \(q=2\) plug-in benchmark closely reproduce fixed-probability Monte Carlo means, while \(q=1\) shows a larger nonlinear discrepancy; full Monte Carlo remains necessary when calibration-parameter uncertainty and simulation intervals are required.

---

# 16. Interpretation cautions

## Event-species candidates are not individual birds

All species frequencies used by the diversity calculations are based on retained event-species candidate records. They are not direct counts of birds and must not be interpreted as individual-bird abundance.

## Classification uncertainty is not detection uncertainty

The calibration framework addresses the probability that an already-exported predicted label is correct.

It does not estimate:

- false-negative detection;
- occupancy;
- detection probability;
- individual abundance;
- species that were present but never exported as candidates.

## Candidate probabilities are not event-normalised multiclass probabilities

More than one candidate species can be correct within the same acoustic event. Candidate-level probabilities are therefore not constrained to sum to one within an event.

## Represented audio hours are not recorder uptime

Coverage proxies are reconstructed only from source recordings represented in the detection export. Zero-detection recordings may be absent.

## Temporal beta is not spatial beta

There is one recorder at one monitoring location. Months are treated as temporal communities. Multiplicative beta diversity therefore represents temporal differentiation among months, not spatial differentiation among sites.

## The analytical benchmark is not a replacement for Output 1

Expected richness has a direct analytical expression under the stated candidate-level conditional-independence assumption. The \(q=1\) and \(q=2\) calculations are plug-in Hill numbers based on expected candidate counts, not exact expectations of the nonlinear Hill numbers.

The analytical benchmark is therefore used to understand when direct calculation is an adequate point-estimate shortcut. Full Monte Carlo remains the main uncertainty-propagation framework.

---

# 17. Recommended reading order for an examiner or supervisor

For the fastest audit of the project:

1. **Read this README.**
2. Open:
   ```text
   results/10_final_outputs/
   ```
   to see the exact thesis figures and tables.
3. Inspect:
   ```text
   log/latest_pipeline_log.txt
   ```
   to confirm that the complete pipeline ran successfully.
4. Inspect:
   ```text
   results/03_model_validation/tables/event_cv_model_summary.csv
   ```
   for the calibration-model selection evidence.
5. Inspect:
   ```text
   results/05_monte_carlo/tables/
   ```
   for Output 1.
6. Inspect:
   ```text
   results/06_method_comparison/tables/
   ```
   for conventional-workflow comparisons.
7. Inspect:
   ```text
   results/07_sampling_coverage/
   results/08_coverage_robustness/
   ```
   for sampling-coverage diagnostics.
8. Inspect:
   ```text
   results/09_analytical_benchmark/
   ```
   for Output 2.
9. Follow the corresponding numbered script if a calculation needs to be traced back to code.

For a complete reconstruction, run:

```bash
Rscript run_all.R
```

---

# 18. File manifest

`FILE_MANIFEST.csv` provides an archive-level inventory and checksum record for the project files. It can be used to verify that files have not changed unintentionally between archive copies.

---

# 19. Final reproducibility statement

The computational part of the final thesis can be regenerated from the archived raw detector export, frozen manual-validation decisions, project scripts, and documented R environment.

The only non-computational stage is the expert manual review of audio and associated ecological / biogeographic evidence. Those judgements are preserved as frozen inputs rather than presented as something that can be recreated by code.

The archive is therefore designed to distinguish clearly between:

```text
original supplied data
-> reproducible cleaning
-> frozen human validation decisions
-> reproducible calibration
-> reproducible probability assignment
-> reproducible full Monte Carlo uncertainty propagation
-> reproducible Hill diversity / temporal partitioning
-> reproducible workflow comparison and robustness checks
-> reproducible analytical benchmark
-> reproducible final thesis figures and tables
```

This is the intended end-to-end audit trail for the project.
