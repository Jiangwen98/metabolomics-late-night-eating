# Late-night eating: circadian postprandial metabolomics

Semi-targeted LC-MS/MS plasma metabolomics of the **Honours** cohort, from my PhD thesis
*Metabolomic phenotyping of different nutritional factors: focusing on meal timing and dietary
composition* (Department of Medicine, National University of Singapore, defended June 2026).

**Study design.** Three-way intra-subject crossover, an identical low-glycaemic-index meal consumed at 08:00, 20:00 and 00:00. The trial was designed and run by
collaborators, who provided plasma samples and de-identified clinical data for metabolomic
analysis.

This repository documents **methods and code**. It contains **no participant data**.

## Layout

```
analysis/    my analysis scripts for this cohort, plus shared_functions_Honours.R,
             the utility functions they call
data/        synthetic quant table documenting the expected input format
data-raw/    script that generates the synthetic table
```

## Assay and input

Agilent 1290 Infinity LC coupled to a 6495C triple quadrupole, SeQuant ZIC-cHILIC column,
dynamic MRM in positive and negative mode, roughly 300 to 500 analyte transitions spanning amino
acids, carnitines and acylcarnitines, nucleotides, organic acids and other polar metabolites.
Peak integration in [MRMkit](https://github.com/HuiSongLab/MRMkit).

MRMkit emits a tab-delimited quant table with three information columns, ten metric rows above
the sample rows (row 7 CoV, row 9 D-ratio, row 10 S/N), and pooled QC samples marked `BQC`.
`data/synthetic_quant_table.txt` is a synthetic table in that format, containing no real
measurements, so the reshaping logic can be checked without data.

## Quality control

| Metric | Threshold | Computed on |
|---|---|---|
| Signal-to-noise (S/N) | `>= 5` | uncorrected data |
| Coefficient of variation (CoV) | `<= 30%` | batch-corrected data |
| Dispersion ratio (D-ratio) | `<= 50%` | batch-corrected data |

An analyte is retained if `S/N >= 5 AND (CoV <= 30% OR D-ratio <= 50%)`. The CoV/D-ratio
criterion is deliberately a disjunction: CoV measures technical dispersion in absolute terms
while D-ratio expresses non-biological variance relative to biological variance, so an analyte
can be usable by one criterion and not the other and requiring both would discard signal.

Every surviving analyte is then inspected manually on four criteria, and analytes failing
manual review are excluded whether or not they passed the automated gate: peak shape, baseline
stability, retention-time consistency across samples, and agreement between subject samples and
pooled QC.

Intra-batch correction uses Gaussian kernel regression in MRMkit. Because the assay alternates
two HPLC columns to reduce cycle time, each column is treated as its own batch, so the batch
count passed to MRMkit is double the number of runs.



## Response measure

The response measure is baseline-normalised area under the curve (nAUC), each postprandial timepoint divided by that participant's own fasting baseline before integration.

Throughout, the preprandial fasting metabolome is analysed separately from the postprandial or
post-intervention response. These answer different questions, and a difference present at
baseline is not evidence of a different response.

## Statistics

```
Baseline ~ Time + (1 | Subject)
nAUC     ~ Time + (1 | Subject)
```

where `Time` is meal timing, and `(1 | Subject)` adjusts for repeated measures within participant.

Linear mixed models with a participant random effect (`lme4`, `lmerTest`), pairwise contrasts across the three timepoints via `emmeans`, Benjamini-Hochberg FDR.

**Small-sample handling.** This cohort has nine participants, which materially limits power over a panel of roughly 300 analytes. For the nAUC analysis both raw and BH-adjusted p-values are reported, ranking is done on raw p-values with that limitation stated explicitly, and every candidate signal is cross-checked against its time-course profile plot rather than being accepted on a p-value alone. The baseline analysis is restricted to BH-significant analytes. Findings from this cohort are exploratory pending replication.

## Provenance and attribution

I analysed this cohort independently in R. The utility functions my scripts call are in `analysis/shared_functions_Honours.R`. They **originate in earlier
analysis code from the laboratory**, inherited with the B&B project, and are not my own work. Where I
modified a function for this analysis its header is tagged "adapted from earlier lab code"; the rest
are as inherited.

See [ATTRIBUTION.md](ATTRIBUTION.md) for detail.

## Scope and limits

**No participant data is included.** Trial participant identifiers that appeared in the original
scripts have been redacted and replaced with `<..._REDACTED...>` placeholders, and absolute local
paths made relative. Affected files carry a header noting this.

**Results are not reported here.** Findings from this cohort are unpublished and belong to collaborative work in progress, so this repository describes how the analysis was done and not what it found.

**The scripts will not run end to end without the input data**, which is not distributed. The
shared laboratory functions they call are included as `analysis/shared_functions_Honours.R`; source that file
first, since in the original working environment these functions were loaded into the R session
directly rather than via `source()`.

## Dependencies

R, with `tidyverse`, `lme4`, `lmerTest`, `emmeans`, `qvalue` (Bioconductor), `mixOmics`,
`ggplot2`, `ggrepel`, `ggpubr`, `rstatix`, `car`, `gridExtra`. Peak integration and intra-batch
correction are performed externally in MRMkit; acquisition used Agilent MassHunter Workstation
10.0.127.

## Acknowledgements

The parent trial was **designed and conducted by Dr Gloria Leung and Prof. Maxine Bonham** at
Monash University, who provided the plasma samples and de-identified clinical data analysed here.
Ethics approval was granted by the Monash University Human Research Ethics Committee (project code
CF15/1301-2015000620) and the trial is registered with the Australian New Zealand Clinical Trials
Registry (ACTRN12616000164493).

I did not design or conduct the trial, recruit participants, or collect samples. My contribution
is the metabolomic analysis and the computational work in this repository.

## A note on tooling

The documentation in this repository was drafted with AI assistance: this README,
[ATTRIBUTION.md](ATTRIBUTION.md), and the synthetic data generator in `data-raw/`.

The analysis code in `analysis/` was not. It is my own work, other than the laboratory's
pre-established utility functions identified in [ATTRIBUTION.md](ATTRIBUTION.md).

## Licence

Code I authored is released under the MIT Licence, see [LICENSE](LICENSE). `analysis/shared_functions_Honours.R` originates in earlier laboratory code and is not my own work.

## Contact

Dong Jiangwen. [Google Scholar](https://scholar.google.com/citations?user=NE3oTFMAAAAJ) ·
[LinkedIn](https://www.linkedin.com/in/jiangwen-dong)
