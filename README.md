# Time-varying Relationships Between Salmon Productivity and PDO and NPGO

## Purpose
The purpose of this analysis is identify to what extent the relationship between salmon productivity (recruits-per-spawner) and the **PDO** and **NPGO** indices has changed over time, and to what extent these changes are shared among regions. 


## Contents

### Directories

Name      | Description
----------|------------------------
/Figs     | Figures produced by specific analysis scripts.
/Data     | Input salmon and environmental data files.
/doc      | Pertinent documents and manuscript drafts.
/Output   | Model output files, mostly as .rds or .csv. Note: Some of the STAN .rds output files are **large**.
/R        | Compilation of `.R` and `.stan` files containing individual models, helper functions, and plotting functions.

### Analysis Scripts

Name                           | Description
-------------------------------|----------------------
Fit-DLM-Ricker.R               | Fits *individual* dynamic linear Ricker model to stock-recruitment data from each stock, with time-varying PDO and NPGO effects.
Fit-DLM-Ricker-plus.R          | Same as `Fit-DLM-Ricker.R` but with time-varying PDO and NPGO effects treated as an AR(1) process.
Fit-MARSS-Ricker.R             | Fits *stock-specific* DLM Ricker with MARSS package in R.
Fit-DLM-Region.R               | Fits *combined* DLM Ricker, where PDO and NPGO effects are common within regions (EBS, GOA, CA, WC).
Fit-DLM-Region-ARMA.R          | Same as `Fit-DLM-Region.R`, but with time-varying PDO and NPGO effects treated as ARMA(1) process.
Get-NPGO-PDO.R                 | Downloads up-to-date PDO and NPGO indices.



