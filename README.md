# Hybrid pilot designs for sparse functional data

This repository accompanies the paper **Pilot-Study Design for Functional Data Analysis: A Novel Hybrid Approach with Application to Mediterranean Fruit Fly Reproduction**, currently under review.

## Repository structure

``` text
paper-code/
├── README.md
├── LICENSE
├── .gitignore
│
├── R/
│   ├── hybrid_design.R
│   │   # Core Hybrid-design implementation and constructors
│   └── generate_hybrid_design.R
│       # Example for generating a hybrid design
│
└── paper-analysis/
    ├── README.md
    │
    ├── functions/
    │   ├── workflow_functions.R
    │   ├── evaluation_functions.R
    │   ├── plotting_functions.R
    │   └── pss_functions.R
    │
    ├── simulation-fourier/
    │   # Simulation
    │
    ├── case-study-medfly/
    │   # Medfly case-study
    │
    ├── supplementary-cluster-size/
    │   # Supplementary cluster-size analysis
    │
    └── data/
        ├── simulation-fourier/
        ├── case-study-medfly/
        └── supplementary-cluster-size/
            # Supplied sparse FPCA inputs and dense/reference data
```

-   **`R/`** contains the reusable code a researcher needs to generate the proposed Hybrid design for a new study;
-   **`paper-analysis/`** contains reference analyses from the publication.

## Paper analysis code:

``` text
supplied sparse datasets
        ↓
01_estimate_fpca_*.R
        ↓
02_evaluate_*.R
        ↓
03_plot_results.R
```

## R dependencies

The supplied code uses:

-   `lpSolve` for Hybrid-design construction;
-   `fdapace` for the paper's FPCA analyses;
-   `dplyr` and `ggplot2` for paper-analysis summaries and figures.

Generative AI was used to assist with reorganizing portions of the code for public release.

## Generate a Hybrid pilot design for your study

The reusable implementation is:

``` text
R/hybrid_design.R
```

The recommended interface is `generate_hybrid_design()`. Supply four design parameters:

-   `N`: number of subjects;
-   `gsize`: number of candidate observation locations on the common grid;
-   `num_obs`: number of observations collected per subject;
-   `W`: hybrid weight, interpreted here as the proportion of subjects assigned to the snippet component of the design.

An example is:

``` r
source("R/hybrid_design.R")

design <- generate_hybrid_design(
  N = 30,
  gsize = 25,
  num_obs = 5,
  W = 0.5,
  seed = 3429
)
```

When construction succeeds, `design` is an `N`-by-`gsize` binary incidence matrix:

-   each **row** represents a subject;
-   each **column** represents a candidate observation location;
-   an entry of `1` means that the location is observed for that subject;
-   each row contains `num_obs` selected locations.

For example, the selected grid positions for subject 1 are:

``` r
which(design[1, ] == 1)
```

If the candidate grid contains actual time points or other measurement coordinates, map the selected columns back to that grid:

``` r
gsize <- 25
grid <- seq(0, 1, length.out = gsize)
selected_times_subject1 <- grid[design[1, ] == 1]
```

### Choosing the Hybrid weight `W`

`W` controls the balance between the snippet and BIBD component of the Hybrid design. The paper setting uses:

``` r
W <- 0.5
```

The supplied implementation expects `0 < W <= 1`. Researchers can change `W` when a different balance is appropriate for their application.

### Ready-to-edit example

A small example script is provided at:

``` text
R/generate_hybrid_design.R
```

Edit the design parameters at the top of that file and run from the repository root:

``` bash
Rscript R/generate_hybrid_design.R
```
