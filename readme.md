Replication package for *Expectile-based Probabilistic Forecasting for
Spatio-Temporal River Network Data*
================
2025-10-29

This repository includes the code for expectile spatio-temporal
smoothing on river network data and forecasting the response observed on
the network.

-   Authors: Hyungryul Park, Joonpyo Kim, Seoncheol Park

-   Correspondence: [Joonpyo Kim](mailto:joonpyokim@sejong.ac.kr),
    [Seoncheol Park](mailto:pscstat@hanyang.ac.kr)

## Overview

This repositories consist of:

-   an R script (`Miho_code_sources.R`) containing the functions
    required for reproduction;

-   17 R scripts generating 16 figures and 1 table;

-   an RDS file (`Miho_forpost.RDS`) containing the dataset (Miho data)
    used in the analysis;

-   reproduced figures and tables stored in `plots` and `tables`
    folders, respectively; and

-   miscellaneous files necessary to reproduce Figure 1 and 7 (under
    `Figure_01_supp/`).

To reproduce the results, simply execute each R script. For instance, to
reproduce **Figure 9**:

``` r
source("Figure_09_estimatedcurves.R")
```

The resulting plot will be saved as `plots/Figure09.pdf`.

### Remark

” - The following two files should be executed **after** running
`Figure_12_forecastedexpectiles.R`:

-   `Figure_13_estimatedcdf.R`

-   `Figure_14_CRPSscores.R`

-   `Table_1_CRPSscores.R` should be executed **after** the 3 files
    above.

-   `Figure_16_estimatedcdf.R` should be executed **after**
    `Figure_15_forecastedexpectiles.R`.

-   Figure 1 and 7 require registering an API key for Stadia Maps.
    Please refer to the documentation:
    <https://search.r-project.org/CRAN/refmans/ggmap/html/register_stadiamaps.html>.
    After obtaining your API key, register it as follows:

``` r
library(ggmap)
register_stadiamaps(key = "YOUR-API-KEY") 
```

-   Approximate computation times (based on R 4.2.1 on a MacBook Pro
    with M1 Max and 64GB memory):

    -   Figure 12: about 5 minutes

    -   Figure 14: about 20 minutes

    -   Figure 15: about 10 minutes

    -   Remaining figures and the table are generated almost
        immediately, provided all dependencies are met.

## Dataset

The Miho River (or Miho-cheon) is the largest tributary of the Geum
River and is mainly located in North Chungcheong Province
(Chungcheongbuk-do) in South Korea, with a watershed area of 1,854km$^2$
and a total basin area of 9,912km$^2$ (Yu et al., 2024). Several water
quality indicators are measured at these stations including total
nitrogen (TN). In this study, we focus on total nitrogen levels observed
from January 3, 2008, to October 30, 2024, across 28 stations, resulting
in a total of 11,060 spatio-temporal data points. The dataset is
publicly available through the Water Environment Information System of
South Korea (<https://water.nier.go.kr/>), or in this repository
(`Miho_forpost.RDS`).
