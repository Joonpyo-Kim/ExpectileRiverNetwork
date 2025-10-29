Replication package for *Expectile-based Probabilistic Forecasting for
Spatio-Temporal River Network Data*
================
2025-10-28

This repository includes the code for expectile spatio-temporal
smoothing on river network data and forecasting the response observed on
the network.

-   Authors: Hyungryul Park, Joonpyo Kim, Seoncheol Park

-   Correspondence: [Joonpyo Kim](mailto:joonpyokim@sejong.ac.kr),
    [Seoncheol Park](mailto:pscstat@hanyang.ac.kr)

## Overview

This repositories consist of:

-   an R file (`Miho_code_sources.R`) of functions required for
    reproduce;

-   17 R files generating 16 figures and 1 table;

-   an RDS file (`Miho_forpost.RDS`) including the dataset (Miho data)
    used;

-   reproduced figures and tables, under `plots` and `tables` folders,
    respectively;

-   and miscellaneous files necessary to reproduce Figure 1 and 7 (under
    `Figure_01_supp/`).

To reproduce results, one can just execute each R file. For instance, to
reproduce Figure 9,

``` r
source("Figure_09_estimatedcurves.R")
```

The resulting plot is saved as `plots/Figure09.pdf`.

### Remark

-   Following two files should be executed followed by
    `Figure_12_forecastedexpectiles.R`.

    -   `Figure_13_estimatedcdf.R`

    -   `Figure_14_CRPSscores.R`

-   `Table_1_CRPSscores.R` should be executed followed by above 3 files.

-   `Figure_16_estimatedcdf.R` should be executed followed by
    `Figure_15_forecastedexpectiles.R`.

-   Figure 12 requires about 5 minutes to be reproduced. Figure 14
    spends about 20 minutes. Figure 15 is reproduced in 10 minutes.
    Remained figures and table are generated almost immediately,
    provided that all requirements are met.

-   The authors used R 4.2.1 in Macbook Pro with M1 Max, 64GB Memory.

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
