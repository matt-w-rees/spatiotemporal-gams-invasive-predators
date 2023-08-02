
<!-- README.md is generated from README.Rmd. Please edit that file -->

## Overview

This repository contains the data and code used in the manuscript 'Spatial variation in predator diel activity patterns'. 

Here I’ve built generalised additive models (‘GAMs’), using the [‘mgcv’
R-package](https://cran.r-project.org/web/packages/mgcv/mgcv.pdf), of
invasive predator spatiotemporal activity in south-west Victoria,
Australia.

This paper aims to provide inference on the drivers of invasive predator
spatiotemporal activity in temperate, heterogeneous landscapes.
Particularly whether a subordinate mesopredator (feral cat *Felis
catus*) exhibits (signs of) behavioural avoidance of an apex predator
(Red fox *Vulpes vulpes*).

## Repository structure

**Raw data is contained in the “raw\_data/” folder**.

  - *“predator\_records.csv”* has a row for every cat and fox
    camera-trap detection (no time separation between consecutive
    detections).
  - *“camdata.csv”* contains camera-trap deployment information and
    explanatory variables. Note that there are plenty of camera-traps
    which were deployed but did not detect predators - predator\_records
    and camdata need to be merged to account for absences.

**R scripts are contained in the “r\_scripts/” folder**.

  - *“format\_records\_hour\_gams.R”* reformats the raw data files into
    a usable dataframe for spatiotemporal GAMs. This accounts for solar
    times using [average double anchorage
    methods](https://doi.org/10.1111/2041-210X.13290), removes
    non-indepedent detections (within 30 mins) and reshapes into a
    dataframe where, for each camera-trap, there are 0 - 23 rows for
    hour of the day, a column for cat and fox filled with the number of
    ‘independent’ detections in each hour, as well as columns for
    deployment information and explanatory variables. This dataframe is
    saved at **“derived\_data/predator\_counts\_hour.csv”**.
  - *“models.R”* uses mgcv to fit GAMs
  - scripts beginning with *“plot\_”* contain code used to create the
    figures. These scripts require the *“models”* script to be run first
    in the same session.


**Figures for the manuscript are contained in the “figs/” folder**.
