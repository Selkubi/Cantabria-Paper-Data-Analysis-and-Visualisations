# Cantabria Paper Data Analysis and Visualizations

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13354231.svg)](https://doi.org/10.5281/zenodo.13354231)

## Overview

This repository contains the data analysis and visualization code accompanying our paper, "Riverine dissolved organic matter responds differently to alteration of two hydrological regimes from Northern Spain." The provided scripts facilitate the replication of analyses and figures presented in the publication.

## Installation

To utilize the analysis tools, install the package directly from this repository using the following R command:

```r
remotes::install_github("Selkubi/Cantabria_optical_Final", force = TRUE)
```
After installation, load the package with:

```r
library(getting.statistics)
```

## Usage
The analysis pipeline employed in our paper is encapsulated in the wrapper script located at analysis/00_analysis_pipeline_paper.R. This script orchestrates the data processing and visualization steps as described in the methods section of the publication.

## Example
To execute the analysis pipeline, run:

```r
source("analysis/00_analysis_pipeline_paper.R")
```
This script will reproduce the analyses and generate the plots  and results in the tables featured in the paper.

## Contributors
For inquiries or further information, please contact: kubilay.selin@gmail.com

## Citation

If you utilize this code or data in your research, please cite our paper:

> Kubilay, S. (2024). Riverine dissolved organic matter responds differently to alteration of two hydrological regimes from Northern Spain. *Journal Name*, *Volume*(Issue), pages. DOI: [10.5281/zenodo.13354231](https://doi.org/10.5281/zenodo.13354231)

## License
This project is licensed under the MIT License.

