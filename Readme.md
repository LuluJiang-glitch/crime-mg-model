# Metric-Graph Models for UK Crime Data
This repository contains the R code used to analyse street-level crime data in the City of London using latent Gaussian spatial models defined on different spatial supports.

## Repository Structure

| File                 | Description                                                                                                                                                    |
| -------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `00_Pre.R`           | Prepares annual UK crime data and constructs the City of London street network from OpenStreetMap.                                                             |
| `01_Model_MG.R`      | Fits the metric-graph aggregated model.                                                                                                                        |
| `01_Model_Spatial.R` | Fits the planar spatial aggregated model.                                                                                                                      |
| `01_Model_Point_.R`  | Fits the planar point model.                                                                                                                        |
| `02_plot.R`          | Reads the fitted models and produces comparison figures and summaries.                                                                                         |
| `Function_Model.R`   | Contains supporting functions used for data preparation, graph operations, covariate extraction, distance calculation, and numerical integration.              |


Generated data objects and fitted models are stored locally in the following directories:
```text
RData/
├── Prepare/
├── MG/
├── Spatial/
└── Point/
```
These generated files are not included in the repository.


## Data
Street-level crime data are obtained from the UK Police data archive: https://data.police.uk/data/
The street network and amenity locations are obtained from OpenStreetMap using the `osmdata` R package.


## External Code
The file `metric_graph.R` used in this project was obtained from: https://github.com/karinalilleborge/LinesMetricGraph/blob/main/metric_graph.R

