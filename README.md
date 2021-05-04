ROAM: Hands on session
=============

# Revisiting Abrahms et al 2017

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Twitter](https://img.shields.io/twitter/follow/CexyNature?style=social)](https://twitter.com/cexynature?lang=en)

This repository is an experimental sandbox which contains a collection of scripts to reproduce the manuscript of [Abrahms et al 2017](https://movementecologyjournal.biomedcentral.com/articles/10.1186/s40462-017-0104-2#Sec14). 

The R code from this repo comes from [Abrahms et al 2017](https://movementecologyjournal.biomedcentral.com/articles/10.1186/s40462-017-0104-2#Sec14) and from [Dana Paige Seidel MovEco-R-Workshop repository](https://github.com/dpseidel/MovEco-R-Workshop). Althought, some modifications were done to their code.

# Repository structure:

```
Project_name/
|--- README.md
|--- roam_syndromes.Rproj
|--- data/
	|--- README.md
	|--- simulated_groups/
        |--- group_CH1.Rds
        |--- group_CH2.Rds
        |--- ...
	|--- simulated_individuals/
        |--- central_CH1.Rds
        |--- central_CH2.Rds
        |--- ...
        |--- territorial_CH2.Rds
|--- figures/
  |--- simulated_individuals/
        |--- central_hist_CH1.png
        |--- central_reloc_CH1.png
        |--- ...
	|--- simulated_group_pca/
	      |--- pca_group_CH1.png
|--- 01_Abrahams_et_al_2017.R
|--- 02_calculate_metrics.R
|--- 03_do_syndromes_cluster.R


```

# Some required reading to understand analysis:

[Van Moorter et al 2015](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2656.12394)

[Dray et al 2010](https://esj-journals.onlinelibrary.wiley.com/doi/pdf/10.1007/s11284-010-0701-7)