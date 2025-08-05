# INFORM_DSP_MultiOmics_Manuscripts

This repository contains code used in the analysis for the manuscript:

**“Integrative Multi-Omics and Drug Sensitivity Profiling Reveals Potential Biomarkers and Therapeutic Strategies in Pediatric Solid Tumors”**  

## Overview

This project integrates ex vivo drug screening with multi-omics profiling across pediatric solid tumors. The repository includes the final drug response matrix (sDSS_asym) used in the analyses, the trained MOFA+ model for latent factor inference, and scripts for generating visualizations presented in the manuscript.

## Repository Structure

- `R/` – Analysis and visualization scripts used in the manuscript.
- `data/` – Manuscript-associated data including processed input matrices, supplementary tables.
- `model/` – Pre-trained MOFA+ model used for latent factor inference.
- `environment.yml` – Conda environment file listing required R packages and dependencies.
- `LICENSE` – GNU General Public License v3.0 (GPLv3).
- `README.md` – Project description and usage instructions.

## Requirements

Install required packages using the `environment.yml`:

```bash
conda env create -f environment.yml
conda activate inform-multiomics
