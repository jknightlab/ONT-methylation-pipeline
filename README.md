# ONT methylation data processing pipeline

Author:		Kiki Cano-Gamez

Email:		kiki.canogamez@well.ox.ac.uk


## Overview

This repository contains a collection of codes to process sequencing data and modified based calling generated using Oxford Nanopore (ONT) sequencing.


## Repository structure

The codes contained within this repository correspond to the main data processing steps followed. They are written in bash and ordered as follows:

```
./
 |-- 0_merge-base-calls.sh			Merges all base calls giles (in BAM format) generated per barcode using the MinKNOW software into a single BAM file.
 `-- 1_trim-adapters.sh				Trims adapter and native barcode sequences from ONT reads using Dorado.
```

