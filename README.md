# ONT methylation data processing pipeline

Author:		Kiki Cano-Gamez

Email:		kiki.canogamez@well.ox.ac.uk


## Overview

This repository contains a collection of codes to process sequencing data and modified based calling generated using Oxford Nanopore (ONT) sequencing.


## Repository structure

The codes contained within this repository correspond to the main data processing steps followed. They are written in bash and ordered as follows:

```
./
 |-- 0_merge-base-calls.sh			Merges all base calls files (in BAM format) generated per barcode by the MinKNOW software into a single BAM file.
 |-- 1_trim-adapters.sh				Trims adapter and native barcode sequences from ONT reads using porechop.
 |-- 2_align-reads.sh				Aligns trimmed ONT reads to the genome using minimap2
 |-- 3_repair-MM-rags.sh			Repairs methylation information encoded in the MM tags of the aligned BAM file using modkit.
 |-- 4_clip-reads.sh				Removes methylation information from the end of reads (i.e. read clipping) using modkit.
 |-- 5_get-methylation-calls.sh			Recovers pileup at CpGs and quantifies 5mC and 5hmC events at each position of interest using modkit.
 |-- 6_deconvolute-with-nanomix.sh		Performs tissue of origin deconvolution on cfDNA methyltomes using the Nanomix software.
 `-- get-mapping-statistics.sh			Computes general mapping and alignment statistics from an aligned BAM file.
```

