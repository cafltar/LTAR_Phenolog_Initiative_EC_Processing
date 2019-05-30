# Phenology EC Processing 

## Introduction

This repo contains the code used to process the eddy covariance data acquired for the LTAR-Network synthesis around the PhenoCam and EC data sets for 2017-2018 timeframe. A brief description of each script is below. All scripts written in Python V3. This READ_ME is still a work in progress and will be updated over time.

## Script Summaries

### JORN_Combine

This script combines the data from the different NEON datastreams after the required columns have been pulled using the LTAR_NEON_Data_Extraction and NEON_Extra_Data_Pull scripts.

### UCB_Format

Forms the time information for the UCB site data into a quasi-ISO8601 format for use within Python.

### Reddy_Format

Function to format data to match with the input format needed for REddyProc R code.

### NEON_Extra_Data_Pull

Script to concat data from multiple sets of NEON *.csv files and truncate columns to specified headers; precursor to JORN_Combine.

### LTAR_NEON_Data_Extraction

Script to extra specific data from the *.h5 flux data files downloaded from the NEON Data Portal.

### LTAR_Pheno_QC_Functions

### LTAR_Phenology_Final_Processing

### LTAR_Phenology_Initiative_QC_REddy

## Contacts

Eric Russell: eric.s.russell@wsu.edu
