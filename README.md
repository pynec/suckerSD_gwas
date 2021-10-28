# suckerSD_gwas
Scripts for conducing a GWAS using GEMMA 

This repository contains a few Java and R scripts for using GEMMA.

These files are still a work in progress but are working!

genotype_file.java is a java file that simulates random data for the three file inputs for GEMMA with no sex specific loci.

sal_sim.java is a java file that simulates random data for three file inputs for GEMMA with sex specific loci. The output of running these files on GEMMA can be used to plot a manhattan plot.

imputation.R is an R file that imputes missing genotypes using genotype likelihoods **As of 2021, use test_imputation.R as it uses allele information 
