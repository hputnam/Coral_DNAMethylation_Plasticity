Coral_DNAMethylation_Plasticity

This repository includes data and analysis scripts to accompany: Ocean acidification influences host DNA methylation and phenotypic plasticity in environmentally susceptible corals

Authors: Putnam, Hollie M., Davidson, Jennifer M., Gates, Ruth D.
Journal: Evolutionary Applications
Link: 
Description: This repository provides data and scripts to analyze the influence of ocean acidification on host DNA methylation and phenotypic plasticity. Clonal coral fragments from Montipora capitata and Pocillopora damicornis were exposed to ambient and low fluctuating pH treatments for 6 weeks. Coral growth was measured every 2 weeks and samples were taken after 6 weeks of exposure for analysis of host DNA methylation and HNMR metabolomic profiling. 

Contents: RAnalysis file containing three folders (Scripts, Data, Output)

Scripts:
	Coral_DNA_Methylation.R: A R script that imports all data and runs all analyses for the maniscript and generates all figures and tables
	opls.R: A R script required for the metabolomic profiling analysis (orthogonal partial least squares discriminate analysis). OPLS DA analysis script by Paul Anderson published in Sogin et al 2014 (http://birg.cs.cofc.edu/index.php/O-PLS)

Data:
	BM_Light_Calibration_Data.csv
	BM_Acclimation.csv
	BM_Field_Temp.csv
	BM_Tank_Temp.csv
	BM_Tank_light.csv
	BM_Daily_Measurements.csv
	BM_SWChem.csv
	BM_NBS_pH.csv
	BM_NMR_Data_0.04binned_truncated_nowater.csv
	BM_NMR_sample_info
	pH probe calibration files : Coral_DNAMethylation_Plasticity/R_Analysis/Data/pH_Calibration_Files/
	BM_Methylation.csv
	BM_Buoyant_Weight.csv

Output: 
	The directory containing the analysis results, figures, tables, and supplementary material.

