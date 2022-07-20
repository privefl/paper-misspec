## Structure of the code

### Prepare data

- `prepare-genotypes-simu-chr22.R`: Prepare genotypes for simulations using chromosome 22 + figure S1

- `prepare-corr-simu-chr22.R`: Prepare the LD reference used in simulations + figure S3

- `prepare-corr-simu-altpop.R`: Prepare the "alternative" LD reference used in simulations

- `prepare-genotypes.R`: Prepare genotypes for real data applications

- `prepare-corr.R`: Prepare the main LD reference used in real data analyses

- `prepare-phecodes.R`: Prepare phecodes (+ a few continuous traits) for real data analyses in the UK Biobank

- `prepare-sumstats/*`: prepare GWAS summary statistics for the main real data analyses + figures 4 & S17-S20 & S23-S27 & S33-S34

- `prepare-sumstats-finngen/*` + `prepare-sumstats-bbj/*`: prepare GWAS summary statistics from FinnGen and Biobank Japan

- `prepare-genotypes-1000G.R`: Prepare genotypes from the 1000 Genomes data for variants in common with the HapMap3 ones used in the UKBB, to be used to compute other LD references

- `prepare-corr-1000G-EUR.R` + `prepare-corr-FIN.R`: prepare two of the three LD references used with the FinnGen GWAS summary statistics (the last one is the UKBB one from the validation set)

- `prepare-corr-1000G-EAS.R` + `prepare-corr-JPN.R` + `prepare-corr-UKBB-EAS.R`: prepare the three LD references used with the BBJ GWAS summary statistics

- `prepare-genotypes-EAS.R`: Prepare genotypes for testing prediction in the UK Biobank East Asian individuals


### Simulations results

- `run-methods.R`: code to run all methods in simulations

- `simu-misspec-N.R`: figures 1 & S6-S7

- `investigate-misspec-N.R`: figures S4-S5 & S35

- `simu-dosage-info.R`: figure 2 & S14

- `simu-altpop.R`: figure 3 & S15-S16

- `simu-rounded.R`: figure S39 (an example of such real GWAS summary statistics is given in `investigate-T2D-rounded.R`)

- `investigate-GWAS-from-imputed.R`: figures S8-S12

- `investigate-PCwithLD-GWAS.R`: figure S21

- `investigate-lassosum2-auto.R`: figure S2

- `investigate-corr-imp.R`: figure S38


### Real data results

- `run-PRS.R` + `process-PRS.R`: run all main real data analyses + figures 5-6 & S22 & S28-S32

- `investigate-CAD.R`: figure S36

- `run-PRS-FIN.R` + `process-PRS-FIN.R`: run real data analyses using FinnGen GWAS summary statistics + figure 7

- `run-PRS-BBJ.R` + `process-PRS-BBJ.R`: run real data analyses using Biobank Japan GWAS summary statistics + figure S37


### Export data

- `prepare-ldref-with-blocks.R`: Add blocks to LD reference matrices provided in the LDpred2 paper

- `compute-all-info-eur.R`: Recompute INFO scores (and MAF) for the NW European subset in the UK Biobank + figure S13
