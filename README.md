# :technologist: MAGs for the data analysis of hydrocarbon‑degrading bacteria. 

In this repository you can find the 42 MAGs used to test rbims and the userguide to replicate our results. The MAGs are submitted under project ID PRJNA816150. 
But you can download the data directly from the directories here displayed. 

- Raw Data: Contains the raw genomes and proteomes, to be functionally annotated within the rbims conda environment.
- Functional annotation: proteomes annotated with KEGG and InterProScan, ready to be used by rbims in R.
- Install: Instructions to generate annotations based on Kofamscan, InterproScan, dbCAN and MEROPS using the rbims conda environment: rbimsenv. Make sure to download only the Raw data.

# :Clipboard: Steps for replicating our analysis

1. Download the repository
2. Download the files from the "Functional annotation" directories and place them in the root of the “rbims_hydrocarbon-main” directory
3. Create an R project inside the root of the “rbims_hydrocarbon-main” directory
4. Run or knit the file “rbims hydroacarbon case study.rmd”
5. Check the "results" new directory and explore the figures and tables
