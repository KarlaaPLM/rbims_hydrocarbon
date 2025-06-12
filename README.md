# :technologist: MAGs for the data analysis of hydrocarbon‑degrading bacteria. 

In this repository you can find the 42 MAGs used to test rbims and the userguide to replicate our results. The MAGs are submitted under project ID PRJNA816150. 
But you can download the data directly from the directories here displayed. 

- Raw Data: Contains the raw genomes and proteomes, to be functionally annotated within the rbims conda environment.
- Functional annotation: proteomes annotated with KEGG and InterProScan, ready to be used by rbims in R.
- Install: Instructions to generate annotations based on Kofamscan, InterproScan, dbCAN and MEROPS using the rbims conda environment: rbimsenv. Make sure to download only the Raw data.


# :computer: Replicate our data analysis 
## Install packages
```r
install.packages("devtools")
library(devtools)
install.packages("tidyverse")
library(tidyverse)
devtools::install_github("mirnavazquez/RbiMs", force = T)
library(rbims)
library(readxl)
library(dplyr)
```
## Read data and metadata
```r
metadata <- read_xlsx("metadata_SIPH.xlsx")
interpro_pfam_profile<-read_interpro(data_interpro = "Interpro_Hidro", 
                                     database="Pfam", profile = T)
interpro_pfam_profile_renamed <- interpro_pfam_profile %>%
  rename(
  "g_Flavobacterium_5m_16" =  "5mSIPHEX2_16",
  "g_Flavobacterium_5m_26" =  "5mSIPHEX1_26",
  "g_Henriciella_5m_15" =  "5mSIPHEX1_15",
  "g_Hyphomonas_5m_32" =  "5mSIPHEX1_32",
  "g_Hyphomonas_5m_33" =  "5mSIPHEX1_33",
  "g_Celeribacter_5m_10" =  "5mSIPHEX2_10",
  "g_Celeribacter_5m_0" =  "5mSIPHEX1_0",
  "s_Planktomarina_temperata_5m_1" =  "5mSIPHEX1_1",
  "s_Lentibacter_algarum_5m_13" =  "5mSIPHEX1_13",
  "s_Lentibacter_algarum_5m_7" =  "5mSIPHEX2_7",
  "g_Tateyamaria_5m_25" =  "5mSIPHEX2_25",
  "g_Tateyamaria_5m_8" =  "5mSIPHEX1_8",
  "o_Pseudomonadales_5m_2" =  "5mSIPHEX1_2",
  "s_Thalassolituus_oleivorans_5m_5" =  "5mSIPHEX2_5",
  "s_Thalassolituus_oleivorans_5m_19" =  "5mSIPHEX1_19",
  "g_Pseudophaeobacter_5m_3" =  "5mSIPHEX2_3",
  "g_Pseudophaeobacter_5m_37" =  "5mSIPHEX1_37",
  "g_Pseudophaeobacter_700m_8" =  "700mSIPHEX1_8",
  "g_Pseudophaeobacter_700m_13" =  "700mSIPHEX2_13",
  "g_Glaciecola_5m_9" =  "5mSIPHEX1_9",
  "g_Glaciecola_700m_16" =  "700mSIPHEX2_16",
  "g_Glaciecola_700m_18" =  "700mSIPHEX1_18",
  "s_Alcanivorax_jadensis_5m_11" =  "5mSIPHEX1_11",
  "g_Alcanivorax_700m_20" =  "700mSIPHEX1_20",
  "g_Alcanivorax_5m_25" =  "5mSIPHEX1_25",
  "g_Alcanivorax_5m_18" =  "5mSIPHEX2_18",
  "s_Marinobacter_salarius_5m_14" =  "5mSIPHEX2_14",
  "s_Marinobacter_salarius_700m_3" =  "700mSIPHEX1_3",
  "s_Marinobacter_salarius_5m_10" =  "5mSIPHEX1_10",
  "s_Marinobacter_salarius_700m_24" =  "700mSIPHEX2_24",
  "g_Oleibacter_5m_18" =  "5mSIPHEX1_18",
  "g_Oleibacter_700m_21" =  "700mSIPHEX2_21",
  "g_Oleibacter_700m_15" =  "700mSIPHEX1_15",
  "g_Olleya_700m_17" =  "700mSIPHEX1_17",
  "g_Olleya_700m_14" =  "700mSIPHEX2_14",
  "g_Dokdonia_700m_23" =  "700mSIPHEX2_23",
  "g_Dokdonia_700m_12" =  "700mSIPHEX1_12",
  "g_Paracoccus_700m_9" =  "700mSIPHEX2_9",
  "g_Paracoccus_700m_1" =  "700mSIPHEX1_1",
  "g_Sulfitobacter_700m_0" =  "700mSIPHEX1_0",
  "g_Alteromonas_700m_22" =  "700mSIPHEX2_22",
  "g_Alteromonas_700m_2" =  "700mSIPHEX1_2"
)
```
## Subset important features
```r
important_PFAMs<-get_subset_pca(tibble_rbims=interpro_pfam_profile_renamed, 
                                    cos2_val=0.97,
                                    analysis="Pfam")
```
## Plot with the visualization tools
### Lets order the taxa names and filter per environment. 
```r
important_PFAMS_5m <- important_PFAMs %>%
  select("Pfam","domain_name",
         "g_Flavobacterium_5m_16",
         "g_Flavobacterium_5m_26",
         "g_Henriciella_5m_15",
         "g_Hyphomonas_5m_32",
         "g_Hyphomonas_5m_33",
         "g_Celeribacter_5m_10",
         "g_Celeribacter_5m_0",
         "s_Planktomarina_temperata_5m_1",
         "s_Lentibacter_algarum_5m_13",
         "s_Lentibacter_algarum_5m_7",
         "g_Tateyamaria_5m_25",
         "g_Tateyamaria_5m_8",
         "o_Pseudomonadales_5m_2",
         "s_Thalassolituus_oleivorans_5m_5",
         "s_Thalassolituus_oleivorans_5m_19",
         "g_Pseudophaeobacter_5m_3",
         "g_Pseudophaeobacter_5m_37",
         "g_Glaciecola_5m_9",
         "g_Alcanivorax_5m_25",
         "g_Alcanivorax_5m_18",
         "s_Marinobacter_salarius_5m_14",
         "s_Marinobacter_salarius_5m_10",
         "g_Oleibacter_5m_18",
         "s_Alcanivorax_jadensis_5m_11")

important_PFAMs_700m <- important_PFAMs %>%
  select("Pfam","domain_name", 
         "g_Pseudophaeobacter_700m_8",
         "g_Pseudophaeobacter_700m_13",
         "g_Glaciecola_700m_16",
         "g_Glaciecola_700m_18",
         "g_Alcanivorax_700m_20",
         "s_Marinobacter_salarius_700m_3",
         "s_Marinobacter_salarius_700m_24",
         "g_Oleibacter_700m_21",
         "g_Oleibacter_700m_15",
         "g_Olleya_700m_17",
         "g_Olleya_700m_14",
         "g_Dokdonia_700m_23",
         "g_Dokdonia_700m_12",
         "g_Paracoccus_700m_9",
         "g_Paracoccus_700m_1",
         "g_Sulfitobacter_700m_0",
         "g_Alteromonas_700m_22",
         "g_Alteromonas_700m_2")

  order_taxa_5m <- c(
          "g_Celeribacter_5m_0",
          "g_Celeribacter_5m_10",
          "g_Henriciella_5m_15",
          "g_Hyphomonas_5m_32",
          "g_Hyphomonas_5m_33",
          "g_Pseudophaeobacter_5m_3",
          "g_Pseudophaeobacter_5m_37",
          "g_Tateyamaria_5m_25",
          "g_Tateyamaria_5m_8",
          "s_Lentibacter_algarum_5m_13",
          "s_Lentibacter_algarum_5m_7",
          "s_Marinobacter_salarius_5m_10",
          "s_Marinobacter_salarius_5m_14",
          "s_Planktomarina_temperata_5m_1",
          "g_Flavobacterium_5m_16",
          "g_Flavobacterium_5m_26",
          "g_Alcanivorax_5m_18",
          "g_Alcanivorax_5m_25",
          "g_Glaciecola_5m_9",
          "g_Oleibacter_5m_18",
          "o_Pseudomonadales_5m_2",
          "s_Alcanivorax_jadensis_5m_11",
          "s_Thalassolituus_oleivorans_5m_19",
          "s_Thalassolituus_oleivorans_5m_5"
          )
order_taxa_700m <- c(  
          "g_Paracoccus_700m_1",
          "g_Paracoccus_700m_9",
          "g_Pseudophaeobacter_700m_13",
          "g_Pseudophaeobacter_700m_8",
          "g_Sulfitobacter_700m_0",
          "s_Marinobacter_salarius_700m_24",
          "s_Marinobacter_salarius_700m_3",
          "g_Dokdonia_700m_12",
          "g_Dokdonia_700m_23",
          "g_Olleya_700m_14",
          "g_Olleya_700m_17",
          "g_Alcanivorax_700m_20",
          "g_Alteromonas_700m_2",
          "g_Alteromonas_700m_22",
          "g_Glaciecola_700m_16",
          "g_Glaciecola_700m_18",
          "g_Oleibacter_700m_15",
          "g_Oleibacter_700m_21")
```
#### Bubble plot

```r
metadata_5m_renamed <- read_xlsx("5m_renamed.xlsx")
bubble_5m <- 
            plot_bubble(important_PFAMS_5m, 
            y_axis=Pfam, 
            x_axis=Bin_name, 
            calc = "Binary",
            range_size = 3, 
            text_x = 9,
            text_y = 10.5,
            analysis = "INTERPRO", 
            order_bins = order_taxa_5m,
            data_experiment = metadata_5m_renamed,
            color_character = Class)

metadata_700m_renamed <- read_xlsx("Hidrocaruros/700m_renamed.xlsx")
bubble_700m <- 
            plot_bubble(important_PFAMs_700m, 
            y_axis=Pfam, 
            x_axis=Bin_name, 
            range_size = 3,
            calc = "Binary",
            text_x = 9,
            text_y = 10.5,
            analysis = "INTERPRO", 
            order_bins = order_taxa_700m,
            data_experiment = metadata_700m_renamed,
            color_character = Class) 
```
#### Heatmap plot
If we set the distance option as TRUE, we can plot to show how the samples could cluster based on the protein domains.

```r
plot_heatmap(important_PFAMs, 
             y_axis=Pfam, 
             analysis = "INTERPRO", 
             distance = T)
```
If we set that to FALSE, we observed the presence and absence of the domains across the genome samples.

```r
plot_heatmap(important_PFAMs, 
             y_axis=Pfam, 
             analysis= "INTERPRO", 
             distance = F)
```
