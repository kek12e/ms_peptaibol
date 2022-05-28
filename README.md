# ms_peptaibol
Code for biology analyses for publication: Trachymyrmex septentrionalis ants promote fungus garden hygiene using Trichoderma-derived metabolite cues
- ForPub_ms-peptaibol_TSinf.R
  -  R code for processing the phyloseq object (created with unpairedITS_DADA2.r) for the *Time-course Infection ITS* dataset
- ForPub_ms-peptaibol_maits.R
  - R code for processing the phyloseq object (created with unpairedITS_DADA2.r) for the *Environmental ITS* dataset
- alldata.tidy.csv
  - all waste quantification data 
  - used by weeding_quant.R
- master_Su2020antExperiments-alldata.tidy.crude.ncdups.csv
  - crude/dmso/NT only waste quantification data 
  - used by weeding_quant.R
- unpairedITS_DADA2.r
  - R code for running modified DADA2 ITS pipeline on fungal community amplicon sequencing FASTQ files
  - used by ForPub_ms-peptaibol_TSinf.R and ForPub_ms-peptaibol_maits.R
- weeding_quant.R
  - R code for wrangling and plotting the *Extract Bioassays* waste quantification data
