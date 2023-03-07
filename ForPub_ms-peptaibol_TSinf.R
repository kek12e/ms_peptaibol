# 
#   First run the my slightly modified (unpaired Trimmomatic, ITSx) DADA2 ITS pipeline from bash: 
#   `"./unpairedITS_DADA2.R > unpairedITS_DADA2.R.out 2>&1"`
# 
# Adjust ITS2 fasta for manipulation: 
#   `cat ITSx_psASV.ITS2.fasta | paste - - | sed -E 's/>ASV(.*)\|.\|/ASV\1\t/' > ITSx_psASV.ITS2.tsv`
# 
# Check if any ASVs did not pass ITSx and filter those out of phyloseq object:
#   `cat ITSx_psASV_no_detections.fasta`


#install.packages("devtools")
library("devtools")
devtools::install_github("benjjneb/dada2", ref="v1.16") # change the ref argument to get other versions
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# BiocManager::install("phyloseq")
library(dada2)
library(phyloseq)
packageVersion("phyloseq")
library(ShortRead)
library(Biostrings)
library(dplyr)
library(tidyr)
library("data.table")
library(ggplot2)
# install.packages("remotes")
# remotes::install_github("vmikk/metagMisc")
library(metagMisc)
library(ggtext)
library(stringr)
#install.packages("decontam")
# BiocManager::install("decontam")
library(decontam); packageVersion("decontam")
library(tibble)
library(hrbrthemes)

set.seed(14)

# paths setup
ps.path="ps.RDS"
itsx.path="ITSx_psASV.ITS2.tsv"
md.path="../TSinf.itsx.ps.sampledata.csv"
unite.ref="../../sh_general_release_dynamic_s_04.02.2020.fasta"

# read in metdata file... should come up with standard format for this. 
#                         Mainly just the rownames/samplenames that matter for phyloseq
metadat=read.csv(md.path, header=TRUE, colClasses="character",strip.white=TRUE, row.names=1)

# colorblind friendly pallete
cborange="#E69F00";cbblue="#0072B2";cbdorange="#D55E00";cbgreen="#009E73";cbyellow="#F0E442";cblblue="#56B4E9";cbpink="#CC79A7";cbblack="#000000";cbgray="#999999"

#############################
### Time-series Infection ITS
#############################

# load phyloseq object generated in unpairedITS_dada2.R
ps = readRDS(ps.path)
# add sample_data
if( length( which( rownames(metadat) != sample_names(ps) ) ) == 0 ) {
  sample_data(ps) = metadat
} else { 
  print("ERROR: metadata rownames not equivalent to otu_tab rownames (samplenames)")
}

# get ITSx treated asvs
asvs.its2=read.csv( itsx.path, 
                    sep = "\t", 
                    header=FALSE, 
                    colClasses="character", 
                    strip.white=TRUE
                  )

# Check if any ASVs did not pass ITSx and filter those out of phyloseq object:
setdiff(taxa_names(ps),asvs.its2$V1)
ps=prune_taxa(asvs.its2$V1,ps)

### merge ITSx caused ASV dups and remove duplicate ASVs from ITSx refseq table
# duplicated() determines which elements of a vector or data frame are duplicates of elements with smaller subscripts
w.dup=which(duplicated(asvs.its2$V3))           # these indices are duplicates of smaller indices
dup.m=match(asvs.its2[w.dup,]$V3,asvs.its2$V3)  # these are the indices that the dups match to
asvs.its2[w.dup,]$V3 == asvs.its2[dup.m,]$V3    # check

if( length(w.dup) > 0 ) {
  for(i in seq_along(w.dup)){
    # print(taxa_names(ps)[c(dup.m[i],w.dup[i])])          # match results must come first so duplicates are merged 
    # into the smaller indexed ASV
    merge.ps=merge_taxa(merge.ps,taxa_names(ps)[c(w.dup[i],dup.m[i])])
  }
  seqs.itsx = Biostrings::DNAStringSet(asvs.its2[-w.dup,]$V3)
} else {
  seqs.itsx = Biostrings::DNAStringSet(asvs.its2$V3)
}
otu.merge=otu_table(merge.ps)
names(seqs.itsx) = taxa_names(otu.merge)

# recall taxonomy for ITSx treated ASVs
x1 = otu.merge
colnames(x1) = as.character(seqs.itsx)
x1.tax=assignTaxonomy(x1,unite.ref,multithread = T, tryRC = T)
rownames(x1.tax) = colnames(otu.merge)

# remake phyloseq object with new ITSx seqs and merged ASV3/5
ps.itsx <- phyloseq(tax_table(x1.tax), 
                    otu_table(otu.merge, taxa_are_rows = F), 
                    refseq(seqs.itsx)
                   )
taxa_names(ps.itsx) = taxa_names(otu.merge)

# save new phyloseq obj
saveRDS(ps.itsx,"ps.itsx.RDS")


# set genus name to "Cultivar Fungus" if unclassified f__Agaricaceae
w.Fagar=which(tax_table(ps.itsx)[,"Family"] == "f__Agaricaceae")
w.gna=which(is.na(tax_table(ps.itsx)[,"Genus"]))
w.int=intersect(w.Fagar,w.gna)

tax_table(ps.itsx)[w.int,"Genus"]="Cultivar Fungus"

  
  
## DECONTAM
# metadat=read.csv("../TSinf.itsx.ps.sampledata.csv",header = T,row.names = 1)  # moved to setup block
sample_data(ps.itsx)=metadat
sample_data(ps.itsx)

w.quant0=which(sample_data(ps.itsx)$quant_reading == 0)
sample_data(ps.itsx)$quant_reading[w.quant0] = 0.00001
contamdf.freq=isContaminant(ps.itsx, method="frequency", conc = "quant_reading")
table(contamdf.freq$contaminant)  # no contaminants this way

contamdf.prev <- isContaminant(ps.itsx, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)  # also none

contamdf.both.batch=isContaminant(ps.itsx, method="combined", batch="batch", conc="quant_reading", neg="is.neg")
table(contamdf.both.batch$contaminant)
# contamdf.both.batch    # also also none

contamdf.both=isContaminant(ps.itsx, method="combined", conc="quant_reading", neg="is.neg")
table(contamdf.both$contaminant)
# contamdf.both    # also also also none

contamdf.min.batch=isContaminant(ps.itsx, method="minimum", batch="batch", conc="quant_reading", neg="is.neg")
table(contamdf.min.batch$contaminant)
# contamdf.min.batch    # also also also also none

contamdf.eith.batch=isContaminant(ps.itsx, method="either", batch="batch", conc="quant_reading", neg="is.neg")
table(contamdf.eith.batch$contaminant)
# contamdf.min.batch    # also also also also also none

contamdf.eith=isContaminant(ps.itsx, method="either", conc="quant_reading", neg="is.neg")
table(contamdf.eith$contaminant)
# contamdf.min.batch    # also also also also also also none

### No contaminants so we can remove the negative control samples. 
tf.blank=sample_data(ps.itsx)$sample_type == "blank"
sample_data(ps.itsx)[tf.blank,]
ps.itsx.no_nc=prune_samples(!tf.blank, ps.itsx)


# ##----- rarefy -----
track=read.csv("track.csv",row.names = 1,header=T)
sort(sample_sums(ps.itsx.no_nc))
barplot(sort(sample_sums(ps.itsx.no_nc))); abline(h=4099, col='red')
#       # rarefy to 4099 reads
set.seed(14); ps.itsx.no_nc.raref = rarefy_even_depth(ps.itsx.no_nc, sample.size = 4099)
#                   2 samples removedbecause they contained fewer reads than `sample.size`.
#                   Up to first five removed samples are: 
#                       TSe26TSe04
sort(sample_sums(ps.itsx.no_nc.raref))


# ## filter highly unclass ASVs
which(tax_table(ps.itsx.no_nc.raref)[,"Kingdom"] != "k__Fungi")   # all classified as fungi
w.pna=which(is.na(tax_table(ps.itsx.no_nc.raref))[,"Phylum"])     # two ASVs unclassified except for k__Fungi (29, 35)
tax_table(ps.itsx.no_nc.raref)[w.pna,]
# ps.itsx.no_nc.raref.filt=subset_taxa(ps.itsx.no_nc.raref, !is.na(Phylum)) # now 31 taxa
# set to Other instead of removing
tax_table(ps.itsx.no_nc.raref)[w.pna,1:7] = "Other"

# ## threshold filtering i.e. ASVs <1% abundance in at least one sample}
library(metagMisc)
ps.itsx.no_nc.raref.filt.gt1perc=phyloseq_filter_sample_wise_abund_trim(ps.itsx.raref.filt, minabund = 41, relabund = FALSE, rm_zero_OTUs = TRUE)  # 40.99 reads is 1% of 4099 (which is what smaples are raref to)
tax_table(ps.itsx.no_nc.raref.filt.gt1perc) # 14 ASVs now
# list of ASVs that were filtered by sample_wise...
gt1perc.filt.ASVs=setdiff(rownames(tax_table(ps.itsx.no_nc.raref)),rownames(tax_table(ps.itsx.no_nc.raref.filt.gt1perc)))
# set all those ASVs genus to "Other" in nonfilt so will glom together
Other.gt1perc.ps.itsx.no_nc.raref=ps.itsx.no_nc.raref
tax_table(Other.gt1perc.ps.itsx.no_nc.raref)[gt1perc.filt.ASVs,1:7] = "Other"

# 1. adjust unclassified genera
# 2. combine cultivar
# 3. glom by genus
# 4. filter based on gt1perc genera (bc glommed, prev filter was by ASV?)
# - dont think this is really necessary... should be fine to just plot since already filtered on 1% thresh
# 5. calc and check prevalence
# 6. plot

# ## 1. adjust unclassified genera
# # define unclassified genera by family so no taxa lost
gen.na=is.na(tax_table(Other.gt1perc.ps.itsx.no_nc.raref)[,"Genus"])
if(length(which(gen.na)) > 0) {
  unclass.string=paste("unclassified",tax_table(Other.gt1perc.ps.itsx.no_nc.raref)[gen.na,"Family"])
  tax_table(Other.gt1perc.ps.itsx.no_nc.raref)[gen.na,"Genus"] = unclass.string
}

# ## 2. combine cultivar
# # combine leucocoprinus and unclassified agaricaecae
unique(tax_table(Other.gt1perc.ps.itsx.no_nc.raref)[,"Genus"])
w.cult=which(tax_table(Other.gt1perc.ps.itsx.no_nc.raref)[,"Genus"] == "g__Leucocoprinus" |
               tax_table(Other.gt1perc.ps.itsx.no_nc.raref)[,"Genus"] == "g__Leucoagaricus" |
               tax_table(Other.gt1perc.ps.itsx.no_nc.raref)[,"Genus"] == "unclassified f__Agaricaceae") 
if(length(w.cult) > 0) {
  tax_table(Other.gt1perc.ps.itsx.no_nc.raref)[w.cult,"Genus"] = "Cultivar Fungus"
  tax_table(Other.gt1perc.ps.itsx.no_nc.raref)[w.cult,] 
}

# ## 3. glom by genus
Other.gt1perc.ps.itsx.no_nc.raref.gglom=tax_glom(Other.gt1perc.ps.itsx.no_nc.raref, "Genus" , NArm=TRUE)
Other.gt1perc.ps.itsx.no_nc.raref.gglom     # now 5 taxa
tax_table(Other.gt1perc.ps.itsx.no_nc.raref.gglom)[,"Genus"]


# # # prevalence 
# prev.Other.gt1perc.ps.itsx.no_nc.raref.gglom = apply(X = otu_table(Other.gt1perc.ps.itsx.no_nc.raref.gglom),
#                 # MARGIN = 2,   # this needs to be 2 if taxa are columns instead of rows!
#                  MARGIN = ifelse(taxa_are_rows(Other.gt1perc.ps.itsx.no_nc.raref.gglom), 
#                                  yes = 1, no = 2),# FANCY https://www.bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/doc/MicrobiomeWorkshopII.html#prevalence-filtering
#                  FUN = function(x){sum(x > 0)})
# # # Add taxonomy and total read counts to this data.frame
# prev.Other.gt1perc.ps.itsx.no_nc.raref.gglom = data.frame(Prevalence = prev.Other.gt1perc.ps.itsx.no_nc.raref.gglom,
#                       TotalAbundance = taxa_sums(Other.gt1perc.ps.itsx.no_nc.raref.gglom),
#                       # RelGlomAbund.perc = 100*taxa_sums(ps.itsx.raref.filt.gt1perc.gglom)/sum( taxa_sums( ps.itsx.raref.filt.gt1perc.gglom )),
#                       # RelNonCultAbund.perc=100*taxa_sums(filtNot0out.glom.genus)/sum(taxa_sums(nonCultivar.raref.itsx.ps)),
#                       RelRarefAbund.perc=100*taxa_sums(Other.gt1perc.ps.itsx.no_nc.raref.gglom)/sum(taxa_sums(Other.gt1perc.ps.itsx.no_nc.raref.gglom)),
#                       tax_table(Other.gt1perc.ps.itsx.no_nc.raref.gglom))
# prev.Other.gt1perc.ps.itsx.no_nc.raref.gglom[order(-prev.Other.gt1perc.ps.itsx.no_nc.raref.gglom$Prevalence),]
# prev.Other.gt1perc.ps.itsx.no_nc.raref.gglom[order(-prev.Other.gt1perc.ps.itsx.no_nc.raref.gglom$RelRarefAbund.perc),]


## plot stacked bar
# calc REL abundance
Other.gt1perc.ps.itsx.no_nc.raref.gglom.prop <- transform_sample_counts(Other.gt1perc.ps.itsx.no_nc.raref.gglom, function(otu) otu/sum(otu))
# stacked bar plot
plot_bar(Other.gt1perc.ps.itsx.no_nc.raref.gglom.prop, fill = "Genus") +
  aes(x=sample_descrip) +
  facet_wrap(~treatment,scales = "free_x") +
  scale_fill_manual(values=c(cborange,cbpink,cbblue,cbgreen,cblblue,cbyellow,"white"))


##----- investigate trichoderma only----

trichonly.Other.gt1perc.ps.itsx.no_nc.raref.gglom=subset_taxa(Other.gt1perc.ps.itsx.no_nc.raref.gglom, Genus=="g__Trichoderma")
trichonly.Other.gt1perc.ps.itsx.no_nc.raref.gglom.prop=subset_taxa(Other.gt1perc.ps.itsx.no_nc.raref.gglom.prop, Genus=="g__Trichoderma")

trichonly.Other.gt1perc.ps.itsx.no_nc.raref.gglom
sort(sample_sums(trichonly.Other.gt1perc.ps.itsx.no_nc.raref.gglom))
sum(sample_sums(trichonly.Other.gt1perc.ps.itsx.no_nc.raref.gglom))   # 22958


plot_bar(trichonly.Other.gt1perc.ps.itsx.no_nc.raref.gglom, fill = "Genus") +
  theme(axis.text=element_text(size=10),legend.text = element_text(size=9)) +
  facet_wrap(~treatment, scales="free_x") +
  ggtitle("Time Series Trich. Infection - Trich only ITSx, filtered, rarefy 4009 -- Absolute Abund.") +
  scale_fill_manual(values=cbgreen)

plot_bar(trichonly.Other.gt1perc.ps.itsx.no_nc.raref.gglom.prop, fill = "Genus") +
  theme(axis.text=element_text(size=10),legend.text = element_text(size=9)) +
  facet_wrap(~treatment, scales="free_x") +
  ggtitle("Time Series Trich. Infection - Trich only ITSx, filtered, rarefy 4009 -- Rel. Abund.") +
  scale_fill_manual(values=cbgreen)

## plot trich heatmap

# sample sums is in same order as sample_data
which(names(sample_sums(trichonly.Other.gt1perc.ps.itsx.no_nc.raref.gglom.prop)) != rownames(sample_data(trichonly.Other.gt1perc.ps.itsx.no_nc.raref.gglom.prop)))
df.trichonly.Other.gt1perc.ps.itsx.no_nc.raref.gglom.prop=sample_data(trichonly.Other.gt1perc.ps.itsx.no_nc.raref.gglom.prop)
df.trichonly.Other.gt1perc.ps.itsx.no_nc.raref.gglom.prop$sample.sums=sample_sums(trichonly.Other.gt1perc.ps.itsx.no_nc.raref.gglom.prop)
df.trichonly.Other.gt1perc.ps.itsx.no_nc.raref.gglom.prop


df.trichonly.Other.gt1perc.ps.itsx.no_nc.raref.gglom.prop %>% 
  ggplot(aes(day,sample_descrip)) + geom_tile(aes(fill = sample.sums))

write.csv(df.trichonly.itsx.ps.filt.raref.prop, "df.trichonly.itsx.ps.filt.raref.prop.csv")

## finetune trich heatmap - Fig 2b
hm=read.csv("heatmap.df.trichonly.itsx.ps.filt.raref.prop.csv",header = T)

hm %>% 
  ggplot(aes(heatmap.x.day,heatmap.y.sampTypeXtreatXantsXrep)) + 
  geom_tile(aes(fill = heatmap.z.trichPerc), color='black') + 
  scale_fill_gradient(low = "lightyellow", high = cbgreen, na.value = "black") +
  geom_text(aes(label=sprintf("%0.2f", round(heatmap.z.trichPerc, digits = 2)))) +
  theme(
    axis.title.x = element_text(size = 12, color='black'),
    panel.background = element_rect(fill = 'grey90',color='black'),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color="darkgrey"),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size=12,color = 'black'),
    legend.text = element_text(size=12,color='black'),
  )

