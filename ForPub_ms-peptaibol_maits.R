# First run the my slightly modified (unpaired Trimmomatic, ITSx) DADA2 ITS pipeline from bash: 
#   `"./unpairedITS_DADA2.R > unpairedITS_DADA2.R.out 2>&1"`
# 
# Adjust ITS2 fasta for manipulation: 
#   `cat ITSx_psASV.ITS2.fasta | paste - - | sed -E 's/>ASV(.*)\|.\|/ASV\1\t/' > ITSx_psASV.ITS2.tsv`
# 
# Check if any ASVs did not pass ITSx and filter those out of phyloseq object:
#   `cat ITSx_psASV_no_detections.fasta`


md.path="../maits.metadata.csv"
unite.ref="../../sh_general_release_dynamic_s_04.02.2020.fasta"
ps.path="ps.RDS"
itsx.path="ITSx_psASV.ITS2.tsv"




library("devtools")
#devtools::install_github("benjjneb/dada2", ref="v1.16") # change the ref argument to get other versions
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
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

set.seed(14)

cborange="#E69F00";cbblue="#0072B2";cbdorange="#D55E00";cbgreen="#009E73";cbyellow="#F0E442";cblblue="#56B4E9";cbpink="#CC79A7";cbblack="#000000";cbgray="#999999";
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", cbblack, cbgray) # can also add black or gray


## MAITS

# load phyloseq object generated in unpairedITS_dada2.R
ps = readRDS(ps.path)
asvs.its2=read.csv(itsx.path,sep = "\t",header=FALSE,colClasses="character",strip.white=TRUE)

metadat=read.csv(md.path, header=TRUE, colClasses="character",strip.white=TRUE)
rownames(metadat) = metadat$JKS_treat
sample_data(ps) = metadat

# filter out MAITS weird sample issues
g=grepl("^3656$|^3501$|^3060$|^1822$|^1823$|^1828$|^1409$|^642$|^288$",rownames(sample_data(ps)))
ps=prune_samples(!g,ps)


# Check if any ASVs did not pass ITSx and filter those out of phyloseq object:
setdiff(taxa_names(ps),asvs.its2$V1)
ps=prune_taxa(asvs.its2$V1,ps)


# ITSx processing
# duplicated() determines which elements of a vector or data frame are duplicates of elements with smaller subscripts
w.dup=which(duplicated(asvs.its2$V3))           # these indices are duplicates of smaller indices
dup.m=match(asvs.its2[w.dup,]$V3,asvs.its2$V3)  # these are the indices that the dups match to
# cbind(w.dup,dup.m)
asvs.its2[w.dup,]$V3 == asvs.its2[dup.m,]$V3

# merged taxa/otu tables
merge.ps=ps
for(i in seq(w.dup)){
  # print(taxa_names(ps)[c(dup.m[i],w.dup[i])])          # match results must come first so duplicates are merged 
  # into the smaller indexed ASV
  merge.ps=merge_taxa(merge.ps,taxa_names(ps)[c(w.dup[i],dup.m[i])])
}
otu.merge=otu_table(merge.ps)

# remove duplicate ASVs from ITSx refseq table
seqs.itsx = Biostrings::DNAStringSet(asvs.its2[-w.dup,]$V3)
names(seqs.itsx) = taxa_names(otu.merge)

# recall taxonomy for ITSx treated ASVs
x1 = otu.merge
colnames(x1) = as.character(seqs.itsx)
x1.tax=assignTaxonomy(x1,unite.ref,multithread = T, tryRC = T)
rownames(x1.tax) = colnames(otu.merge)

# remake phyloseq object with new ITSx seqs and merged ASV3/5
ps.itsx <- phyloseq(tax_table(x1.tax), otu_table(otu.merge, taxa_are_rows = F))
ps.itsx <- merge_phyloseq(ps.itsx,seqs.itsx,sample_data(merge.ps))
taxa_names(ps.itsx) = taxa_names(otu.merge)


# ##----- rarefy -----
ps.itsx.no_nc = ps.itsx
track=read.csv("track.csv",row.names = 1,header=T)
sort(sample_sums(ps.itsx.no_nc))
# biggest read disparity in this range... 
#    ID: 1243   3013    827   3450   3113   3196    614
# reads: 249    978   1456   5701  12531  12556  15995
## rarefy to 12531 reads? else 15995?
barplot(sort(sample_sums(ps.itsx.no_nc))); abline(h=15000, col='red'); 
			abline(h=12531, col='blue'); abline(h=15995, col='purple')


set.seed(14); ps.itsx.no_nc.raref = rarefy_even_depth(ps.itsx.no_nc, sample.size = 12531)
#                   7 samples removedbecause they contained fewer reads than `sample.size`.
#                   Samples removed:
#                     3440   3156   3634   1243   3013    827   3450               
#                   211 OTUs were removed because they are no longer 
#                   present in any sample after random subsampling


## filter highly unclass ASVs
w.kna=which(is.na(tax_table(ps.itsx.no_nc.raref)[,"Kingdom"]))
w.pna=which(is.na(tax_table(ps.itsx.no_nc.raref))[,"Phylum"])     
length(w.kna)          # one not classified as fungi
length(w.pna)          # 16 ASVs unclassified except for k__Fungi (this includes the unclassified kingdom taxa)
tax_table(ps.itsx.no_nc.raref)[c(w.pna),]
# if NA at kingdom will also be NA at phylum
ps.itsx.no_nc.raref.filt=subset_taxa(ps.itsx.no_nc.raref, !is.na(Phylum)) # now 436 taxa 


## threshold filtering i.e. ASVs <1% abundance in at least one sample}
library(metagMisc)
ps.itsx.no_nc.raref.filt.gt1perc=phyloseq_filter_sample_wise_abund_trim(ps.itsx.no_nc.raref.filt, 
																		minabund = 125, relabund = FALSE, 
																		rm_zero_OTUs = TRUE)  # 125.31 reads is 1% of 12,531 (which is what smaples are raref to)
tax_table(ps.itsx.no_nc.raref.filt.gt1perc) # 35 ASVs now


# 1. adjust unclassified genera
# 2. combine cultivar
# 3. glom by genus
# 4. filter based on gt1perc genera (bc glommed, prev filter was by ASV?)
# - dont think this is really necessary... should be fine to just plot since already filtered on 1% thresh
# 5. calc and check prevalence
# 6. plot


# ## 1. adjust unclassified genera
# # define unclassified genera by family so no taxa lost
gen.na=is.na(tax_table(ps.itsx.no_nc.raref.filt.gt1perc)[,"Genus"])
if(length(which(gen.na)) > 0) {
  unclass.string=paste("unclassified",tax_table(ps.itsx.no_nc.raref.filt.gt1perc)[gen.na,"Family"])
  tax_table(ps.itsx.no_nc.raref.filt.gt1perc)[gen.na,"Genus"] = unclass.string
}

# ## 2. combine cultivar
# # combine leucocoprinus and unclassified agaricaecae
unique(tax_table(ps.itsx.no_nc.raref.filt.gt1perc)[,"Genus"])
w.cult=which(tax_table(ps.itsx.no_nc.raref.filt.gt1perc)[,"Genus"] == "g__Leucocoprinus" |
               tax_table(ps.itsx.no_nc.raref.filt.gt1perc)[,"Genus"] == "g__Leucoagaricus" |
               tax_table(ps.itsx.no_nc.raref.filt.gt1perc)[,"Genus"] == "unclassified f__Agaricaceae") 
if(length(w.cult) > 0) {
  tax_table(ps.itsx.no_nc.raref.filt.gt1perc)[w.cult,"Genus"] = "Cultivar Fungus"
#  tax_table(ps.itsx.no_nc.raref.filt.gt1perc)[w.cult,] 	# prints if you wanna see
}

# ## 3. glom by genus
ps.itsx.no_nc.raref.filt.gt1perc.gglom=tax_glom(ps.itsx.no_nc.raref.filt.gt1perc, "Genus" , NArm=TRUE)
ps.itsx.no_nc.raref.filt.gt1perc.gglom     # now 10 taxa
tax_table(ps.itsx.no_nc.raref.filt.gt1perc.gglom)[,"Genus"]


# # prevalence 
prev.itsx.no_nc.raref.filt.gt1perc.gglom = apply(X = otu_table(ps.itsx.no_nc.raref.filt.gt1perc.gglom),
                                                 # MARGIN = 2,   # this needs to be 2 if taxa are columns instead of rows!
                                                 MARGIN = ifelse(taxa_are_rows(ps.itsx.no_nc.raref.filt.gt1perc.gglom), 
                                                                 yes = 1, no = 2),# FANCY https://www.bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/doc/MicrobiomeWorkshopII.html#prevalence-filtering
                                                 FUN = function(x){sum(x > 0)})
# # Add taxonomy and total read counts to this data.frame
prev.itsx.no_nc.raref.filt.gt1perc.gglom = 
		data.frame(Prevalence = prev.itsx.no_nc.raref.filt.gt1perc.gglom,
        			TotalAbundance = taxa_sums(ps.itsx.no_nc.raref.filt.gt1perc.gglom),
            		# RelGlomAbund.perc = 100*taxa_sums(ps.itsx.raref.filt.gt1perc.gglom)/sum( taxa_sums( ps.itsx.raref.filt.gt1perc.gglom )),
  					# RelNonCultAbund.perc=100*taxa_sums(filtNot0out.glom.genus)/sum(taxa_sums(nonCultivar.raref.itsx.ps)),
 					RelRarefAbund.perc=100*(taxa_sums(ps.itsx.no_nc.raref.filt.gt1perc.gglom)/
                                        sum(taxa_sums(ps.itsx.no_nc.raref.filt.gt1perc.gglom))),
           	 		tax_table(ps.itsx.no_nc.raref.filt.gt1perc.gglom)[,"Genus"]
        )
prev.itsx.no_nc.raref.filt.gt1perc.gglom[order(-prev.itsx.no_nc.raref.filt.gt1perc.gglom$Prevalence),]
prev.itsx.no_nc.raref.filt.gt1perc.gglom[order(-prev.itsx.no_nc.raref.filt.gt1perc.gglom$RelRarefAbund.perc),]


## add "Other" taxa for ASVs <1% abundant
# list of ASVs that were filtered by sample_wise...
gt1perc.filt.ASVs=setdiff(rownames(tax_table(ps.itsx.no_nc.raref.filt)),rownames(tax_table(ps.itsx.no_nc.raref.filt.gt1perc)))
# set all those ASVs genus to "Other" in nonfilt so will glom together
ps.itsx.no_nc.raref.filt.gt1percOther=ps.itsx.no_nc.raref.filt
tax_table(ps.itsx.no_nc.raref.filt.gt1percOther)[gt1perc.filt.ASVs,1:7] = "Other"
# ## 1. adjust unclassified genera
# # define unclassified genera by family so no taxa lost
gen.na=is.na(tax_table(ps.itsx.no_nc.raref.filt.gt1percOther)[,"Genus"])
if(length(which(gen.na)) > 0) {
  unclass.string=paste("unclassified",tax_table(ps.itsx.no_nc.raref.filt.gt1percOther)[gen.na,"Family"])
  tax_table(ps.itsx.no_nc.raref.filt.gt1percOther)[gen.na,"Genus"] = unclass.string
}
# ## 2. combine cultivar
# # combine leucocoprinus and unclassified agaricaecae
unique(tax_table(ps.itsx.no_nc.raref.filt.gt1percOther)[,"Genus"])
w.cult=which(tax_table(ps.itsx.no_nc.raref.filt.gt1percOther)[,"Genus"] == "g__Leucocoprinus" |
               tax_table(ps.itsx.no_nc.raref.filt.gt1percOther)[,"Genus"] == "g__Leucoagaricus" |
               tax_table(ps.itsx.no_nc.raref.filt.gt1percOther)[,"Genus"] == "unclassified f__Agaricaceae") 
if(length(w.cult) > 0) {
  tax_table(ps.itsx.no_nc.raref.filt.gt1percOther)[w.cult,"Genus"] = "Cultivar Fungus"
#  tax_table(ps.itsx.no_nc.raref.filt.gt1percOther)[w.cult,] 	# prints
}
# ## 3. glom by genus
ps.itsx.no_nc.raref.filt.gt1percOther.gglom=tax_glom(ps.itsx.no_nc.raref.filt.gt1percOther, "Genus" , NArm=TRUE)
ps.itsx.no_nc.raref.filt.gt1percOther.gglom     # now 11 taxa
tax_table(ps.itsx.no_nc.raref.filt.gt1percOther.gglom)[,"Genus"]


# relative abundance
ps.itsx.no_nc.raref.filt.gt1percOther.gglom.prop <- 
			transform_sample_counts(ps.itsx.no_nc.raref.filt.gt1percOther.gglom, 
			function(otu) otu/sum(otu))
# top20 = names(sort(taxa_sums(ps.itsx), decreasing = T))[1:20]
# ps.itsx.prop.top20=prune_taxa(top20,ps.itsx.prop)

### Stacked bar chart for Figure S2:
plot_bar(ps.itsx.no_nc.raref.filt.gt1percOther.gglom.prop, fill="Genus", 
         title = "Filtered/Rarefied Wild Fungus Garden ITS w/ Other cat.") + 
  theme(axis.text=element_text(size=10)) +
  scale_fill_manual(values=c(cborange,cbpink,cblblue,cbblue,cbgray,cbyellow,cbblack,cbgreen,"grey35", "white",cbdorange)) +
  facet_grid(~State,scales="free_x", space="free_x", switch = "x") +
  labs(x=NULL,y="Relative Abundance") +
  theme(#aspect.ratio = 0.5,
    panel.background = element_blank(),    
    axis.text.x = element_blank(),
    legend.text = element_markdown(),
    #        legend.position = "none",
    #        legend.background = element_rect(color="black", fill = NA, size=0.25),
    #        legend.margin = margin(t=-5, r=3, b=3)
    panel.grid  = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_blank(),
    strip.switch.pad.grid = unit('00.1',"cm")
  )

### Wrangle data for Fig 1A:

## DATAFRAMES #3 (Other cat)
metadat3=as.data.frame(as.matrix(sample_data(ps.itsx.no_nc.raref.filt.gt1percOther.gglom.prop)))
colnames(metadat3)[1]="sampleID"

otucounts3=as.data.frame(otu_table(ps.itsx.no_nc.raref.filt.gt1percOther.gglom))
otucounts3$sampleID=rownames(otucounts3)
otucounts3 <- otucounts3 %>% dplyr::select(sampleID,starts_with("ASV"))  %>% 
					pivot_longer(-sampleID, values_to = "count", names_to = "ASV")

tax_table(ps.itsx.no_nc.raref.filt.gt1percOther.gglom)[,"Species"]="NA"
taxonomy3 <- as.data.frame(tax_table(ps.itsx.no_nc.raref.filt.gt1percOther.gglom))
taxonomy3$ASV=rownames(taxonomy3)
taxonomy3 <- dplyr::select(taxonomy3, ASV, 1:7)

oturelabund3 <- inner_join(metadat3, otucounts3, by="sampleID") %>% inner_join(., taxonomy3, by="ASV") %>% 
  group_by(sampleID) %>%
  mutate(wiSample_rel_abund = count / sum(count)) %>%
  ungroup() %>%
  dplyr::select(-count) %>%
  pivot_longer(
    c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "ASV"),
    names_to="level",
    values_to="taxon") #%>%

taxonrelabund3 <- 
  oturelabund3 %>%
  filter(level=="Genus") %>%
  group_by(sampleID, taxon) %>%
  summarize(wiSample_rel_abund = 100*sum(wiSample_rel_abund), .groups="drop") %>%  
  mutate(taxon = str_replace_all(taxon,".__", "*"),
         taxon = str_replace(taxon, "unclassified (.*)", "Unclassified<br>\\1*"),
         taxon = str_replace(taxon,"(^\\*.*)$", "\\1*"),
  ) %>% 
  inner_join(.,metadat3, by="sampleID")

# all ~12531 (what we raref to)
nseqs_per_sample3 <- otucounts3 %>%
  group_by(sampleID) %>%
  summarize(N = sum(count), .groups="drop") 

lod.approx3 <- 100* 1/12531

## filter out weird maits samples
# this should already be out though
taxonrelabund3.mafilt <-taxonrelabund3 %>%
  group_by(sampleID) %>% 
  filter(!grepl("^3656_|^3501_|^3060_|^1822_|^1823_|^1828_|^1409_|^642_",sampleID))


## save this final values!
# final values with Other taxon
finalTableForPlotWOther <- taxonrelabund3.mafilt %>%
  group_by(sampleID, taxon) %>% 
  mutate(wiSample_rel_abund = if_else(wiSample_rel_abund == 0,
                             2/3 * lod.approx3,
                             wiSample_rel_abund)) %>% 
  group_by(taxon) %>% 
  mutate(fullDataset_mean_rel_abund = mean(wiSample_rel_abund)) %>% 
  mutate(fullDataset_med_rel_abund = median(wiSample_rel_abund))
write.csv(finalTableForPlotWOther,"finalMAITSvaluesWOther_Fig1A.csv")
  
# final values with Other removed (for plot)
finalTableForPlot <- taxonrelabund3.mafilt %>%
  group_by(sampleID, taxon) %>% 
  filter(!grepl("Other",taxon)) %>%           # remove genus 'Other'
  mutate(wiSample_rel_abund = if_else(wiSample_rel_abund == 0,
                             2/3 * lod.approx3,
                             wiSample_rel_abund)) %>% 
  group_by(taxon) %>% 
  mutate(fullDataset_mean_rel_abund = mean(wiSample_rel_abund)) %>% 
  mutate(fullDataset_med_rel_abund = median(wiSample_rel_abund))
write.csv(finalTableForPlotWOther,"finalMAITSvalues_Fig1A.csv")


########## Fig 1 A FINAL
## jitter plot
set.seed(14); finalTableForPlot %>%
  ggplot(aes(y=wiSample_rel_abund, 
             #            x=taxon,
             x=reorder(taxon,-fullDataset_mean_rel_abund),
             #    x=reorder(taxon,-mean_rel_abund),
             #            x=taxon,
             #             fill=State
             color=taxon
  ) ) +
  geom_boxplot(outlier.colour = NA, 
               #               size =1
  ) +
  geom_jitter( 
    #    aes(fill=State),
    shape=20,
    stroke=NA,
    width = 0.2,
    #   stroke=0,
    size=0.001,
    alpha=0.5 
  ) +
  geom_hline(yintercept = lod.approx3, size=0.2, lty=2, alpha =0.2) +
  # geom_hline(yintercept = 0.0334505, color='red') + # clado top whisker
  # geom_hline(yintercept = 0.0133802, color='blue') +  # clado top notch
  
  #  geom_violin(scale = "count") +
  # stat_summary(fun=mean, geom = "crossbar",
  #              width=0.7,
  #              size=0.5,
  #              alpha=0.5,
  #              color="gray40",
  #              show.legend=FALSE) +
  # stat_summary(fun=median, geom = "crossbar",
#            width=0.4,
#            size=0.25,
#            #alpha=0.5,
#            color="blue",
#            show.legend=FALSE) +
coord_trans(y="log10") +
  scale_y_continuous(limits=c(NA, 110),
                     breaks=c(0.1, 1, 10, 100),
                     labels=c(0.1, 1, 10, 100)
  ) +
  #   scale_fill_manual(values=c(cbPallete),
  #                      name = NULL,
  # #                     breaks = c("NY", "NJ", "NC", "GA", "FL", "LA"), #fl ga la nc nj ny 
  #                                                                       # 9 12  9 11 48  1
  #                      labels=c("FL<br>(n = 9)", "GA<br>(n = 12)", "LA<br>(n = 9)", 
  #                               "NC<br>(n = 11)", "NJ<br>(n = 48)", "NY<br>(n = 1)")) +
  scale_color_manual(values=c(cbpink,cblblue, cbblue,cbgray, cbyellow,cbblack,cbgreen,"grey35", cborange, cbdorange),
                     name = NULL,
                     #                     breaks = c("NY", "NJ", "NC", "GA", "FL", "LA"), #fl ga la nc nj ny 
                     # 9 12  9 11 48  1
                     # labels=c("FL<br>(n = 9)", "GA<br>(n = 12)", "LA<br>(n = 9)", 
                     #          "NC<br>(n = 11)", "NJ<br>(n = 48)", "NY<br>(n = 1)")
  ) +
  labs(x=NULL,
       y="Relative Abundance (%)") +
  theme_classic() +
  guides(color = guide_legend(nrow=1)) +
  theme(#aspect.ratio = 0.5,
    axis.text.x = element_markdown(angle = 45, hjust=1),
    legend.text = element_markdown(),
    legend.position = "none",
    legend.background = element_rect(color="black", fill = NA, size=0.25),
    legend.margin = margin(t=-5, r=3, b=3)
    #        panel.grid.minor.x = element_line(color="gray")
  )


########## Fig 1 A final WITH OTHER cat
## jitter plot
set.seed(14); finalTableForPlotWOther %>%
  ggplot(aes(y=wiSample_rel_abund, 
             #            x=taxon,
             x=reorder(taxon,-fullDataset_mean_rel_abund),
             #    x=reorder(taxon,-mean_rel_abund),
             #            x=taxon,
             #             fill=State
             color=taxon
  ) ) +
  geom_boxplot(outlier.colour = NA, 
               #               size =1
  ) +
  geom_jitter( 
    #    aes(fill=State),
    shape=20,
    stroke=NA,
    width = 0.2,
    #   stroke=0,
    size=0.001,
    alpha=0.5 
  ) +
  geom_hline(yintercept = lod.approx3, size=0.2, lty=2, alpha =0.2) +
  # geom_hline(yintercept = 0.0334505, color='red') + # clado top whisker
  # geom_hline(yintercept = 0.0133802, color='blue') +  # clado top notch
  
  #  geom_violin(scale = "count") +
  # stat_summary(fun=mean, geom = "crossbar",
  #              width=0.7,
  #              size=0.5,
  #              alpha=0.5,
  #              color="gray40",
  #              show.legend=FALSE) +
  # stat_summary(fun=median, geom = "crossbar",
#            width=0.4,
#            size=0.25,
#            #alpha=0.5,
#            color="blue",
#            show.legend=FALSE) +
coord_trans(y="log10") +
  scale_y_continuous(limits=c(NA, 110),
                     breaks=c(0.1, 1, 10, 100),
                     labels=c(0.1, 1, 10, 100)
  ) +
  #   scale_fill_manual(values=c(cbPallete),
  #                      name = NULL,
  # #                     breaks = c("NY", "NJ", "NC", "GA", "FL", "LA"), #fl ga la nc nj ny 
  #                                                                       # 9 12  9 11 48  1
  #                      labels=c("FL<br>(n = 9)", "GA<br>(n = 12)", "LA<br>(n = 9)", 
  #                               "NC<br>(n = 11)", "NJ<br>(n = 48)", "NY<br>(n = 1)")) +
#  scale_color_manual(values=c(cbpink,cblblue, cbblue,cbgray, cbyellow,cbblack,cbgreen,"grey35", cborange, cbdorange,"grey55"),
#                     name = NULL,
#                     #                     breaks = c("NY", "NJ", "NC", "GA", "FL", "LA"), #fl ga la nc nj ny 
#                     # 9 12  9 11 48  1
#                     # labels=c("FL<br>(n = 9)", "GA<br>(n = 12)", "LA<br>(n = 9)", 
#                     #          "NC<br>(n = 11)", "NJ<br>(n = 48)", "NY<br>(n = 1)")
#  ) +
  labs(x=NULL,
       y="Relative Abundance (%)") +
  theme_classic() +
  guides(color = guide_legend(nrow=1)) +
  theme(#aspect.ratio = 0.5,
    axis.text.x = element_markdown(angle = 45, hjust=1),
    legend.text = element_markdown(),
    legend.position = "none",
    legend.background = element_rect(color="black", fill = NA, size=0.25),
    legend.margin = margin(t=-5, r=3, b=3)
    #        panel.grid.minor.x = element_line(color="gray")
  )
