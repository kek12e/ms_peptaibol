#!/usr/bin/env Rscript

## suggested command to run: 
##		./unpairedITS_DADA2.R > unpairedITS_DADA2.R.out 2>&1


# rm(list=ls())		# clear environment variables

######################################################
######### SET THESE VARIABLES TO YOUR SPECIFICS ######
######################################################
#setwd("/Users/kkyle/ms_peptaibol/TSinfITS/test")		# set your working directory, where you want output and .RData

# getting input files
fqpath="./fastq"	# folder containing fastq files
	# ^don't put trailing '/' for this path!
R1.pattern="_R1.fastq.gz"	# unique pattern to match your forward reads

# for primer detection before/after cutadapt
FWD="GCCGGCTGCGACGTGARTCATCGAATCTTTG"   # MARS version (fITS7: GTGARTCATCGAATCTTTG)
REV="AGGCAGTCAGCCTCCTCCGCTTATTGATATGC"  # MARS version (ITS4: TCCTCCGCTTATTGATATGC)

# path to programs wrapped in this script
#cutadapt="/home/kkyle/usr/bin/cutadapt"
cutadapt="/usr/bin/cutadapt"
TrimmomaticSE="/usr/bin/TrimmomaticSE"
ITSx="/usr/local/bin/ITSx"


# path to UNITE fungal db for taxa calling
unite.ref="./sh_general_release_dynamic_s_04.02.2020.fasta"

######################################################

# get latest version of dada2 if don't have it
#library("devtools")
#devtools::install_github("benjjneb/dada2", ref="v1.20") # change the ref argument to get other versions

# load libraries
library(dada2); packageVersion("dada2")
library(phyloseq); packageVersion("phyloseq")
library(ShortRead); packageVersion("ShortRead")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")

#library(dplyr)		# not necessary for DADA2 pipeline, but useful for dataframe manipulation

## ----get input files setup--------------------------------------------------------------------------------------------------------------
fnFs=sort(list.files(fqpath,pattern=R1.pattern,full.names = T))

## ----get sample names-------------------------------------------------------------------------------------------------------------------
# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(fnFs, get.sample.name))
# make sure looks hunky dory
head(cbind(sample.names,fnFs))

## ----Trimmomatic--------------------------------------------------------------------------------------------------------------
trimm.flags = paste("SLIDINGWINDOW:5:20", "MINLEN:125", "-threads 4")
#for i in $(ls /sda/data/KEKdecon3TimeSeriesTestInfITS/Samples/TimeSeriesInfection/fastq/*R1.fastq.gz); do j=`basename -s .fastq.gz $i`; TrimmomaticSE $i ${j}_Trimm.SW5.20.min125.fq.gz SLIDINGWINDOW:5:20 MINLEN:125 -threads 4 -trimlog TSinfITS_r1_trimmSW5.20.min125.trimlog; done
for(i in seq_along(fnFs)) {
	j=paste0(fqpath, "/", sample.names[i],"_R1_TrimmSE_SW5.20_min125.fastq.gz")
	system2(TrimmomaticSE, 
		args = c(	fnFs[i], 		# input
  					j, 				# output
  					trimm.flags
  		),
  		stdout = paste0(sample.names[i],".TrimmSE_SW5.20_min125.out"), 
  		stderr = paste0(sample.names[i],".TrimmSE_SW5.20_min125.err")
  ) 
}

## ----identify/remove primers------------------------------------------------------------------------------------------------
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}


## ---- get all primer orients and save to objects ------------------------------------------------------------------------------------------------
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)


## ---- filter out any reads with Ns (necessary for dada function) -----------------------------------------------------------------------
fnFs.filtN <- file.path(fqpath, "filtN", paste0(sample.names, "_filtN_R1.fastq.gz")) # Put N-filterd files in filtN/ subdirectory
filterAndTrim(fnFs, fnFs.filtN, maxN = 0, multithread = TRUE,compress=T)


## ----check where primers are detected in reads------------------------------------------------------------------------------------------------
# may need to reverse complement and repeat if finding in the wrong places
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]))


## ----remove primers with cutadapt ------------------------------------------------------------------------------------------------
system2(cutadapt, args = "--version") # Run shell commands from R
# setup cutadapt output dir and output filenames
fqpath.cut <- file.path(fqpath, "cutadapt")
if(!dir.exists(fqpath.cut)) dir.create(fqpath.cut)
fnFs.cut <- file.path(fqpath.cut, paste0(sample.names, "_filtN_clip_R1.fastq"))
if(!dir.exists(fqpath.cut)) dir.create(fqpath.cut)


# setting reverse complements of primers
FWD.RC <- dada2:::rc(FWD)	
REV.RC <- dada2:::rc(REV)

# flags for cutadapt
R1.flags <- paste("-g", FWD, "-a", REV.RC) # Trim FWD and the reverse-complement of REV off of R1 (forward reads)
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, 
  	args = c(R1.flags, 
  				"-n", 2, # -n 2 required to remove FWD and REV from reads
                "-m", 1, # prevent zero-length seqs in output which mess up plotQualityProfile
                "-o", fnFs.cut[i],  # output files
                fnFs.filtN[i] # input files
                #,"-j 0" # auto-detect cores... apparently parallel not supported on python 2... :(
  	), 
  	stdout = paste0(fnFs.cut[i],".cutadapt.out"), 
  	stderr = paste0(fnFs.cut[i],".cutadapt.err")
  ) 
}

## ----check for no more primers detected after cutadapt------------------------------------------------------------------------------------
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]))

## ----setting clipped filenames and get sample names-------------------------------------------------------------------------------------
# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(fqpath.cut, pattern = "_R1.fastq$", full.names = TRUE))


## ----forward clipped qual plots----------------------------------------------------------------
pdf("cutFs.qualProf.pdf")
plotQualityProfile(cutFs)
dev.off()
# Error in density.default(qscore) : 'x' contains missing values
# error now prevented with `-m 1` in cutadapt (see https://github.com/benjjneb/dada2/issues/159)


## ----filter and trim --------------------------------------------------------------------------------------------------------
filtFs <- file.path(fqpath.cut, "filtered", paste0(sample.names, "_filtN_clip_filt_R1.fastq.gz"))	# end .gz because setting                                                                                          # compress=T in filterAndTrim here
names(filtFs) <- sample.names

# written to file and commented out so dont have to redo a bunch, just read in filterAndTrim.out
out <- filterAndTrim(cutFs, filtFs, maxN = 0, maxEE = 2, 
                      truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE
write.table(out, "filterAndTrim.out")
#out = as.matrix(read.table("filterAndTrim.out"))

head(out)  


## ----learn error rates and plot---------------------------------------------------------------------------------------------
# setting randomize=T bc i dont feel like the first 10 samples in order necessarily represent all the data
errF <- learnErrors(filtFs, multithread = TRUE, randomize = T)  

pdf("filt.plotErrors.pdf")          # can look at PDf if don't want to rerun errF/R
plotErrors(errF, nominalQ = TRUE)
dev.off()


## ----derep identical reads--------------------------------------------------------------------------------------------------
# this step is not included in the v1.16 dada2 16S tutorial... see this: https://github.com/benjjneb/dada2/issues/1095
derepFs <- derepFastq(filtFs, verbose = FALSE)
# Name the derep-class objects by the sample names
#names(derepFs) <- sample.names


## ----sample inference (dada function!)----------------------------------------------------------------------------------------
dadaFs <- dada(derepFs, err = errF, multithread = TRUE, verbose = FALSE)
# dadaFs[[1]]


## ----construct sequence table-----------------------------------------------------------------------------------------------------------
seqtab <- makeSequenceTable(dadaFs)
rownames(seqtab) = sample.names
dim(seqtab) # wow, 958 ASVs 
table(nchar(getSequences(seqtab)))  # 251 is most abundant size, avg larger than 16S also way more variable

## ----remove chimeras--------------------------------------------------------------------------------------------------------------------
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=FALSE)
sum(seqtab.nochim)/sum(seqtab)	# 0.9705774
dim(seqtab.nochim)  # ok, now 718 ASVs
table(nchar(getSequences(seqtab.nochim)))


## ----track reads through pipeline and plot---------------------------------------------------------------------------------------
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
head(track)
pdf("track.barplot.pdf",width=12,height=6)
barplot(t(track),beside=T,legend.text=T,las=2,cex.names=1,col=c("black","red","blue","green"),
		space=c(0.2,2),ylim=c(0,max(t(track))),args.legend = list(cex=0.6))
dev.off()

## ----assign taxonomy--------------------------------------------------------------------------------------------------------
# Assign taxonomy 
#unite.ref="/sda/data/UNITEdb/general/sh_general_release_dynamic_s_04.02.2020.fasta"	# moved to top of script
taxa=assignTaxonomy(seqtab.nochim,unite.ref,multithread = T,tryRC = T)
# inspect
taxa.print <- taxa  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
# looks like fungi!


## ----phyloseq-----------------------------------------------------------------------------------------------------------------------------
# phyloseq time
theme_set(theme_bw())

ps=phyloseq(otu_table(seqtab.nochim,taxa_are_rows = F), tax_table(taxa))

# if you want to write these tables to read into excel
write.csv(seqtab.nochim, file = "ps.otutable.csv")
write.csv(taxa, file= "ps.taxtable.csv")

# change reads counts to proportions
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))

top20 = names(sort(taxa_sums(ps), decreasing = T))[1:20]
# ps.prop.top20=transform_sample_counts(ps, function(OTU) OTU/sum(OTU))   # unneccesary bc the same as ps.prop
ps.prop.top20=prune_taxa(top20,ps.prop)
pdf("genus_stackedbar.pdf",width=12,height=6)
plot_bar(ps.prop.top20, fill="Genus", title = "Top 20 taxa")+theme(axis.text=element_text(size=14))
dev.off()

pdf("fam_stackedbar.pdf",width=12,height=6)
plot_bar(ps.prop.top20, fill="Family", title = "Top 20 taxa")+theme(axis.text=element_text(size=14))
dev.off()

## ----ITSx-----------------------------------------------------------------------------------------------------------------------------
sequences <- Biostrings::DNAStringSet(taxa_names(ps))
names(sequences) <- taxa_names(ps)
ps <- merge_phyloseq(ps, sequences)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
Biostrings::writeXStringSet(refseq(ps),"ps.ASV.fasta")

system2(ITSx, 
		args=c("-i","ps.ASV.fasta",		# input fasta
				"-o","ITSx_psASV"		# output prefix
			),
		stdout = "ITSx.out",
		stderr = "ITSx.err"
		)

## ----rarefy-----------------------------------------------------------------------------------------------------------------------------
# set.seed(14)	# set random seed to make this reproducible, otherwise youll get a diff random 10,000 ASVs each time
# raref.ps=rarefy_even_depth(ps, sample.size = 10000)
# raref.ps.prop <- transform_sample_counts(raref.ps, function(otu) otu/sum(otu))

## ----save image-----------------------------------------------------------------------------------------------------------------------------
saveRDS(ps,"ps.RDS")
save.image(paste0(format(Sys.time(), "%b%d%Y_%X%Z"),".RData"))
