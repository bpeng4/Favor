setwd("/work/benson/bpeng4/Favor/Favor_16S")
suppressMessages({
library("tidyverse")
library("dada2")
library("gridExtra")
library("devtools")
})
no_of_cores = 36
metadata = "Meta_Favor.csv"
metadata_df = read.csv(metadata)
fnFs <- metadata_df$fq1

sample.names <- metadata_df$Sample_Name 

filt_path <- paste0(getwd(), '/filtered') # don't change this 
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))

FORWARD_TRUNC <- 250 # determine from quality plots

out <- filterAndTrim(fnFs, filtFs,
                     truncLen=FORWARD_TRUNC, 
                     trimLeft=20, maxEE=2, 
                     multithread=FALSE,
                     matchIDs=TRUE, compress=TRUE, 
                     verbose=TRUE)

derepFs <- derepFastq(filtFs, n = 1e+06, verbose = TRUE)

names(derepFs) <- sample.names

errF <- learnErrors(filtFs, verbose=TRUE, multithread=no_of_cores)

dadaFs <- dada(derepFs, err=errF, pool=TRUE, multithread=no_of_cores, 
               verbose=TRUE)

seqtab <- makeSequenceTable(dadaFs)

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus",
                                    multithread=no_of_cores, verbose=TRUE)

#View(seqtab.nochim)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x)) # getUniques() gets abundance of unique sequences

# Calculate the number of reads at different steps
track <- cbind(out, 
               sapply(dadaFs, getN), 
               rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
track

track.long <- as.data.frame(track, 
                            row.names = row.names(track)) %>% 
  gather(., key = steps, value = counts, input:nonchim, factor_key = TRUE)

Datatracking<-ggplot(track.long, aes(x = steps, y = counts, color = steps)) +
  theme_classic() +
  geom_boxplot() +
  geom_jitter(shape = 16, position = position_jitter(0.3)) +
  scale_y_continuous(labels = scales::comma)
ggsave("/work/benson/bpeng4/Favor/Plots/Datatracking_Forward.png", plot = Datatracking, width = 6, height = 4, dpi = 300)



seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", 
                                    multithread=no_of_cores, verbose=TRUE)

save(metadata_df, seqtab.nochim, file = "./intermediate/FavorForward.rda")
