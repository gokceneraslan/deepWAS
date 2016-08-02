library(GenomicRanges)
library(GenomicFeatures)
library(dplyr)
library(ggplot2)

states <- readRDS('chromhmm_states.Rds')
codes <- read.delim('epigenomes.tsv', header = T, stringsAsFactors = F)

hmm <- readRDS('chromHMM_all_states.Rds')
hmm$group <- codes$GROUP[match(hmm$code, codes$EID)]

#group.hmm <- split(hmm, hmm$name)
# hmm.reduced <- lapply(group.hmm, function(gr){unlist(reduce(split(gr, gr$group)))})
# hmm.reduced <- unlist(GRangesList(hmm.reduced))
# names <- simplify2array(strsplit(names(hmm.reduced), '.', fixed = T))
# hmm.reduced$name <- names[1,]
# hmm.reduced$group <- names[2,]


hits <- readRDS('all_hits_gr.Rds')

lasso.coef <- readRDS('../results_fixedTF_1540ind/lasso_fits_17M_1500inds_coef.Rds')
good.lasso.coef <- lasso.coef[sapply(lasso.coef, length) > 0]
gmhits <- unique(unname(unlist(good.lasso.coef[ grep('GM12878', names(good.lasso.coef))])))
gmhits <- gmhits[substr(gmhits, 1, 2) == 'rs']
gmhits.gr <- hits[hits$rsid %in% gmhits]

hmm.overlap <- findOverlaps(hits, hmm)
hmm.overlap.df <- as.data.frame(mcols(hmm[subjectHits(hmm.overlap)]))
hmm.overlap.df$rsid <- hits$rsid[queryHits(hmm.overlap)]

hmm.overlap.df$tissue <- codes$Standardized.Epigenome.name[match(hmm.overlap.df$code, codes$EID)]
hmm.overlap.df$group <- codes$GROUP[match(hmm.overlap.df$code, codes$EID)]
hmm.overlap.df$chromatin <- states$DESCRIPTION[match(hmm.overlap.df$name, states$STATE)]


interesting.codes <- c('1_TssA', '2_TssAFlnk', '6_EnhG', '7_Enh', '10_TssBiv', '11_BivFlnk', '12_EnhBiv')
uninteresting.codes <- c('13_ReprPC', '14_ReprPCWk', '15_Quies', '5_TxWk')
uninteresting.tissues <- c("ENCODE2012", 'Other', 'Mesench', 'IMR90')

df <- subset(hmm.overlap.df,  !(name %in% uninteresting.codes) & !(group %in% uninteresting.tissues))


# try with dplyr ----------------------------------------------------------

tbl_df(hmm.overlap.df) %>%
  filter(!(name %in% uninteresting.codes) & !(group %in% uninteresting.tissues)) %>%
  group_by(group, chromatin) %>%
  summarise(count=n_distinct(rsid)) %>%
  ggplot(aes(x=group, y=count, fill=chromatin)) +
  labs(x='Roadmap Tissue Group', y='Number of overlapping SNPs') +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_brewer(name='Chromatin state', palette = 'Paired') + theme_bw()
  viridis::scale_fill_viridis(discrete=T)

ggplot(df, aes(x=chromatin, fill=chromatin)) +
  geom_bar() +
  scale_fill_brewer(palette = 'Set1') +
  facet_wrap(~group) +
  coord_flip() +
  theme_minimal()

ggplot(df, aes(x=chromatin, fill=group)) +
  geom_bar(position='dodge') +
  scale_fill_brewer(palette = 'Paired') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


# plot blood and brain separately -----------------------------------------

group.hmm <- split(hmm, hmm$group)

bb.hmm <- hmm[hmm$group %in% c('Blood & T-cell', 'Brain')]
hmm.overlap <- findOverlaps(gmhits.gr, bb.hmm)
hmm.overlap.df <- as.data.frame(mcols(bb.hmm[subjectHits(hmm.overlap)]))
hmm.overlap.df$rsid <- gmhits.gr$rsid[queryHits(hmm.overlap)]

hmm.overlap.df$tissue <- codes$Standardized.Epigenome.name[match(hmm.overlap.df$code, codes$EID)]
hmm.overlap.df$group <- codes$GROUP[match(hmm.overlap.df$code, codes$EID)]
hmm.overlap.df$chromatin <- states$DESCRIPTION[match(hmm.overlap.df$name, states$STATE)]

tbl_df(hmm.overlap.df) %>%
  filter(!(name %in% uninteresting.codes)) %>%
  group_by(tissue, chromatin) %>%
  mutate(count=n_distinct(rsid)) %>%
  select(-rsid) %>%
  ggplot(aes(x=tissue, y=count, fill=chromatin)) +
  labs(x='Roadmap Tissue Group', y='Number of overlapping SNPs') +
  geom_bar(stat='identity') +
  facet_wrap(~group, scales = 'free') +
  scale_fill_brewer(name='Chromatin state', palette = 'Paired') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
  #viridis::scale_fill_viridis(discrete = T)

df <- subset(hmm.overlap.df,  !(name %in% uninteresting.codes))
ggplot(df, aes(x=tissue, fill=chromatin)) +
  geom_bar() +
  scale_fill_brewer(palette = 'Set1') +
  facet_wrap(~group, scales='free_x') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



# Check overlap for each epignome with all 177 SNPs --------------

epigenome.hmm <- split(hmm, hmm$code)

final.percentages <- lapply(codes$EID, function(epi) {

print(epi)
x <- tbl_df(data.frame(w=width(epigenome.hmm[[epi]]), 
                       state=epigenome.hmm[[epi]]$name, 
                       stringsAsFactors = F)) %>%
  group_by(state) %>% 
  summarise(all.states=sum(as.numeric(w)))

overlaps <- findOverlaps(hits, epigenome.hmm[[epi]])
overlaps.states <- epigenome.hmm[[epi]][overlaps@subjectHits]

y <- tbl_df(data.frame(w=width(overlaps.states), 
                       state=overlaps.states$name, 
                       stringsAsFactors = F)) %>%
  group_by(state) %>% 
  summarise(overlap.regions=sum(as.numeric(w)))

z <- merge(x,y)
data.frame(state=z$state, 
           percentage=(z$overlap.regions/z$all.states) * 100, 
           snp.regions=z$overlap.regions,
           all.regions=z$all.states,
           epi=epi, 
           stringsAsFactors = F)
})

names(final.percentages) <- codes$EID

f <- do.call(rbind, final.percentages)

f$mnemonic <- codes$Epigenome.Mnemonic[ match(f$epi, codes$EID)]
f$tissue <-   codes$Standardized.Epigenome.name[ match(f$epi, codes$EID)]
f$group <- codes$GROUP[ match(f$epi, codes$EID)]
f$state.full <- states$DESCRIPTION [ match(f$state, states$STATE)]

f <- subset(f, !(state %in% uninteresting.codes) & !(group %in% uninteresting.tissues))

ggplot(f, aes(x=epi, y=percentage, fill=group)) +
  geom_bar(stat = 'identity') +  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~state.full, scales = 'free') + labs(x='Epigenome', y='Percent of chromatin state overlapping with DeepQTLs')

ggsave('roadmap_overlap_results_only27.pdf', width = 18, height = 11)


ggplot(subset(f, group=='Blood & T-cell'), aes(x=epi, y=percentage, fill=state.full)) +
  geom_bar(stat = 'identity') +  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~state.full, scales = 'free') + labs(x='Epigenome', y='Percent of chromatin state overlapping with DeepQTLs')


ggplot(f, aes(x=epi, y=all.regions, color=state.full)) + geom_point() + facet_wrap(~group, scales = 'free_x') + theme_minimal()

ggplot(f, aes(x=epi, y=snp.regions, color=state.full)) + geom_point() + facet_wrap(~group, scales = 'free_x') + theme_minimal()



# upset plot --------------------------------------------------------------

library(ChIPseeker)

#txdb <- makeTxDbFromUCSC(genome='hg19', tablename='refGene')
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

peakAnno = annotatePeak(hits, TxDb = txdb)#, annoDb="org.Hs.eg.db")
peakAnno
pdf('../Manuscript/Figures/upset177_pie.pdf', 12, 8)
plotAnnoPie(peakAnno)
dev.off()

upsetplot(peakAnno)
hits

pdf('../Manuscript/Figures/upset177.pdf', 12, 8)
upsetplot(peakAnno)
dev.off()

gmhits.gr <- hits[hits$rsid %in% gmhits]
peakAnnoGM <- annotatePeak(gmhits.gr, TxDb = txdb, annoDb = "org.Hs.eg.db")
plotAnnoPie(peakAnnoGM)
upsetplot(peakAnnoGM)

