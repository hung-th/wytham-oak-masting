library(biomaRt)
library(LEA)
library(tidyr)
library(ggplot2)
library(ggthemes)
library(forcats)
library(patchwork)
library(RColorBrewer)
library(vcfR)
library(qvalue)
library(data.table)
library(viridis)
library(ggrepel)
library(ggpubr)
library(qvalue)
library(clusterProfiler)
library(GenomicRanges)
library(rtracklayer)
library(enrichplot)
library(org.At.tair.db)
library(ggsci)
library(tidyverse)

dir.create("output")

#write.table(data[,-1], "lfmm_pheno.env", sep = "\t", quote = F, row.names = F, col.names = F)

#########################################################
# 1. Population structure using PCA
#########################################################

pca <- read.table("../1-bcftools/pca.eigenvec", header = F, sep = " ")
eigenval <- scan("../1-bcftools/pca.eigenval")

pca <- pca[,-1]
pca[,1] <- gsub(".bam", "", pca[,1])
names(pca)[1] <- "sample"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
pve <- data.frame(PC = 1:length(eigenval), pve = eigenval/sum(eigenval)*100, bs = brokenStick(1:length(eigenval), length(eigenval)) * 100)

## Plotting first three PCs
p <- ggplot(pca, aes(PC1, PC2, label = sample, col = sample)) + geom_point(size = 3) +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  theme_bw() +
  theme(legend.title = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 13),
		legend.position = "none") +
  scale_color_viridis(discrete = T) +
  geom_text_repel() +
  stat_ellipse()

ggsave("output/pop.pca12.png", width = 7, height = 7, units = "cm", dpi = 300, p)

p <- ggplot(pca, aes(PC2, PC3, label = sample, col = sample)) + geom_point(size = 3) +
  xlab(paste0("PC2 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC3 (", signif(pve$pve[2], 3), "%)")) +
  theme_bw() +
  theme(legend.title = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 13),
		legend.position = "none") +
  scale_color_viridis(discrete = T) +
  geom_text_repel() +
  stat_ellipse()

ggsave("output/pop.pca23.png", width = 7, height = 7, units = "cm", dpi = 300, p)

## Plotting variance explained by PCs
p <- ggplot(pve, aes(PC, pve, fill = PC)) + geom_bar(stat = "identity", show.legend = F) +
	geom_point(aes(PC, bs), size = 3) +
	geom_line(aes(PC, bs)) +
  xlab("PC") + ylab("Variance explained\n(%)") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        legend.position = "none") +
  scale_fill_viridis()

ggsave("output/pop.pc.var.png", width = 7, height = 7, units = "cm", dpi = 300, p)

#########################################################
# 2. Population structure using sNMF
#########################################################

## Subsampling for 100,000 SNPs
lfmm.file <- "../1-bcftools/WWO_37WGS.lfmm"
lfmm <- fread(lfmm.file, header = F, data.table = F)
lfmm_t <- t(lfmm)
set.seed(225)
lfmm_s <- t(lfmm_t[sample(nrow(lfmm_t), 100000), ])
fwrite(as.data.table(lfmm_s), file = "../1-bcftools/WWO_37WGS.subsampled.lfmm", sep = "\t", col.names = F)

## Running sNMF
project.snmf <- snmf("../1-bcftools/WWO_37WGS.subsampled.lfmm", K = 1:10, entropy = T, CPU = 256, repetitions = 3, project = "new")
#project.snmf <- load.snmfProject("../1-bcftools/WWO_37WGS.subsampled.snmfProject")  # If need to reload

## Plotting cross-entropy
ce <- sapply(1:max(project.snmf@K), function(k) {
  sapply(1:(length(project.snmf@K)/max(project.snmf@K)), function(rep) {
    cross.entropy(project.snmf, K = k, run = rep)
  })
})

ce.df <- as.data.frame(as.table(ce))
colnames(ce.df) <- c("K", "run", "cross.entropy")
ce.df$K <- rep(1:max(project.snmf@K), each = length(project.snmf@K)/max(project.snmf@K))
ce.df$run <- rep(1:(length(project.snmf@K)/max(project.snmf@K)), times = max(project.snmf@K))

p <- ggerrorplot(ce.df, "K", "cross.entropy", add = "mean", error.plot = "errorbar") +
	xlab("K") + ylab("Cross-entropy") +
	theme_bw() +
	theme(axis.title.x = element_text(size = 15),
		axis.title.y = element_text(size = 15),
		axis.text.x = element_text(size = 13),
		axis.text.y = element_text(size = 13),
		)
ggsave("output/snmf.ce.png", width = 7, height = 7, units = "cm", dpi = 300, p)

#########################################################
# 3. Missing genotype imputation
#########################################################

k <- 1
best <- which.min(cross.entropy(project.snmf, K = k))
impute(project.snmf, "../1-bcftools/WWO_37WGS.lfmm", method = 'random', K = k, run = best)

#########################################################
# 4. Preprocessing of variables
#########################################################

## Read data
data <- read.csv("../WWO_37WGS_data.csv", header = T)

## Transforming soil variables to PCs
soil.vars <- data[, c("SAND", "SILT", "CLAY")]
soil.pca <- prcomp(soil.vars, center = TRUE, scale. = TRUE)
soil.pcs <- soil.pca$x
data$SOIL_PC1 <- soil.pcs[,1]
data$SOIL_PC2 <- soil.pcs[,2]
pve <- data.frame(pve = (soil.pca$sdev^2) / sum(soil.pca$sdev^2) * 100)

data$SOIL_CLASS <- gsub("_", " ", data$SOIL_CLASS)

p <- ggplot(data, aes(SOIL_PC1, SOIL_PC2, color = SOIL_CLASS)) +
	geom_point(size = 3) +
	scale_color_npg() +
	xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +
	ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
	theme_bw() +
	theme(
		axis.title.x = element_text(size = 15),
		axis.title.y = element_text(size = 15),
		axis.text.x = element_text(size = 13),
		axis.text.y = element_text(size = 13),
		legend.title = element_blank(),
		text = element_text(size = 10)
	)

ggsave("output/soil.pca.png", width = 10, height = 7, units = "cm", dpi = 300, p)

## Retaining only a priori significant variables
data <- data[, c("GENOME", "TREE_ID", "VIZ_COUNT", "MATURE_ACORNS", "IMMAT_ACORNS", "ENLARGED_CUPS", "FLOWERS", "CANOPY_CLOSURE", "SPRING_PHENO", "SOIL_PC1", "MIDNOV_LAI")]

## Plotting all untransformed variables
myplots <- list()
p_idx <- 1
for (i in c(3:11)){
	
	p <- ggqqplot(data[,i]) +
	xlab("") +
	ylab(colnames(data[i])) +
	theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
		plot.title = element_text(size = 15))

	myplots[[p_idx]] <- p
	p_idx <- p_idx + 1

}

q <- ggarrange(plotlist = myplots, ncol = 3, nrow = 3)

ggsave("output/phenotypes_untrasnformed_qqplots.png", width = 20, height = 20, units = "cm", q)

## Log-transforming skewed variables
cols <- c("MATURE_ACORNS", "IMMAT_ACORNS")
data[cols] <- log1p(data[cols])

## Plotting transformed variables
myplots <- list()
p_idx <- 1
for (i in c(3:11)){
	
	p <- ggqqplot(data[,i]) +
	xlab("") +
	ylab(colnames(data[i])) +
	theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
		plot.title = element_text(size = 15))

	myplots[[p_idx]] <- p
	p_idx <- p_idx + 1

}

q <- ggarrange(plotlist = myplots, ncol = 3, nrow = 3)

ggsave("output/phenotypes_transformed_qqplots.png", width = 20, height = 20, units = "cm", q)

#########################################################
# 5. GWAS using LFMM2
#########################################################

Y <- as.matrix(fread(file = "../1-bcftools/WWO_37WGS.lfmm", sep = "\t")) #1577819
Y_no9 <- Y[, colSums(Y == 9) == 0] #8,602,814 SNPs  - 1304587

X <- as.matrix(data[, -c(1:2)])
lfmm2 <- lfmm2(input = Y_no9, env = X, K = 1)
save(lfmm2, file = "output/lfmm2.RData")
lfmm2.test <- lfmm2.test(object = lfmm2, input = Y_no9, env = X, genomic.control = T)
save(lfmm2.test, file = "output/lfmm2.test.RData")

loci <- fread(file = "../1-bcftools/WWO_37WGS.012.pos", sep = "\t")
colnames(loci) <- c("CHROM", "POS")
loci <- loci[colSums(Y == 9) == 0]
#loci$CHROM <- gsub("Qrob_Chr", "", loci$CHROM)

lfmm.p <- as.data.frame(cbind(loci, t(lfmm2.test$pvalues)))

lfmm.q <- lfmm.p
for (c in 3:ncol(lfmm.q)){
	p <- lfmm.q[, c]
	q <- qvalue(p)$qvalues
	lfmm.q[, c] <- q
}

lfmm.z <- cbind(loci, t(lfmm2.test$zscores))
lfmm.g <- as.data.frame(lfmm2.test$gif)

write.table(lfmm.p, "output/lfmm.p.tsv", col.names = T, row.names = F, sep = "\t")
write.table(lfmm.q, "output/lfmm.q.tsv", col.names = T, row.names = F, sep = "\t")
write.table(lfmm.z, "output/lfmm.z.tsv", col.names = T, row.names = F, sep = "\t")
write.table(lfmm.g, "output/lfmm.g.tsv", col.names = T, row.names = T, sep = "\t")

q_i <- which(rowSums(lfmm.q[,-c(1,2)] < 0.05) >= 1)
lfmm.sig <- q_i  #4,495 significant SNPs

lfmm.q_sig <- lfmm.q[lfmm.sig, ]
lfmm.z_sig <- lfmm.z[lfmm.sig, ]

write.table(lfmm.q_sig, "output/lfmm.q_sig.tsv", col.names = T, row.names = F, sep = "\t")
write.table(lfmm.z_sig, "output/lfmm.z_sig.tsv", col.names = T, row.names = F, sep = "\t")

#########################################################
# 6. GWAS visualisation
#########################################################

## Plotting P-value historgrams
lfmm.q_long <- gather(lfmm.q, trait, q, 3:ncol(lfmm.q), factor_key = TRUE)
lfmm.p_long <- gather(lfmm.p, trait, p, 3:ncol(lfmm.p), factor_key = TRUE)
lfmm.z_long <- gather(lfmm.z, trait, z, 3:ncol(lfmm.z), factor_key = TRUE)

lfmm.full <- cbind(lfmm.p_long, lfmm.q_long$q, lfmm.z_long$z)
colnames(lfmm.full) <- c("CHROM", "POS", "trait", "p", "q", "z")

p <- ggplot(lfmm.p_long, aes(x = p)) +
  geom_histogram(color="darkblue", fill="lightblue") + 
  #scale_x_continuous(breaks = seq(0, 100, 10)) +
  facet_wrap(~trait, ncol = 4) +
  theme_minimal() + labs(x = expression(paste(italic(P), "-value")), y = "Frequency") +
  scale_x_continuous(expand = c(0, 0)) +
  #scale_x_discrete(expand = expand_scale(add = 1)) +
  theme(
    panel.spacing.x = unit(1, "lines"),
    #axis.text.x = element_blank(),
	axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5, hjust=1),
	axis.text.y = element_text(size = 15),
    #panel.grid = element_blank(),
	#legend.position = "none",
	axis.title = element_text(size=18),
	strip.text.x = element_text(size=15),
	strip.placement = "outside",
	#panel.border = element_rect(colour = "black", fill = "transparent", size = 1)
  )

png("output/lfmm.p.histogram.png", width = 25, height = 25, res = 600, units = "cm")
p
dev.off()

## Plotting Q-value historgrams
p <- ggplot(lfmm.q_long, aes(x = q)) +
  geom_histogram(color="darkblue", fill="lightblue") + 
  #scale_x_continuous(breaks = seq(0, 100, 10)) +
  facet_wrap(~trait, ncol = 4) +
  theme_minimal() + labs(x = expression(paste(italic(Q), "-value")), y = "Frequency") +
  scale_x_continuous(expand = c(0, 0)) +
  #scale_x_discrete(expand = expand_scale(add = 1)) +
  theme(
    panel.spacing.x = unit(1, "lines"),
    #axis.text.x = element_blank(),
	axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5, hjust=1),
	axis.text.y = element_text(size = 15),
    #panel.grid = element_blank(),
	#legend.position = "none",
	axis.title = element_text(size=18),
	strip.text.x = element_text(size=15),
	strip.placement = "outside",
	#panel.border = element_rect(colour = "black", fill = "transparent", size = 1)
  )

png("output/lfmm.q.histogram.png", width = 25, height = 25, res = 600, units = "cm")
p
dev.off()

## Plotting Manhattan plots
traits <- unique(lfmm.q_long$trait)
lfmm.q_long_clean <- lfmm.q_long
lfmm.q_long_clean$CHROM <- gsub("Qrob_Chr", "", lfmm.q_long_clean$CHROM)
lfmm.q_long_clean$CHROM <- as.numeric(lfmm.q_long_clean$CHROM)
lfmm.q_long_clean <- lfmm.q_long_clean[!is.na(lfmm.q_long_clean$CHROM), ]

myplots <- list()

for (d in 1:length(traits)){
  
  tmp <- lfmm.q_long_clean[lfmm.q_long_clean$trait == traits[d], ]

  p <- ggplot(tmp, aes(POS/1000000, -log10(q))) +
    geom_point(color = ifelse(tmp$q < 0.1, "blue", "lightgrey")) + 
    scale_x_continuous(breaks = seq(0, 150, 10)) +
    facet_grid(~factor(CHROM), scales = "free", switch = "x", space = "free") +
    theme_minimal() + 
    labs(
      x = "SNP position (Mbp)", 
      y = expression(paste(-log[10], "(", italic(Q),"-",value,")"))
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    ggtitle(traits[d]) +
    theme(
      panel.spacing.x = unit(0.1, "lines"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
      legend.position = "none",
      plot.title = element_text(size=18),
      axis.title = element_text(size=18),
      strip.text.x = element_text(size = 15),
      strip.placement = "outside",
      panel.border = element_rect(colour = "black", fill = "transparent", size = 1)
    )

	png(paste0("output/lfmm.q.manhattan_", traits[d], ".png"), width = 30, height = 10, res = 600, units = "cm")
	print(p)
	dev.off()
}

#########################################################
# 7. GWAS annotations
#########################################################


cds <- import("/local/home/univ4423/oak/ref_Q_robur_v1/Qrob_cds.gff")
mcols(cds)$Parent <- as.character(unlist(mcols(cds)$Parent))

loci.gr <- GRanges(seqnames = loci$CHROM, ranges = IRanges(start = loci$POS, end = loci$POS), SNP_ID = paste0(loci$CHROM, ":", loci$POS))

hits <- findOverlaps(loci.gr, cds)
tmp <- cbind(loci[queryHits(hits), ], as.data.frame(cds[subjectHits(hits)]))
loci.annotated <- left_join(loci, tmp[, c("CHROM", "POS", "Parent", "ID", "strand")], by = c("CHROM", "POS")) %>% as.tibble()

qrob.annotation <- read.table("/local/home/univ4423/oak/ref_Q_robur_v1/Qrob_araport_pep_annotation.tsv", sep = "\t")
qrob.annotation$V1 <- gsub("Qrob_P", "Qrob_T", qrob.annotation$V1)
qrob.annotation <- qrob.annotation %>% filter(V11 < 1e-10) # 23,778 annotated
loci.annotated <- left_join(loci.annotated, qrob.annotation[,1:2], by = c("Parent" = "V1"))
colnames(loci.annotated) <- c("CHROM", "POS", "gene", "cds", "strand", "tair")

write.table(loci.annotated, "loci.annotated.tsv", sep="\t", quote=FALSE, row.names=FALSE)

lfmm.q_annotated <- merge(lfmm.q, loci.annotated, by = c("CHROM", "POS"))
lfmm.z_annotated <- merge(lfmm.z, loci.annotated, by = c("CHROM", "POS"))

write.table(lfmm.q, "output/lfmm.q_annotated.tsv", col.names = T, row.names = F, sep = "\t")
write.table(lfmm.z, "output/lfmm.z_annotated.tsv", col.names = T, row.names = F, sep = "\t")

lfmm.qz_long <- cbind(lfmm.q_long, lfmm.z_long$z)
lfmm.qz_long_sig <- lfmm.qz_long %>% filter(q < 0.05)
lfmm.qz_long_sig_annotated <- merge(lfmm.qz_long_sig, loci.annotated, by = c("CHROM", "POS"), sort = F)
lfmm.qz_long_sig_annotated_genes <- lfmm.qz_long_sig_annotated %>% filter(!is.na(gene))
write.table(lfmm.qz_long_sig_annotated_genes, "output/lfmm.qz_long_sig_annotated_genes.tsv", sep = "\t", row.names = F)

#########################################################
# 8. GWAS post-hoc analyses: VIZ_COUNT
#########################################################

#createKEGGdb::create_kegg_db("ath")
#install.packages("KEGG.db_1.0.tar.gz", repos = NULL, type = "source")

library(KEGG.db)

ek.df <- NULL
ego.df <- NULL

traits <- unique(lfmm.q_long$trait)

for (d in 1:length(traits)){

	tmp <- lfmm.q_long %>% filter(trait == traits[d])
	tmp$z <- lfmm.z_long %>% filter(trait == traits[d]) %>% pull(z)
	tmp <- tmp %>% filter(q < 0.05)
	tmp <- left_join(tmp, loci.annotated, by = c("CHROM", "POS")) %>% as.tibble

	gene.list <- gsub("\\..*", "", tmp$tair[!is.na(tmp$tair)])

	ek <- enrichKEGG(gene = gene.list, organism = "ath", qvalueCutoff = 0.05, use_internal_data = T)
	ego <- enrichGO(gene = gene.list, OrgDb = org.At.tair.db, keyType = "TAIR", ont = "ALL", qvalueCutoff = 0.05, readable = T)

	if(!is.null(ek)){
		ek.res <- ek@result %>% filter(p.adjust < 0.05)
	
		if(nrow(ek.res) > 0){
			ek.res$trait <- traits[d]
			ek.df <- rbind(ek.df, ek.res)
		}
	}
	
	if(!is.null(ego)){
		ego.res <- ego@result %>% filter(p.adjust < 0.05)
	
		if(nrow(ego.res) > 0){
			ego.res$trait <- traits[d]
			ego.df <- rbind(ego.df, ego.res)
		}
	}

}

write.table(ek.df, "ek.tsv", sep = "\t", row.names = F)
write.table(ego.df, "ego.tsv", sep = "\t", row.names = F)

#########################################################
# 9. Genotype-phenotype-environment
#########################################################

lfmm.q_annotated <- left_join(lfmm.q, loci.annotated, by = c("CHROM", "POS")) %>% as.tibble

comp1 <- c("VIZ_COUNT", "CANOPY_CLOSURE")
comp2 <- c("MATURE_ACORNS", "SPRING_PHENO")
comp3 <- c("ENLARGED_CUPS", "SOIL_PC1", "MIDNOV_LAI")

find_snp_overlap <- function(df, traits, threshold = 0.05) {
  sig_rows <- apply(df[traits], 1, function(row) all(row < threshold))
  df[sig_rows, ]
}

comp1.sig <- find_snp_overlap(lfmm.q_annotated, comp1)
comp2.sig <- find_snp_overlap(lfmm.q_annotated, comp2)
comp3.sig <- find_snp_overlap(lfmm.q_annotated, comp3)

snp_gr <- GRanges(
  seqnames = comp2.sig$CHROM,
  ranges = IRanges(start = comp2.sig$POS - 10000, end = comp2.sig$POS + 10000),
  strand = "*",
  SNP_ID = paste0(comp2.sig$CHROM, ":", comp2.sig$POS),
  SNP_POS = comp2.sig$POS
)

# Find CDS that fall within window
hits <- findOverlaps(snp_gr, cds, type = "any")
cds_in_window <- cds[subjectHits(hits)]
snp_matched <- snp_gr[queryHits(hits)]

# Collect matched CDS info
cds_df <- tibble(
  SNP_ID = snp_matched$SNP_ID,
  SNP_POS = snp_matched$SNP_POS,
  CHROM = as.character(seqnames(snp_matched)),
  CDS_start = start(cds_in_window),
  CDS_end = end(cds_in_window),
  CDS_strand = as.character(strand(cds_in_window)),
  Parent = cds_in_window$Parent
) %>%
  left_join(snp_info, by = "SNP_ID")

# Plot
p <- ggplot() +
  # CDS as rectangles
  geom_rect(data = cds_df,
            aes(xmin = CDS_start, xmax = CDS_end, ymin = 0.4, ymax = 0.6,
                fill = CDS_strand),
            alpha = 0.7, color = "black") +

  # Strand direction arrows
  geom_segment(data = cds_df,
               aes(
                 x = ifelse(CDS_strand == "+", CDS_start, CDS_end),
                 xend = ifelse(CDS_strand == "+", CDS_end, CDS_start),
                 y = 0.5, yend = 0.5
               ),
               arrow = arrow(length = unit(0.15, "cm")),
               color = "black", size = 0.5) +

  # Gene name (Parent) as label
  geom_text(data = cds_df,
            aes(x = (CDS_start + CDS_end)/2, y = 0.65, label = Parent),
            size = 3, angle = 45, hjust = 0) +

  # SNPs as red points
  geom_point(data = snp_info,
             aes(x = SNP_POS, y = 0.5),
             color = "red", size = 3) +

  # Facet by SNP ID, one row per plot
  facet_wrap(~ SNP_ID, scales = "free_x", ncol = 1) +

  # Theme and labels
  scale_fill_manual(values = c("+" = "forestgreen", "-" = "firebrick")) +
  labs(
    x = "Genomic position",
    y = "",
  ) +
  guides(fill = "none") +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

png(paste0("output/comp2.png"), width = 10, height = 20, res = 600, units = "cm")
print(p)
dev.off()
#########################

# Count number of sig SNPs for each variable

ensembl <- useMart("plants_mart", 
                   dataset = "athaliana_eg_gene", 
                   host = "https://plants.ensembl.org")

no_of_snps_d <- data.frame(trait = unique(lfmm.q_long$trait), adaptive_SNPs = NA, adaptive_SNPs_in_cds = NA, unique_genes = NA)

for (d in 1:length(unique(lfmm.q_long$trait))){
	tmp.q <- lfmm.q_long[lfmm.q_long$trait == unique(lfmm.q_long$trait)[d], ] ## using p-values
	tmp.z <- lfmm.z_long[lfmm.z_long$trait == unique(lfmm.z_long$trait)[d], ]$z
	tmp <- cbind(tmp.q, tmp.z)
	names(tmp)[4:5] <- c("q", "z")
	lfmm_d_sig <- tmp[tmp$q < 0.05,]
	lfmm_d_sig_in_cds <- merge(lfmm_d_sig, loci.annotated, by = c("CHROM", "POS"))
	lfmm_d_sig_in_cds <- lfmm_d_sig_in_cds %>% filter(!is.na(gene))
	tair_ids <- unique(lfmm_d_sig_in_cds$tair)
	tair_ids <- tair_ids[!is.na(tair_ids)]
	tair_ids <- gsub("\\.\\d+$", "", tair_ids)  # remove version suffix
	annot <- getBM(attributes = c("tair_symbol", "external_gene_name", "description"),
               filters = "tair_symbol",
               values = tair_ids,
               mart = ensembl)
	lfmm_d_sig_in_cds$tair_base <- gsub("\\.\\d+$", "", lfmm_d_sig_in_cds$tair)
	annotated <- merge(lfmm_d_sig_in_cds, annot, 
                   by.x = "tair_base", by.y = "tair_symbol", 
                   all.x = TRUE)
	
	no_of_snps_d[d, ]$adaptive_SNPs <- nrow(lfmm_d_sig)
	no_of_snps_d[d, ]$adaptive_SNPs_in_cds <- nrow(lfmm_d_sig_in_cds)
	no_of_snps_d[d, ]$unique_genes <- length(unique(lfmm_d_sig_in_cds$gene))
}

write.csv(no_of_snps_d, "output/no_of_snps_d.csv", row.names = F)


# to prepare a master list of SNPs with annotations

cds_annotation <- read.table("/local/home/univ4423/oak/ref/Qrob_cds.bed", sep = "\t")

lfmm_in_cds <- merge(lfmm.full, lfmm_all_annotated, by.x = c("CHROM", "POS"), by.y = c("chr", "pos"))
lfmm_in_cds <- merge(lfmm_in_cds, cds_annotation[,2:4], by.x = "cds", by.y = "V4", all.x = T)
lfmm_in_cds_final <- lfmm_in_cds[,c(4,2,3,5,6,7,1,13,14,8,12)]
lfmm_in_cds_final <- lfmm_in_cds[,c(4,2,3,5,6,1,12,13,7,11)]
names(lfmm_in_cds_final) <- c("trait", "CHROM", "POS", "p", "q", "z", "cds", "cds_start", "cds_end", "gene", "tair")
lfmm_in_cds_final <- lfmm_in_cds_final[with(lfmm_in_cds_final, order(trait, CHROM, POS)), ]

write.table(lfmm_in_cds_final, "output/lfmm_full_results.tsv", sep = "\t", quote = F, row.names = F)

lfmm_d_sig <- tmp[tmp$p < 0.05,]




##### To only look at COV

lfmm.q_cov <- lfmm.q[, "CROP_PERCENT"]
lfmm.z_cov <- lfmm.z[, "CROP_PERCENT"]
lfmm_cov <- cbind(loci, lfmm.q_cov, lfmm.z_cov)
names(lfmm_cov)[3:4] <- c("q", "z")
lfmm_cov_sig <- lfmm_cov[lfmm_cov$q < 0.05]
#lfmm_cov_sig <- lfmm_cov[lfmm_cov$q < 0.05 & abs(lfmm_cov$z) > 2,]

lfmm_cov_bed <- lfmm_cov_sig[,1:2]
lfmm_cov_bed$end <- lfmm_cov_sig[,2] + 1

write.table(lfmm_cov_bed, "output/lfmm_cov_sig.bed", row.names = F, col.names = F, quote = F, sep = "\t")

lfmm_cov_annotated <- read.table("output/lfmm_cov_sig_intersect.tsv", sep = "\t")

qrob_annotation <- read.table("/local/home/univ4423/oak/ref/Qrob_araport_pep_annotation.tsv", sep = "\t")
qrob_annotation$V1 <- gsub("Qrob_P", "Qrob_T", qrob_annotation$V1)
lfmm_cov_annotated <- merge(lfmm_cov_annotated, qrob_annotation[,1:2], by.x = "V4", by.y = "V1")

save(lfmm_cov_annotated, file = "output/post-gwas/lfmm_cov_annotated.RData")

lfmm_cov_ek <- enrichKEGG(gene = gsub("\\..*", "", lfmm_cov_annotated$V2.y),
                 organism     = 'ath',
                 pvalueCutoff = 0.05)
				 
lfmm_cov_mk <- gseKEGG(geneList = gsub("\\..*", "", lfmm_cov_annotated$V2.y),
               organism     = 'hsa',
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)

#########################
tmp2 <- lfmm.z_long[lfmm.z_long$trait == unique(lfmm.z_long$trait)[d], ]
	tmp2$CHROM <- as.numeric(tmp2$CHROM)
	tmp2 <- tmp2[!is.na(tmp2$CHROM), ]


 head(tmp[-log10(tmp$p) > 20, ])
 -log10(getmode(tmp$p))
 
 head(tmp2[-log10(tmp$p) > 20,])

########################

## RDA
rda <- rda(as.data.frame(Y) ~ ., data = as.data.frame(X), scale = T)

RsquareAdj(rda)
$r.squared
[1] 0.3789006
$adj.r.squared
[1] 0.0008400423 #Most of the SNPs are neutral

summary(rda)$concont

Importance of components:
                           RDA1      RDA2      RDA3      RDA4      RDA5
Eigenvalue            1.755e+05 1.702e+05 1.672e+05 1.668e+05 1.642e+05
Proportion Explained  7.799e-02 7.563e-02 7.432e-02 7.413e-02 7.298e-02
Cumulative Proportion 7.799e-02 1.536e-01 2.279e-01 3.021e-01 3.751e-01
                           RDA6      RDA7      RDA8      RDA9     RDA10
Eigenvalue            1.639e+05 1.632e+05 1.613e+05 1.597e+05 1.574e+05
Proportion Explained  7.285e-02 7.254e-02 7.168e-02 7.097e-02 6.997e-02
Cumulative Proportion 4.479e-01 5.205e-01 5.921e-01 6.631e-01 7.331e-01
                          RDA11     RDA12     RDA13     RDA14
Eigenvalue            1.560e+05 1.543e+05 1.508e+05 1.394e+05
Proportion Explained  6.934e-02 6.860e-02 6.702e-02 6.197e-02
Cumulative Proportion 8.024e-01 8.710e-01 9.380e-01 1.000e+00

rda.concont <- as.data.frame(summary(rda)$concont$importance)
rda.prop <- as.data.frame(t(rda.concont[2,]))
rda.prop[,1] <- rda.prop[,1] * 100
names(rda.prop)[1] <- "Proportion"
rda.prop$RDA <- as.numeric(gsub("RDA", "", row.names(rda.prop)))

p <- ggplot(rda.prop, aes(RDA, Proportion, fill = RDA)) + geom_bar(stat = "identity", show.legend = F) +
	#geom_point(aes(PC, bs), size = 3) +
  xlab("RDA") + ylab("Variance explained (%)") +
  theme_bw() +
  theme(legend.title = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        legend.position = "none") +
  scale_fill_viridis()

ggsave("output/rda.var.png", width = 10, height = 10, units = "cm", dpi = 300, p)

#####

library(corrplot)

res <- cor(X)

png("output/corr.png", width = 17, height = 17, res = 600, units = "cm")
corrplot(res, type = "upper", order = "hclust", tl.col = "black", tl.srt = 45)
dev.off()

anova.cca(rda, by="axis")
Permutation test for rda under reduced model
Forward tests for axes
Permutation: free
Number of permutations: 999

Model: rda(formula = as.data.frame(Y) ~ SPRING_PHENO + MIDAPRIL_LAI + AUTUMN_PHENO + MIDNOV_LAI + MIDDEC_LAI + AUPPC + VIZ_COUNT + FLOWERS + ABORTED_ACORNS + IMMAT_ACORNS + MATURE_ACORNS + COV + CROP_PERCENT + CAT, data = as.data.frame(X), scale = T)
         Df Variance      F Pr(>F)
RDA1      1   175471 1.0943  0.464
RDA2      1   170162 1.0612  1.000
RDA3      1   167207 1.0428  1.000
RDA4      1   166790 1.0402  1.000
RDA5      1   164187 1.0239  1.000
RDA6      1   163909 1.0222  1.000
RDA7      1   163214 1.0179  1.000
RDA8      1   161262 1.0057  1.000
RDA9      1   159672 0.9958  1.000
RDA10     1   157411 0.9817  1.000
RDA11     1   155993 0.9728  1.000
RDA12     1   154347 0.9626  0.990
RDA13     1   150789 0.9404  0.926
RDA14     1   139428 0.8695  0.754
Residual 23  3687974

vif.cca(rda)
  SPRING_PHENO   MIDAPRIL_LAI   AUTUMN_PHENO     MIDNOV_LAI     MIDDEC_LAI
      5.336914       4.290131      64.288429      29.846459      11.697352
         AUPPC      VIZ_COUNT        FLOWERS ABORTED_ACORNS   IMMAT_ACORNS
      6.083166       7.609899       2.978488       3.575444       4.413398
 MATURE_ACORNS            COV   CROP_PERCENT            CAT
      5.220568       1.752450       5.657333       7.185536

collinear_vars <- c("AUTUMN_PHENO", "MIDNOV_LAI", "MIDDEC_LAI")
Z <- X[, !colnames(X) %in% collinear_vars]

rda2 <- rda(as.data.frame(Y) ~ ., data = as.data.frame(Z), scale = T)
vif.cca(rda2)
  SPRING_PHENO   MIDAPRIL_LAI          AUPPC      VIZ_COUNT        FLOWERS
      3.891355       3.721162       1.575964       7.578599       2.620908
ABORTED_ACORNS   IMMAT_ACORNS  MATURE_ACORNS            COV   CROP_PERCENT
      2.754256       4.198703       4.690946       1.587877       4.382698
           CAT
      6.002184

collinear_vars <- c("MIDNOV_LAI", "MIDDEC_LAI")
Z <- X[, !colnames(X) %in% collinear_vars]

rda3 <- rda(as.data.frame(Y) ~ ., data = as.data.frame(Z), scale = T)
vif.cca(rda3)
  SPRING_PHENO   MIDAPRIL_LAI   AUTUMN_PHENO          AUPPC      VIZ_COUNT
      5.123729       3.946352       6.811459       5.181601       7.603146
       FLOWERS ABORTED_ACORNS   IMMAT_ACORNS  MATURE_ACORNS            COV
      2.876482       3.194543       4.210305       5.058427       1.587970
  CROP_PERCENT            CAT
      4.797237       6.135911


RsquareAdj(rda3)
$r.squared
[1] 0.3789006
$adj.r.squared
[1] 0.0008400423 #Most of the SNPs are neutral

summary(rda3)$concont

Importance of components:
                           RDA1      RDA2      RDA3      RDA4      RDA5
Eigenvalue            1.732e+05 1.691e+05 1.656e+05 1.649e+05 1.634e+05
Proportion Explained  8.994e-02 8.784e-02 8.598e-02 8.566e-02 8.487e-02
Cumulative Proportion 8.994e-02 1.778e-01 2.638e-01 3.494e-01 4.343e-01
                           RDA6      RDA7      RDA8      RDA9     RDA10
Eigenvalue            1.628e+05 1.603e+05 1.597e+05 1.576e+05 1.562e+05
Proportion Explained  8.456e-02 8.323e-02 8.295e-02 8.185e-02 8.112e-02
Cumulative Proportion 5.189e-01 6.021e-01 6.850e-01 7.669e-01 8.480e-01
                          RDA11     RDA12
Eigenvalue            1.522e+05 1.405e+05
Proportion Explained  7.903e-02 7.296e-02
Cumulative Proportion 9.270e-01 1.000e+00

rda3.concont <- as.data.frame(summary(rda3)$concont$importance)
rda3.prop <- as.data.frame(t(rda3.concont[2,]))
rda3.prop[,1] <- rda3.prop[,1] * 100
names(rda3.prop)[1] <- "Proportion"
rda3.prop$RDA <- as.numeric(gsub("RDA", "", row.names(rda3.prop)))

p <- ggplot(rda3.prop, aes(RDA, Proportion, fill = RDA)) + geom_bar(stat = "identity", show.legend = F) +
	#geom_point(aes(PC, bs), size = 3) +
  xlab("RDA") + ylab("Variance explained (%)") +
  theme_bw() +
  theme(legend.title = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        legend.position = "none") +
  scale_fill_viridis()

ggsave("output/rda3.var.png", width = 10, height = 10, units = "cm", dpi = 300, p)

res3 <- cor(Z)

png("output/corr3.png", width = 17, height = 17, res = 600, units = "cm")
corrplot(res3, type = "upper", order = "hclust", tl.col = "black", tl.srt = 45)
dev.off()

signif.full <- anova.cca(rda3, parallel=getOption("mc.cores"))