library(RColorBrewer)
library(ape)
library(gplots)
# run analysis.r on mlgf_on_refs
ref_hmmer_clade_counts<-read.table("results/mlgf_on_refs/clade_mean_sizes.tsv", sep="\t")
ref_hmmer_clade_counts<-ref_hmmer_clade_counts[!apply(is.na(ref_hmmer_clade_counts), 1, any),]
# run analysis.r on mlgf_on_OneKP_prot_cdhit
onekp_hmmer_clade_counts<-read.table("results/mlgf_on_OneKP_prot_cdhit99/clade_mean_sizes.tsv", sep="\t") #hmmer_clade_counts

# sample counts
ref_hmmer_sample_counts <- read.table("results/mlgf_on_refs/sample_count_table.tsv", sep="\t")
onekp_hmmer_sample_counts <- read.table("results/mlgf_on_OneKP_prot_cdhit99/sample_count_table.tsv", sep="\t")

x<-as.vector(as.matrix(ref_hmmer_clade_counts))
y<-as.vector(as.matrix(onekp_hmmer_clade_counts[dimnames(ref_hmmer_clade_counts)[[1]], dimnames(ref_hmmer_clade_counts)[[2]]]))

cor(x, y, method="pearson")
# 0.9469704
cor(x, y, method="spearman")
# 0.9113313
mean(x[x!=0&y!=0]/y[x!=0&y!=0])
# 2.303616

pdf("plots/figure_S1a genomes vs RNA-Seq Inferred gene family sizes per clade.pdf", width=89/25.4, height=89/25.4)
	par(cex=0.5, mar=c(4, 4, 2, 2) + 0.1, bty="l", las=1)
	
	plot(x, y, xlab="from genomes", ylab="from OneKP", pch=20, ylim=c(1, 5e+3), xlim=c(1,5e+3), log="xy", main="Average Gene Family Sizes per Clade\n(for 23 gene families and 13 clades)")
	lines(c(1,5e+3),c(1,5e+3), col="grey")
	abline(lm(y ~ x - 1), untf=T)
	text(4, 3000, labels=paste0("Pearson = ", round(cor(x, y, method="pearson"), digits=2), "\nSpearman = ", round(cor(x, y, method="spearman"), digits=2)))
dev.off()

# control for assemlby/genome sizes
clade_col<-c("#654522", "#872c17", "#ddd300", "#b3446c", "#f6a601", "#604e97", "#f99379", "#0067a5", "#e68fac", "#018856", "#848482", "#c3b280", "#f38401", "#a1caf1", "#be0032", "#875692", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#909200")
names(clade_col)<-c("Chromista", "Rhodophyta", "Glaucophyta", "Chlorophyta", "Streptophyte Algae", "Hornworts", "Liverworts", "Mosses", "Lycophytes", "Monilophytes", "Gymnosperms", "Basal Angiosperms", "Monocots", "Chloranthales", "Magnoliids", "Basal Eudicots", "Asterids", "Caryophyllales", "Rosids", "Saxifragales", "Santalales")

tree<-read.tree("results/mlgf_on_refs/tree.newick")
tmp<-read.table("data/phylogeny/clade_dictionary_edit.tsv", row.names=1, sep="\t")
clade_dictionary<-as.character(tmp[,2])
names(clade_dictionary)<-dimnames(tmp)[[1]]
rm(tmp)

clades<-dimnames(ref_hmmer_clade_counts)[[1]]
ratios<-(onekp_hmmer_clade_counts[dimnames(ref_hmmer_clade_counts)[[1]], dimnames(ref_hmmer_clade_counts)[[2]]]+1)/(ref_hmmer_clade_counts+1)

pdf("plots/figure_S1b assembly vs genome sizes per clade.pdf", width=89/25.4, height=89/25.4)
	par(mar=c(2, 4, 2, 2) + 0.1, bty="n", xaxt="n", cex=0.5, las=1)

	boxplot2(t(ratios)[,clades], log="y", ylim=c(0.1, 10), ylab="Gene family size ratio: OneKP / genomes", names= F, col=clade_col[clades], pch=20)
	legend("topright", clades, fill=clade_col[clades], ncol=2, bty="n", cex=0.8)
dev.off()

pdf("plots/figure_S1c size ratio vs mean size.pdf", width=89/25.4, height=89/25.4)
	x<-as.vector(as.matrix((ref_hmmer_clade_counts+onekp_hmmer_clade_counts[dimnames(ref_hmmer_clade_counts)[[1]], dimnames(ref_hmmer_clade_counts)[[2]]])/2))
	y<-as.vector(as.matrix(ratios))
	par(mar=c(2, 4, 2, 2) + 0.1, cex=0.5, las=1)
	plot(x, y, log="y", ylim=c(0.1, 10), ylab="Gene family size ratio: OneKP / genomes", names= F, xlab="mean size", pch=20)
dev.off()
