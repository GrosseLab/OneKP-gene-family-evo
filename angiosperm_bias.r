# run pipeline on mlgf_on_OneKP_prot_cdhit99
pfamA_sample<-read.table("results/mlgf_on_OneKP_prot_cdhit99/sample_count_table.tsv", sep="\t")
pfamA_clade<-read.table("results/mlgf_on_OneKP_prot_cdhit99/clade_mean_sizes.tsv", sep="\t")
# run pipeline on Pfam-A_noStrepto_on_OneKP_prot_cdhit99
pfamA_noStrepto_sample<-read.table("results/Pfam-A_noStrepto_on_OneKP_prot_cdhit99/sample_count_table.tsv", sep="\t")
pfamA_noStrepto_clade<-read.table("results/Pfam-A_noStrepto_on_OneKP_prot_cdhit99/clade_mean_sizes.tsv", sep="\t")

include<-c("F-box", "AP2_EREBP", "bHLH", "bZIP", "C2H2", "C3H", "WRKY", "ABC", "PP2C", "Receptor-like Kinases", "P450s", "NAC")#, "GRAS")

pfamA_sample<-pfamA_sample[include,]
pfamA_clade<-pfamA_clade[,gsub("[ -]", '.', include)]
pfamA_noStrepto_sample<-pfamA_noStrepto_sample[include,]
pfamA_noStrepto_clade<-pfamA_noStrepto_clade[,gsub("[ -]", '.', include)]

# Korrelationen
x <- as.numeric(as.matrix(pfamA_sample))
y <- as.numeric(as.matrix(pfamA_noStrepto_sample))
cor(x, y, method="pearson")
# 0.9809391
cor(x, y, method="spearman")
# 0.8633884
pdf("plots/figure_4 PfamA vs PfamA noStrepto samples.pdf", width=89/25.4, height=89/25.4)
	par(cex=0.5, mar=c(4, 4, 2, 2) + 0.1, bty="l", las=1)

	smoothScatter(as.vector(as.matrix(pfamA_sample)), as.vector(as.matrix(pfamA_noStrepto_sample)), xlab="PfamA", ylab="PfamA w/o Streptophytes")#, log="xy")#, main="per sample gene family sizes")
	lines(c(1,5e+3),c(1,5e+3), col="grey")
	text(250, 1500, labels=paste0("Pearson = ", round(cor(x, y, method="pearson"), digits=2), "\nSpearman = ", round(cor(x, y, method="spearman"), digits=2)))
dev.off()

x <- as.numeric(as.matrix(pfamA_clade))
y <- as.numeric(as.matrix(pfamA_noStrepto_clade))
cor(x, y, method="pearson")
# 0.9800696
cor(x, y, method="spearman")
# 0.8864
pdf("plots/figure_4 PfamA vs PfamA noStrepto clades.pdf", width=89/25.4, height=89/25.4)
	par(cex=0.5, mar=c(4, 4, 2, 2) + 0.1, bty="l", las=1)
	plot(c(1, max(pfamA_clade, pfamA_noStrepto_clade)),c(1, max(pfamA_clade, pfamA_noStrepto_clade)), col="grey", t="l", xlab="PfamA", ylab="PfamA w/o Streptophytes", main="Average Gene Family Sizes per Clade (for 15 gene families)")
	points(as.vector(as.matrix(pfamA_clade)), as.vector(as.matrix(pfamA_noStrepto_clade)), pch=20)
	lines(c(1,5e+3),c(1,5e+3), col="grey")
	text(80, 650, labels=paste0("Pearson = ", round(cor(x, y, method="pearson"), digits=2), "\nSpearman = ", round(cor(x, y, method="spearman"), digits=2)))
dev.off()
