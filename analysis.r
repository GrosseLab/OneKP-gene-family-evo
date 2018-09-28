library("data.table")
library("ape")
library("geiger")
library("RColorBrewer")
library("gplots")

# reset everything
rm(list=ls())
relative<-NULL
completeness<-NULL

################################
# BEGIN: uncomment desired study
################################

# # study 1: large gene families on OneKP data set
# case_id<-"mlgf_on_OneKP_prot_cdhit99"
# seq_to_species_function <- function(s){return(substr(s, 1, 4))}
# # relative<-"data/sequence_counts_prot_cdhit_99.txt"
# completeness<-"data/CEGMA_BUSCO_100-missing.tsv"

# # study 2: large gene families on reference genomes
# case_id<-"mlgf_on_refs"
# seq_to_species_function <- function(s){
# 	if (substr(s, 1, 5) == "Vocar") {return("Volca_v2.0")}
# 	else if (substr(s, 1, 9) == "gnl|CMER|") {return("Cyame_v1.0")}
# 	else if (substr(s, 1, 3) == "Gb_") {return("Ginbi_v1.0")}
# 	else if (substr(s, 1, 3) == "TnS" || substr(s, 1, 3) == "TnC") {return("Gnemo_v1.0")}
# 	else if (substr(s, 1, 6) == "Mapoly") {return("Marpo_v1.0")}
# 	else {return(sub("\\|.*$", "", sub("^.*?\\|", "", s, perl=T), perl=T))}
# }
# # relative<-"sequence_counts_refs.txt"

# study 3: Subset of Pfam-A HMMs used in this study with sequences from all streptophytes removed on OneKP data set
case_id<-"Pfam-A_noStrepto_on_OneKP_prot_cdhit99"
seq_to_species_function <- function(s){return(substr(s, 1, 4))}
# relative<-"data/sequence_counts_prot_cdhit_99.txt"
completeness<-"data/CEGMA_BUSCO_100-missing.tsv"

################################
# END: uncomment desired study
################################

dir<-case_id
if (!is.null(relative)) {
	dir <- paste0(dir, "_rel")
}

# data
hmmer_raw<-as.data.table(read.table(paste0("results/hmmsearch/hmmsearch_" , case_id, ".txt"), fill=T, flush=T))
setnames(hmmer_raw, c(1, 5, 3), c("subject", "score", "query"))
setkey(hmmer_raw, query, subject, score)

samples<-sort(unique(sapply(hmmer_raw[,subject], seq_to_species_function)))

hmmer_dom <- split(hmmer_raw[,subject], hmmer_raw[,query])
hmmer_dom<-lapply(hmmer_dom, as.vector)

hmmer_families <- list(
	"F-box" = unique(unlist(hmmer_dom[c("F-box", "F-box-like", "F-box-like_2")])),
	"MADS-box" = hmmer_dom[["MADS"]],
	"GT1" = hmmer_dom[["GT1"]],
	"GH28" = hmmer_dom[["GH28"]],
	"ABI3_VP1" = setdiff(hmmer_dom[["B3"]], c(hmmer_dom[["AP2"]], hmmer_dom[["Auxin_resp"]], hmmer_dom[["WRKY"]])),
	"AP2_EREBP" = hmmer_dom[["AP2"]],
	"bHLH" = hmmer_dom[["HLH"]],
	"bZIP" = setdiff(unlist(hmmer_dom[c("bZIP_1", "bZIP_2")]), c(hmmer_dom[["HLH"]], hmmer_dom[["Homeobox"]])),
	"C2H2" = setdiff(hmmer_dom[["zf-C2H2"]], hmmer_dom[["zf-MIZ"]]),
	"C3H" = setdiff(hmmer_dom[["zf-CCCH"]], c(hmmer_dom[["AP2"]], hmmer_dom[["SRF-TF"]], hmmer_dom[["Myb_DNA-binding"]], hmmer_dom[["zf-C2H2"]])),
	"HB" = setdiff(hmmer_dom[["Homeobox"]], c(hmmer_dom[["EIN3"]], hmmer_dom[["KNOX1"]], hmmer_dom[["KNOX2"]], hmmer_dom[["HALZ"]], hmmer_dom[["bZIP_1"]])),
	"WRKY" = hmmer_dom[["WRKY"]],
	"ABC" = hmmer_dom[["ABC_tran"]],
	"PP2C" = hmmer_dom[["PP2C"]],
	"PPR P-class" = setdiff(hmmer_dom[["P"]], c(hmmer_dom[["L1"]], hmmer_dom[["E1"]], hmmer_dom[["S1"]])),
	"PPR PLS-class" = unique(unlist(hmmer_dom[c("L1", "E1", "S1")])),
	"Receptor-like Kinases" = unique(unlist(hmmer_dom[c("Pkinase", "Pkinase_Tyr")])),
	"P450s" = hmmer_dom[["p450"]],
	"MYB" = setdiff(hmmer_dom[["Myb_DNA-binding"]], c(hmmer_dom[["ARID"]], hmmer_dom[["Response_reg"]], hmmer_dom[["G2-like_Domain"]], hmmer_dom[["trihelix"]])),
	"NAC" = hmmer_dom[["NAM"]],
	"NBS-LRR" = hmmer_dom[["NB-ARC"]],
	"PHD" = setdiff(hmmer_dom[["PHD"]], c(hmmer_dom[["Myb_DNA-binding"]], hmmer_dom[["Alfin-like"]], hmmer_dom[["ARID"]], hmmer_dom[["DDT"]], hmmer_dom[["Homeobox"]], hmmer_dom[["JmjC"]], hmmer_dom[["JmjN"]], hmmer_dom[["SWIB"]], hmmer_dom[["zf-TAZ"]], hmmer_dom[["zf-MIZ"]], hmmer_dom[["zf-CCCH"]])),
	"GRAS" = setdiff(hmmer_dom[["GRAS"]], hmmer_dom[["zf-TAZ"]])
)

# remove empty families
if(length(which(lapply(hmmer_families, length) == 0)) > 0 ) {
	hmmer_families<-hmmer_families[-which(lapply(hmmer_families, length) == 0)]
}

sample_counts<-do.call(rbind, lapply(hmmer_families, function(fam){
	fam_tmp<-sapply(fam, seq_to_species_function)
	return(
		sapply(samples, function(sample){return(
			sum(sample == fam_tmp)
		)}
	)
)}))

meta_sample_counts <- matrix(0, nrow=3, ncol=dim(sample_counts)[2], dimnames=list(c("aa-sequences", "BUSCO", "all_gene_families"), dimnames(sample_counts)[[2]]))

if (!is.null(relative)) {
	sample_size<-read.table(relative, row.names=1, sep=":")
	dimnames(sample_size)[[1]]<-sapply(dimnames(sample_size)[[1]], substr, 10, 13)
	sample_counts <- t(t(sample_counts)/sample_size[dimnames(sample_counts)[[2]],1]) * 1e4#mean(sample_size[,1])
	meta_sample_counts["aa-sequences",] <- sample_size[dimnames(sample_counts)[[2]],1]
	rm(sample_size)
}
if (!is.null(completeness)) {
	completeness_data <- read.table(completeness, header=T)[,2]
	names(completeness_data) <- sapply(read.table(completeness, header=T)[,1], substr, 1, 4)
	meta_sample_counts["BUSCO",] <- completeness_data[dimnames(sample_counts)[[2]]]
}

weights<-1/apply(sample_counts, 1, max)
# weights<-weights/sum(weights)

meta_sample_counts["all_gene_families", ] <-  apply(sample_counts * weights, 2, mean)
meta_sample_counts["all_gene_families", ] <- meta_sample_counts["all_gene_families", ] * mean(apply(sample_counts, 1, max))

rm(weights)


# Phylogeny
## Clade Dictionary
tmp<-read.table("data/phylogeny/clade_dictionary_edit.tsv", row.names=1, sep="\t")
clade_dictionary<-as.character(tmp[,2])
names(clade_dictionary)<-dimnames(tmp)[[1]]
rm(tmp)

## Version 17-05
tree<-read.nexus("data/phylogeny/astral-33-rooted.nex.tre")
tree<-root(tree, outgroup=names(clade_dictionary[clade_dictionary=="Chromista"]), resolve.root=T)
tree$tip.label<-gsub("'", "", tree$tip.label)

## Species names
tmp<-read.table("data/phylogeny/species_dictionary.tsv", row.names=1, sep="\t")
species_names<-as.character(tmp[,1])
names(species_names)<-dimnames(tmp)[[1]]
rm(tmp)

## if reference study rename some samples
if (case_id == "mlgf_on_refs") {
	tree$tip.label[tree$tip.label == "QFND"] <- "Cyapa_v1.0"
	tree$tip.label[tree$tip.label == "JPYU"] <- "Marpo_v1.0"
	tree$tip.label[tree$tip.label == "SGTW"] <- "Ginbi_v1.0"
	tree$tip.label[tree$tip.label == "GTHK"] <- "Gnemo_v1.0"
}

## transform tree
tree<-ladderize(tree, right=F)
tree<-rotate(tree, 1405)
tree<-rotate(tree, 1415)
tree<-rotate(tree, 1416)
write.tree(tree, file="tmp_tree.newick")
tree<-read.tree("tmp_tree.newick")
tree<-rotate(tree, 2122)
tree<-rotate(tree, 2123)

write.tree(tree, file="tmp_tree.newick")
tree<-read.tree("tmp_tree.newick")

file.remove("tmp_tree.newick")

## update clade dictionary and clade list
ordered_tips <- tree$tip.label[tree$edge[tree$edge[,2] <= length(tree$tip.label), 2]]
clades<-unique(clade_dictionary[ordered_tips])

## write out results for later use
dir.create(paste0("results/", dir))
write.tree(tree, file=paste0("results/", dir, "/tree.newick"))

get_data<-function(counts, tree, indices=tree$tip.label) {
	i<-indices
	i<-intersect(i, tree$tip.label)
	i<-setdiff(i, name.check(tree, counts)[[1]])
	i<-setdiff(i, as.character(read.table("data/contamination.txt")[,1]))
	if (!is.null(completeness)) {
		i<-setdiff(i, names(which(completeness_data < 57.5)))
	}
	return(counts[i])
}
clade_wise_count<-function(data, fun){
	return(sapply(unique(clade_dictionary), function(clade){
		fun(get_data(data, tree, indices=names(which(clade_dictionary==clade))))
	}))
}
write.table(apply(sample_counts, 1, clade_wise_count, mean)[clades,], file=paste0("results/", dir, "/clade_mean_sizes.tsv"), quote=F, sep="\t")

# set colors: color scheme B
clade_col<-c("#654522", "#872c17", "#ddd300", "#b3446c", "#f6a601", "#604e97", "#f99379", "#0067a5", "#e68fac", "#018856", "#848482", "#c3b280", "#f38401", "#a1caf1", "#be0032", "#875692", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#909200")
names(clade_col)<-clades

# expansion test method
expansions <- matrix(NA, ncol=19, nrow=dim(sample_counts)[1])
all_expansions <- matrix(NA, ncol=19, nrow=dim(sample_counts)[1])
num_expansions <- matrix(NA, ncol=19, nrow=dim(sample_counts)[1])
pvalues <- matrix(NA, ncol=19, nrow=dim(sample_counts)[1])
dimnames(expansions) <- list(
	dimnames(sample_counts)[[1]],
	c("Chromista -> Rhodophyta", "Rhodophyta -> Glaucophyta, Chlorophyta", "Chlorophyta -> Streptophyte Algae", "Streptophyte Algae -> Hornworts, Liverworts, Mosses", "Streptophyte Algae -> Hornworts, Liverworts", "Streptophyte Algae, Hornworts -> Liverworts", "Hornworts, Liverworts -> Mosses", "Hornworts, Liverworts -> Lycophytes", "Hornworts, Liverworts, Mosses -> Lycophytes", "Lycophytes -> Monilophytes", "Monilophytes -> Gymnosperms", "Gymnosperms -> ANA grade, Monocots", "ANA grade, Monocots -> Chloranthales, Magnoliids", "Chloranthales, Magnoliids -> CRPT grade", "CRPT grade -> Asterids", "Asterids -> Caryophyllales", "Asterids, Caryophyllales -> Rosids, Saxifragales, Santalales", "Rosids, Saxifragales -> Santalales", "Rosids -> Saxifragales")
)
dimnames(all_expansions)<-dimnames(expansions)
dimnames(num_expansions)<-dimnames(expansions)
dimnames(pvalues)<-dimnames(expansions)

get_data_for_plotting<-function(counts, tree, fam, indices) {
	i<-indices
	i<-setdiff(i, name.check(tree, t(counts))[[1]])
	res<-rep(0, length(indices))
	names(res)<-indices
	res[i]<-counts[fam, i]
	
	return(res)
}

get_color_for_plotting<-function(clade_col, tree, fam, indices) {
	i<-indices
	i<-setdiff(i, name.check(tree, t(sample_counts))[[1]])
	col<-rep("white", length(indices))
	names(col)<-indices
	col[i]<-clade_col[clade_dictionary[i]]
	col[intersect(indices, as.character(read.table("data/contamination.txt")[,1]))]<-"white"
	if (!is.null(completeness)) {
		col[intersect(indices, names(which(completeness_data < 0.575)))]<-"white"
		col[intersect(indices, names(which(is.na(completeness_data))))]<-"white"
	}
	return(col)
}

test_clades <- function(a, b, i=NULL, family=fam, test=ks.test, edge=NULL){
	a_data <- get_data(sample_counts[family,], tree, a)
	b_data <- get_data(sample_counts[family,], tree, b)

	pval <- test(a_data, b_data)$p.value

	fc <- mean(b_data, trim=0.05) / mean(a_data, trim=0.05)

	print(paste0("p-value= ", pval, " , fold-change= ", fc, " from= ", mean(a_data), " to= ", mean(b_data)))
	
	if (!is.null(i)) {
		if(is.infinite(fc)) {
			all_expansions[family, i] <<- paste0("[", round(mean(b_data), digits=1), "]")
		}
		else {
			all_expansions[family, i] <<- as.character(round(fc, digits=1))
		}
		pvalues[family, i] <<- as.character(pval)
	}
	if(max(fc, 1/fc, na.rm=T) >= 1.5 && pval <= 1e-6) {
		if (!is.null(edge)) {
			edgelabels(text=format(fc, digits=1, nsmall=1), edge=edge, frame="c", cex=0.7)
		}
		if (!is.null(i)) {
			num_expansions[family, i] <<- round(fc, digits=1)
			if(is.infinite(fc)) {
				expansions[family, i] <<- paste0("[", round(mean(b_data), digits=1), "]")
			}
			else {
				expansions[family, i] <<- as.character(round(fc, digits=1))
			}
			pvalues[family, i] <<- as.character(pval)
		}
	}
}

cl_mem <- function(cl){
	ret <- c()
	for (i in 1:length(cl)) {
		ret <- c(ret, names(which(clade_dictionary == cl[i])))
	}
	return(unique(ret))
}

plot_tree_hist<-function(tree, data, main="", ylab="Gene family sizes", mark_exp="none", tips=NULL, legend="topright", bars="auto", cex=0.7) {#, tree.file=NULL) {
	
	# remove tips
	if(!is.null(tips)) {
		tree <- drop.tip(tree, setdiff(tree$tip.label, tips))
		data <- data[tree$tip.label]
		ordered_tips <- tree$tip.label[tree$edge[tree$edge[,2] <= length(tree$tip.label), 2]]
		clades<-unique(clade_dictionary[ordered_tips])
	}

	# precomputations
	edge.color<-apply(tree$edge, 1, function(x){
		if(length(table(clade_dictionary[tips(tree, x[2])])) == 1) {
			return(clade_col[clade_dictionary[tips(tree, x[2])[1]]])
		} else if( length(table(clade_dictionary[tips(tree, x[1])])) == 1 ) {
			return(clade_col[clade_dictionary[tips(tree, x[1])[1]]])
		} else {
			return("black")
		}
	})

	# plotting
	## upper half: phylogentic tree
	layout(matrix(c(1,2), 2, 1, byrow = TRUE), heights=c(1,1))
	old_par<-par(mar=c(0,5,2,0)+0.1, cex=cex)#, cex.axis=0.6

	if(bars == "single" | (bars == "auto" & length(tree$tip.label) <= 215) ) {
		tmp_tree<-tree
		tmp_tree$tip.label<-paste(species_names[tree$tip.label], tree$tip.label)
		plot(tmp_tree, edge.color=edge.color, show.tip.label=T, direction="downwards", no.margin=F, main=main, cex=0.3, align.tip.label=T)#, edge.width=edge.width
		rm(tmp_tree)
	} else {
		plot(tree, edge.color=edge.color, show.tip.label=F, direction="downwards", no.margin=F, y.lim=c(0,73), adj=0.5, main=main)#, edge.width=edge.width
	}
	if(!is.na(legend) && legend %in% c("bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right", "center")) {
		legend(legend, clades, fill=clade_col[clades], ncol=2, bty = "n")#, cex=0.7	
	}

	if(mark_exp == "curated" && is.null(tips)) {
		get.edge<-function(x, y){
			intersect(ceiling(which(t(tree$edge) == x)/2), ceiling(which(t(tree$edge) == y)/2))
		}

		test_clades(cl_mem(c("Chromista")), cl_mem("Rhodophyta"), 1, edge=66)
		test_clades(cl_mem("Rhodophyta"), cl_mem(c("Glaucophyta", "Chlorophyta")), 2, edge=124)
		test_clades(cl_mem("Chlorophyta"), cl_mem("Streptophyte Algae"), i=3, edge=372)
		test_clades(cl_mem("Streptophyte Algae"), cl_mem(c("Hornworts", "Liverworts", "Mosses")), 4, edge=470)
		test_clades(cl_mem("Streptophyte Algae"), cl_mem(c("Hornworts", "Liverworts")), 5)
		test_clades(cl_mem(c("Streptophyte Algae", "Hornworts")), cl_mem(c("Liverworts")), 6)
		test_clades(cl_mem(c("Hornworts", "Liverworts")), cl_mem(c("Mosses")), 7)
		test_clades(cl_mem(c("Hornworts", "Liverworts")), cl_mem(c("Lycophytes")), 8)
		test_clades(cl_mem(c("Hornworts", "Liverworts", "Mosses")), cl_mem(c("Lycophytes")), 9, edge=618)
		test_clades(cl_mem(c("Lycophytes")), cl_mem(c("Monilophytes")), 10, edge=662)
		test_clades(cl_mem(c("Monilophytes")), cl_mem(c("Gymnosperms")), 11, edge=808)
		test_clades(cl_mem(c("Gymnosperms")), cl_mem(c("ANA grade", "Monocots")), 12, edge=972)
		test_clades(cl_mem(c("ANA grade", "Monocots")), cl_mem(c("Chloranthales", "Magnoliids")), 13, edge=1202)
		test_clades(cl_mem(c("Chloranthales", "Magnoliids")), cl_mem(c("CRPT grade")), 14, edge=1256)
		test_clades(cl_mem(c("CRPT grade")), cl_mem(c("Asterids")), 15, edge=1328)
		test_clades(cl_mem(c("Asterids")), cl_mem(c("Caryophyllales")), 16, edge=1761)
		test_clades(cl_mem(c("Asterids", "Caryophyllales")), cl_mem(c("Rosids", "Saxifragales", "Santalales")), 17, edge=1884)
		test_clades(cl_mem(c("Rosids", "Saxifragales")), cl_mem(c("Santalales")), 18, edge=2344)
		test_clades(cl_mem(c("Rosids")), cl_mem(c("Saxifragales")), 19, edge=229/9)

	} else if(mark_exp == "auto") {
		for (edge_i in 1:dim(tree$edge)[1]) {
			edge <- tree$edge[edge_i,]
			if (length(tips(tree, edge[1])) > length(tips(tree, edge[2]))) {
				parent <- edge[1]
				child <- edge[2]
			}else{
				parent <- edge[2]
				child <- edge[1]
			}
			a_data <- get_data(sample_counts[fam,], tree, tips(tree, child))
			b_data <- get_data(sample_counts[fam,], tree, setdiff(tips(tree, parent), tips(tree, child)))
			if(length(a_data) > 2 && length(b_data) > 2) {
				fc <- mean(a_data, trim=0.05) / mean(b_data, trim=0.05)

				pval <- ks.test(a_data, b_data)$p.value
				
				if(fc >= 1.5 && pval <= 1e-6) {
					print(paste0(parent, "-", child, " ", fc))
					edgelabels(text=format(fc, digits=1, nsmall=1), edge=edge_i, frame="c", cex=cex)
				}
			}
			
		}
	}

	## bottom half: barplot
	par(mar=c(c(0.2, 5, 0.1 , 0))+0.1, mgp=c(4, 1, 0), cex=cex, las=2)#, cex.axis=0.6

	if(bars == "single" | (bars == "auto" & length(tree$tip.label) <= 215) ) {
		bar_col<-get_color_for_plotting(clade_col, tree, fam, tree$tip.label)
		bar_count<-get_data_for_plotting(sample_counts, tree, fam, tree$tip.label)
		barplot2(bar_count, col=bar_col, ylab=ylab, names.arg=NA)
	} else {
		cladewise_counts <- clade_wise_count(data, mean)
		se<-function(x, na.rm=T) sqrt(var(x, na.rm=na.rm)/length(x[!is.na(x)]))
		cladewise_se <- clade_wise_count(data, se)
		cladewise_sd <- clade_wise_count(data, sd)
		width<-table(clade_dictionary[tree$tip.label])[clades]
		barplot2(cladewise_counts[clades], plot.ci=T, ci.l=cladewise_counts[clades]-cladewise_se[clades], ci.u=cladewise_counts[clades]+cladewise_se[clades], col=clade_col[clades], width=width, space=0.0, ylab=ylab, axisnames=F)

		if(is.na(legend)) {
			midpoints<-rep(0, length(width))
			names(midpoints)<-clades
			midpoints[1]<-width[1]/2
			for (i in 2:length(width)) {
				midpoints[i]<-sum(width[1:i-1])+width[i]/2
			}
			for(clade in clades) {
				if(cladewise_counts[clade] < 0.25*max(cladewise_counts)){
					text(y=0.3*max(cladewise_counts), x=midpoints[clade], label=clade, adj=0, srt=90)
				} else if(cladewise_counts[clade] < 0.6*max(cladewise_counts)){
					text(y=0.65*max(cladewise_counts), x=midpoints[clade], label=clade, adj=0, srt=90)
				} else {
					text(y=1/40*max(cladewise_counts), x=midpoints[clade], label=clade, adj=0, srt=90)
				}
			}
		}
	}
}

## save total count table
tmp<-t(sapply(dimnames(sample_counts)[[1]], function(fam){get_data(sample_counts[fam,], tree)}))
tmp<-rbind(Species=species_names[colnames(tmp)], Clade=clade_dictionary[colnames(tmp)], tmp)

write.csv(tmp, file=paste0("results/", dir, "/sample_count_table.csv"))
write.table(tmp, file=paste0("results/", dir, "/sample_count_table.tsv"), sep="\t")

## plotting
dir.create(paste0("plots/", dir))
pdf(paste0("plots/", dir, "/full dendrogram.pdf"), width=200, height=28)
	edge.color<-apply(tree$edge, 1, function(x){
		if(length(table(clade_dictionary[tips(tree, x[2])])) == 1) {
			return(clade_col[clade_dictionary[tips(tree, x[2])[1]]])
		}
		else if( length(table(clade_dictionary[tips(tree, x[1])])) == 1 ) {
			return(clade_col[clade_dictionary[tips(tree, x[1])[1]]])
		}
		else {
			return("black")
		}
	})
	par(mar=rep(0.1, 4))
	tmp_tree<-tree
	tmp_tree$edge.length<-tmp_tree$edge.length+1
	plot(tmp_tree, direction="downwards", no.margin=T, edge.color=edge.color, edge.width=6)
	nodelabels()
	rm(tmp_tree)
dev.off()

if (!is.null(relative)) {
	fam<-"aa-sequences"
	pdf(paste0("plots/", dir, "/", fam, ".pdf"), 183/25.4, 150/25.4)
	plot_tree_hist(tree, meta_sample_counts[fam,], main="", ylab="Number of inferred protein sequences", legend=NA)
	dev.off()
}

if (!is.null(completeness)) {
	fam<-"BUSCO"
	pdf(paste0("plots/", dir, "/", fam, ".pdf"), 183/25.4, 150/25.4)
	plot_tree_hist(tree, meta_sample_counts[fam,], main="", ylab="Completeness [%]", legend=NA)
	dev.off()
}

dir.create(paste0("plots/", dir, "/KS-test_sigFC/"))
fam<-"all_gene_families"
pdf(paste0("plots/", dir, "/KS-test/", fam, ".pdf"), 183/25.4, 150/25.4)
plot_tree_hist(tree, meta_sample_counts[fam,], main=paste0("Weighted Mean of Size from ", dim(sample_counts)[1] , " Gene Families"), ylab="", legend=NA, bars="folded", mark_exp="none")
dev.off()

for (fam in dimnames(sample_counts)[[1]]) {
	print(fam)
	pdf(paste0("plots/", dir, "/KS-test_sigFC/", fam, ".pdf"), 183/25.4, 150/25.4)
	plot_tree_hist(tree, sample_counts[fam,], main=paste(fam, "genes, fold-changes (at least 50% gain or 33% loss and p < 1e-6)"), mark_exp="curated", bars="folded", legend=NA)
	dev.off()
}

dir.create(paste0("results/", dir, "/KS-test/"))
write.table(expansions, file=paste0("results/", dir, "/KS-test/sig_expansions.csv"), sep=";", dec=",", na="", col.names=NA)
write.table(all_expansions, file=paste0("results/", dir, "/KS-test/all_expansions.csv"), sep=";", dec=",", na="", col.names=NA)
write.table(pvalues, file=paste0("results/", dir, "/KS-test/all_pvalues.csv"), sep=";", dec=",", na="", col.names=NA)

dir.create(paste0("plots/", dir, "/no_FC/"))
for (fam in dimnames(sample_counts)[[1]]) {
	print(fam)
	pdf(paste0("plots/", dir, "/no_FC/", fam, ".pdf"), 183/25.4, 150/25.4)
	plot_tree_hist(tree, sample_counts[fam,], main=paste(fam, "genes, fold-changes (at least 50% gain or 33% loss and p < 1e-6)"), mark_exp="none", bars="folded", legend=NA)
	dev.off()
}

dir.create(paste0("plots/", dir, "/allFC/"))
for (fam in dimnames(sample_counts)[[1]]) {
	print(fam)
	pdf(paste0("plots/", dir, "/allFC/", fam, ".pdf"), 183/25.4, 150/25.4)
	
	plot_tree_hist(tree, sample_counts[fam,], main=paste(fam, "proteins, fold-changes"), mark_exp="auto")

	dev.off()
}

dir.create(paste0("plots/", dir, "/allFC_ga-mosses/"))
for (fam in dimnames(sample_counts)[[1]]) {
	print(fam)
	pdf(paste0("plots/", dir, "/allFC_ga-mosses/", fam, ".pdf"), 183/25.4, 130/25.4)
	
	plot_tree_hist(tree, sample_counts[fam,], tips=names(clade_dictionary[clade_dictionary %in% c("Chlorophyta", "Streptophyte Algae", "Streptophyte Algae-Desmidiales", "Hornworts", "Liverworts", "Mosses")]), mark_exp="auto", main=paste(fam, "Genes, fold-changes (at least 50% gain or 33% loss and p < 1e-6)"), legend="topright" )

	dev.off()
}

dir.create(paste0("plots/", dir, "/allFC_ga-bryophytes/"))
for (fam in dimnames(sample_counts)[[1]]) {
	print(fam)
	pdf(paste0("plots/", dir, "/allFC_ga-bryophytes/", fam, ".pdf"), 183/25.4, 130/25.4)
	
	plot_tree_hist(tree, sample_counts[fam,], tips=names(clade_dictionary[clade_dictionary %in% c("Chlorophyta", "Streptophyte Algae", "Streptophyte Algae-Desmidiales", "Hornworts", "Liverworts")]), mark_exp="auto", main=paste(fam, "Genes, fold-changes (at least 50% gain or 33% loss and p < 1e-6)"), legend="topright" )

	dev.off()
}

dir.create(paste0("plots/", dir, "/allFC_gymnosperms/"))
for (fam in dimnames(sample_counts)[[1]]) {
	print(fam)
	pdf(paste0("plots/", dir, "/allFC_gymnosperms/", fam, ".pdf"), 183/25.4, 130/25.4)
	
	plot_tree_hist(tree, sample_counts[fam,], tips=names(clade_dictionary[clade_dictionary %in% c("Gymnosperms")]), mark_exp="auto", main=paste(fam, "Genes, fold-changes (at least 50% gain or 33% loss and p < 1e-6)"), legend="topright" )

	dev.off()
}

dir.create(paste0("plots/", dir, "/allFC_monilophytes/"))
for (fam in dimnames(sample_counts)[[1]]) {
	print(fam)
	pdf(paste0("plots/", dir, "/allFC_monilophytes/", fam, ".pdf"), 183/25.4, 130/25.4)
	
	plot_tree_hist(tree, sample_counts[fam,], tips=names(clade_dictionary[clade_dictionary %in% c("Monilophytes")]), mark_exp="auto", main=paste(fam, "Genes, fold-changes (at least 50% gain or 33% loss and p < 1e-6)"), legend="topright" )

	dev.off()
}

dir.create(paste0("plots/", dir, "/allFC_Chloranthales_to_Asterids/"))
for (fam in dimnames(sample_counts)[[1]]) {
	print(fam)
	pdf(paste0("plots/", dir, "/allFC_Chloranthales_to_Asterids/", fam, ".pdf"), 183/25.4, 130/25.4)

	plot_tree_hist(tree, sample_counts[fam,], tips=names(clade_dictionary[clade_dictionary %in% c("Chloranthales", "Magnoliids", "CRPT grade", "Asterids")]), mark_exp="auto", main=paste(fam, "Genes, fold-changes (at least 50% gain or 33% loss and p < 1e-6)"), legend="topright", bars="single" )

	dev.off()
}

dir.create(paste0("plots/", dir, "/allFC_Chromista_to_Chlorophyta/"))
for (fam in dimnames(sample_counts)[[1]]) {
	print(fam)
	pdf(paste0("plots/", dir, "/allFC_Chromista_to_Chlorophyta/", fam, ".pdf"), 183/25.4, 130/25.4)

	plot_tree_hist(tree, sample_counts[fam,], tips=names(clade_dictionary[clade_dictionary %in% c("Chromista", "Rhodophyta", "Glaucophyta", "Chlorophyta")]), mark_exp="auto", main=paste(fam, "Genes, fold-changes (at least 50% gain or 33% loss and p < 1e-6)"), legend="topright", bars="single" )

	dev.off()
}
