#!/bin/R

rm(list=ls())
closeAllConnections()

library(openxlsx)

dir.create('./metabolon_metabolomics/interval/gwas/conditional_analysis/snptest_out_all_locus', showWarnings=F)
setwd("./metabolon_metabolomics/interval/gwas/conditional_analysis/snptest_out_all_locus")

start.time <- Sys.time()
system('clear')
cat("\nStarted on ", format(Sys.time(), "%a %b %d %X %Y"))

##Load all required programs
snptest.exe <- "/usr/local/Cluster-Apps/ceuadmin/snptest_new/2.5.2/snptest"
plink.exe <- "/usr/local/Cluster-Apps/ceuadmin/plink_linux_x86_64_beta3.32/plink"
qctool.exe <- "/usr/local/Cluster-Apps/ceuadmin/qctool_v1.4-linux-x86_64/qctool"
interval.dir <- "./metabolon_metabolomics/interval/gwas/conditional_analysis/interval_locus_bgen_files/"
epicn.dir <- "./metabolon_metabolomics/interval/gwas/conditional_analysis/epic_norfolk_data"

	args <- commandArgs(trailing=T)
	locus.no <- args[1]
	load('./metabolon_metabolomics/interval/gwas/4_locus_definition/2_distance_based_definition/Conditional_Analysis_Loci_Reference.RData')

	dir.out <- condition.ref[locus==locus.no, paste0(chromosome, "_", coordinates_conditional_analysis)]
	unlink(dir.out, recursive=T)
	dir.create(dir.out, showWarnings=F)
	setwd(dir.out)

	interval.files <- list.files(paste0(interval.dir,gsub("-","_",dir.out)), full.names=T)
	interval.gen.fname <- interval.files[grep('.gen$',interval.files)]
	interval.samples.fname <- interval.files[grep('.samples$',interval.files)]

	epicn.gen.fname <- list.files(epicn.dir, full.names=T)[grep(paste0('chr',gsub("-","_",dir.out),".gen.gz"), list.files(epicn.dir, full.names=T))]
	epicn.samples.fname <- list.files(epicn.dir, full.names=T)[grep(paste0('chr',gsub("-","_",dir.out),".sample$"), list.files(epicn.dir, full.names=T))]


		cat("\n\nReading INTERVAL gen and sample files...\n")
		cat(paste0("\n\t", interval.gen.fname, "\n\t", interval.samples.fname, "\n"))
			interval.gen <- fread(interval.gen.fname,header=F,colClasses="character",showProgress=F)
		cat("\nReading EPIC-Norfolk gen and sample files...\n")
		cat(paste0("\n\t", epicn.gen.fname, "\n\t", epicn.samples.fname, "\n"))
			epicn.gen <- fread(cmd=paste0('zcat ',epicn.gen.fname),header=F,colClasses="character",showProgress=F)
			epicn.gen[,c('V2','V3'):=list(V3,V2)]

		cat("\n\tChecking files for consistency...\n")

			var.common <- intersect(interval.gen$V2, epicn.gen$V2)
			
			cat("\n\t\tNumber of variants in INTERVAL:", nrow(interval.gen),"\n")
			cat("\t\tNumber of variants in EPIC-Norfolk:", nrow(epicn.gen),"\n")
			cat("\t\tNumber of variants in INTERVAL but not in EPIC-Norfolk:", length(setdiff(interval.gen$V2, epicn.gen$V2)),"\n")
			cat("\t\tNumber of variants in EPIC-Norfolk but not in INTERVAL:", length(setdiff(epicn.gen$V2, interval.gen$V2)),"\n")
			cat("\t\tNumber of variants common between the two datasets:", length(var.common),"\n")
			if(nrow(interval.gen)==length(var.common) & nrow(epicn.gen)==length(var.common)) cat("\t\tNumber of variant mismatch between INTERVAL and EPIC-Norfolk: None\n")

			interval.gen.out <- interval.gen[match(var.common, V2),]
			epicn.gen.out <- epicn.gen[match(var.common, V2),]


				##Create dummy rsIDs for consistency and also sort the variant information as required by qctools!!!!!!!!!!!!!!!!!
				interval.gen.out[,V3:=paste0("rs",1:nrow(interval.gen.out))]
				interval.gen.out <- interval.gen.out[order(V1,V2,V3,V4),]
				epicn.gen.out <- epicn.gen.out[match(interval.gen.out$V2, V2),]
				epicn.gen.out[,V3:=interval.gen.out$V3]	
				
				cat("\n\tWriting formatted INTERVAL and EPIC-Norfolk files for merging...\n")

				write.table(interval.gen.out, 'interval.gen', col.names=F, row.names=F, quote=F, sep=" ")
				write.table(epicn.gen.out, 'epicn.gen', col.names=F, row.names=F, quote=F, sep=" ")

				##Fix the issue with ID_1 column being all 0s in EPIC-Norfolk sample files
				epic.info <- fread(epicn.samples.fname, header=T, colClasses="character", showProgress=F)
				epic.info[,ID_1:=ID_2]

				write.table(epic.info, 'epicn.samples', col.names=T, row.names=F, quote=F, sep=" ")

					##Combine data using qctool
					cat("\tMerging formatted data...\n")
					sys.cmd <- paste0(qctool.exe, " -g interval.gen -s ",interval.samples.fname, " -g epicn.gen -s epicn.samples -match-alleles-to-cohort1 -og interval_epicn.bgen -os interval_epicn.samples 2> /dev/null")
                                        system(sys.cmd)

					##Create snp summary to check the total number of variants in the final merged dataset
					cat("\tCalculating snp summaries in the final merged dataset...\n")
					sys.cmd <- paste0(qctool.exe, " -g interval_epicn.bgen -s interval_epicn.samples -snp-stats interval_epicn_snp_summary.stats 2> /dev/null")
					system(sys.cmd)

					snp.stats.combined <- fread('interval_epicn_snp_summary.stats', header=T, colClasses="character")
					cat("\n\t\tNumber of variants in interval_epicn.bgen:", nrow(snp.stats.combined),"\n")
                                	##Create VCF files to later calculate LD for locus plots
					#cat("\tCreating vcf file to calcualte LD...\n")
                                	#sys.cmd <- paste0(qctool.exe," -g interval_epicn.bgen -og interval_epicn.vcf 2> /dev/null")
                                	#system(sys.cmd)	
	
		cat("\nCreating phenotype file for SNPTEST analysis...\n")
			
			metabolites.ref <- unlist(strsplit(condition.ref[locus==locus.no, metabolites_condition_final], split=","))

			interval.pheno <- fread('./metabolon_metabolomics/interval/gwas/conditional_analysis/reference_files/INTERVAL_Metabolon_selectedForConditionalAnalysis_PS_280717.pheno',header=T,colClasses="character",showProgress=F)
			interval.pheno[,study:="1"]
			epicn.pheno <- fread('./metabolon_metabolomics/interval/gwas/conditional_analysis/epic_norfolk_data/epicn_metabolon_Apr2016_gwascohort_residuals_5841_767metabolites_9renamedCOMPIDs.txt',header=T,colClasses="character",showProgress=F)
			epicn.pheno[,study:="2"]

				epicn.pheno.sub <- epicn.pheno[,c("FID", "IID", "study", paste0("EPICN_OMICS_PC", 1:5), paste0("resw_", metabolites.ref)), with=F]
				setnames(epicn.pheno.sub, c('FID','IID', "study", paste0('PC', 1:5), metabolites.ref))
				##Scale PCs by multiplying the PCs x 100 as Isobel suggested
				epicn.pheno.sub[, c(match(paste0('PC', 1:5), colnames(epicn.pheno.sub))) := lapply(.SD, function(x) as.character(as.numeric(x)*100)), .SDcols=c(match(paste0('PC', 1:5), colnames(epicn.pheno.sub)))]

				interval.pheno.sub <- interval.pheno[-1, c("id_1", "id_2", "study", paste0("PC_", 1:5), metabolites.ref), with=F]
				setnames(interval.pheno.sub, c('FID','IID', "study", paste0('PC', 1:5), metabolites.ref))
		
				combined.pheno <- rbind(interval.pheno.sub, epicn.pheno.sub)
				combined.pheno[,FID:=IID]

				pc.tmp <- NULL
				pc.tmp <- cbind(combined.pheno[,FID],setnames(combined.pheno[,grep('PC',colnames(combined.pheno)),with=F], paste0('EPICN_PC_',1:5)),
							setnames(combined.pheno[,grep('PC',colnames(combined.pheno)),with=F], paste0('INTERVAL_PC_',1:5)))
				pc.tmp[match(epicn.pheno.sub$IID, V1), c(7:11) := "0"]
				pc.tmp[match(interval.pheno.sub$IID, V1), c(2:6) := "0"]


					pheno.tmp <- cbind(combined.pheno[,.(FID, IID, study)], pc.tmp, combined.pheno[,metabolites.ref,with=F])
					pheno.tmp[,V1:=NULL]

					combined.pheno <- NULL
					combined.pheno <- pheno.tmp

			pheno.ref <- fread('interval_epicn.samples',header=T,colClasses="character",showProgress=F)
			pheno.ref[,ID_1:=ID_2]
			pheno.ref <- pheno.ref[,.(ID_1,ID_2,missing)]
			pheno.ref <- pheno.ref[-1,]

			#ID_1 ID_2 missing study M01417
			#0 0 0 D P

			pheno.ref <- cbind(pheno.ref, combined.pheno[match(pheno.ref[,ID_1], FID),])
			pheno.ref[,c('FID','IID'):=NULL]
				string.dat <- data.table(matrix(c("0","0","0","D", rep("C", 10), rep("P", ncol(pheno.ref)-14)), nrow=1))
				setnames(string.dat, colnames(pheno.ref))
			pheno.ref <- rbind(string.dat, pheno.ref)
			pheno.ref[ pheno.ref == "-9" ] <- NA
			write.table(pheno.ref, 'SNPTEST.pheno', col.names=T, row.names=F, quote=F, sep=" ")

			###Perform conditional analysis of each metabolite using SNPTEST

				cat("\nPerforming conditional analysis for each metabolite...")
				metabolite.run <- colnames(pheno.ref)[15:ncol(pheno.ref)]
				n.metabolite <- length(metabolites.ref)
				vars.qcd <- NULL
				for(metabolite in metabolite.run) vars.qcd <- c(vars.qcd, fread(paste0('./metabolon_metabolomics/interval/gwas/conditional_analysis/metabolites_variants_list/', metabolite, '/qcd.snps'), header=F,sep="\n")[[1]])
				vars.qcd <- unique(vars.qcd)
				n.variants <- length(intersect(vars.qcd, var.common))
				cat("\nNumber of variants used for correction:", n.variants)

				###
				#The pvalue below is based on the pruning done as described in the script - /home/ps629/scripts/metabolomics/3_0_conditional_analysis_v1_exact_SNPTEST_getPValueForReference.R
				logp.conditional <- -log10(1.24741348813236e-08)
				#logp.conditional <- -log10((0.05/(n.metabolite*n.variants)))
				cat("\n-log(Pvalue) threshold used for conditional analysis:", logp.conditional)

				index <- 0

				for(metabolite in colnames(pheno.ref)[15:ncol(pheno.ref)]) {
				index <- index+1

					assoc.tmp <- fread(paste0('./metabolon_metabolomics/interval/gwas/2_metal_significant_results/',metabolite,'/',metabolite,'_primary_associations.txt'),header=T,colClasses="character")

					chromosome.ref <- as.numeric(condition.ref[locus==locus.no, chromosome])
					start.locus <- as.numeric(unlist(strsplit(condition.ref[locus==locus.no, coordinates_conditional_analysis], split="-")))[1]	
					end.locus <- as.numeric(unlist(strsplit(condition.ref[locus==locus.no, coordinates_conditional_analysis], split="-")))[2]

					assoc.tmp[,c('chromosome','position'):=tstrsplit(MarkerName, split=':', fixed=T)[1:2]]
					assoc.tmp[,chromosome:=as.numeric(gsub('chr','',chromosome))]
					assoc.tmp[,match(c('Pvalue_log10','chromosome','position'), colnames(assoc.tmp)) := lapply(.SD, as.numeric), .SDcols=match(c('Pvalue_log10','chromosome','position'), colnames(assoc.tmp))]

					assoc.tmp <- assoc.tmp[chromosome==chromosome.ref & position>start.locus & position<end.locus,]
					variants.qcd <- fread(paste0('./metabolon_metabolomics/interval/gwas/conditional_analysis/metabolites_variants_list/', metabolite, '/qcd.snps'), header=F, sep="\n")[[1]]
					assoc.tmp <- assoc.tmp[MarkerName%in%variants.qcd==T,]

					if(nrow(assoc.tmp)==0) stop("\n\nCheck what's going on with qcd list and the most significant variant...\n\n")
					assoc.tmp <- assoc.tmp[order(Pvalue_log10, decreasing=T),]

					variant.independent <- assoc.tmp[Pvalue_log10==max(Pvalue_log10), MarkerName]
					condition.list <- assoc.tmp[Pvalue_log10==max(Pvalue_log10), paste(chromosome, position, sep=":")]

					#variant.independent <- condition.ref[locus==locus.no, variant_condition]
					#condition.list <- paste(gsub("chr","",unlist(lapply(strsplit(condition.ref[locus==locus.no, variant_condition], split=":"), head, 2))), collapse=":")
					unlink(metabolite, recursive=T)
					dir.create(metabolite, showWarnings=F)
					setwd(metabolite)

					cat("\n\n\tProcessing metabolite ",index, "/", n.metabolite, "\n\n")
					chr.no <- as.numeric(unlist(strsplit(unlist(lapply(strsplit(getwd(), split="/"), tail, 2))[1], split="_"))[1])
					variants.qcd <- fread(paste0('./metabolon_metabolomics/interval/gwas/conditional_analysis/metabolites_variants_list/', metabolite, '/qcd.snps'), header=F)
					variants.qcd[,V1:=as.numeric(gsub('chr','',V1))]
					variants.qcd[,V2:=as.numeric(V2)]
					variants.qcd <- variants.qcd[order(V1, V2),]

					variants.qcd <- variants.qcd[V1==chr.no,]
					variants.qcd[,pos:=paste0("0",V1,":",V2)]
					if(chr.no>9) variants.qcd[,pos:=paste0(V1,":",V2)]
					write(variants.qcd$pos, 'metabolite_variants.txt')

					sys.cmd <- paste0(qctool.exe,
							" -g ../interval_epicn.bgen -incl-positions metabolite_variants.txt -og interval_epicn_sub.bgen 2> /dev/null")
					system(sys.cmd)
					
					iteration.index <- 0
					while(1) {

						cat("\r",paste0("\t\tIteration ", iteration.index))

						if(iteration.index==0) {
								snptest.cmd <- paste0(snptest.exe,
									" -data interval_epicn_sub.bgen ../SNPTEST.pheno",
									" -o SNPTEST_iteration_",iteration.index,".out",
									" -method expected",
									" -frequentist 1",
									" -cov_names study EPICN_PC_1 EPICN_PC_2 EPICN_PC_3 EPICN_PC_4 EPICN_PC_5 INTERVAL_PC_1 INTERVAL_PC_2 INTERVAL_PC_3 INTERVAL_PC_4 INTERVAL_PC_5",
									#" -cov_all",
									" -pheno ",metabolite,
									" -use_raw_phenotypes",
									" -chunk 1000 > SNPTEST_iteration_",iteration.index,".log")
								system(snptest.cmd)
								iteration.index <- iteration.index+1
						} else {
								 snptest.cmd <- paste0(snptest.exe,
									" -data interval_epicn_sub.bgen ../SNPTEST.pheno",
									" -o SNPTEST_iteration_",iteration.index,".out",
									" -method expected",
									" -frequentist 1",
									" -cov_names study EPICN_PC_1 EPICN_PC_2 EPICN_PC_3 EPICN_PC_4 EPICN_PC_5 INTERVAL_PC_1 INTERVAL_PC_2 INTERVAL_PC_3 INTERVAL_PC_4 INTERVAL_PC_5",
									#" -cov_all",
									" -pheno ",metabolite,
									" -condition_on ",paste(paste0("position=",condition.list), collapse=" "),
									" -use_raw_phenotypes",
									" -chunk 1000 > SNPTEST_iteration_",iteration.index,".log")

								system(snptest.cmd)
						
							assoc <- fread(paste0("grep chr SNPTEST_iteration_",iteration.index,".out"),header=T,showProgress=F,colClasses="character")
							col.convert <- match(c("chromosome","position","info","all_maf","all_total","frequentist_add_pvalue","frequentist_add_beta_1","frequentist_add_se_1"),colnames(assoc))
							assoc[, c(col.convert) := lapply(.SD, as.numeric), .SDcols=col.convert]
							assoc[, log10p:=-pchisq((frequentist_add_beta_1/frequentist_add_se_1)^2, df=1, lower.tail=F, log.p=T)/log(10)]
							assoc[,mac:=all_maf*(all_total*2)]
							assoc <- assoc[mac>10 & info>0.3,]
							assoc <- assoc[order(log10p, decreasing=T),]
						
							if(nrow(assoc[log10p>logp.conditional,])==0) break;
							variant.independent <- append(variant.independent, assoc[1,alternate_ids])
							condition.list <- append(condition.list, assoc[1,paste(chromosome, position, sep=":")])
							iteration.index <- iteration.index+1

						}

					}
					
					write(variant.independent, 'all_forward_independent.snps')
					setwd('..')

				}

###Use SNPTEST to calculate dosage for each variant
cat("\n\nCalculating dosages for conditionally independent variants...\n\n")
snps.conditional.list <- NULL
fname.list <- list.files(recursive=T, full.names=T)[grep(".snps$", list.files(recursive=T, full.names=T))]
for(fname in fname.list) {
	snps.conditional.list <- rbind(snps.conditional.list, fread(fname, header=F, sep=":", colClasses="character"))
}
snps.conditional.list <- unique(snps.conditional.list)
snps.conditional.list <- unique(snps.conditional.list[,V1:=gsub("chr","",V1)])

dosage.ref <- fread('SNPTEST.pheno',header=T,colClasses="character")
dosage.ref <- dosage.ref[-1,]

for(index in 1:nrow(snps.conditional.list)) {

	cat("\r", paste0("\tProcessing variant ",index, "/", nrow(snps.conditional.list)))
	condition.ref <- snps.conditional.list[index, paste0("pos~", paste(V1, V2, sep=":"))]
	qctool.cmd <- paste0(qctool.exe, " -g interval_epicn.bgen -s interval_epicn.samples -os temp.samples -condition-on ", condition.ref, " 2> /dev/null")
	system(qctool.cmd)

	dosage.temp <- fread('temp.samples',header=T,colClasses="character")
	dosage.ref <- cbind(dosage.ref, dosage.temp[match(dosage.ref[,ID_1], ID_1), ncol(dosage.temp), with=F])
	colnames(dosage.ref)[ncol(dosage.ref)] <- snps.conditional.list[index, paste(paste0("chr",V1), V2, V3, V4, sep="_")]
	unlink('temp.samples')

}

dosage.ref[,c(5:ncol(dosage.ref)):=lapply(.SD, as.numeric), .SDcols=c(5:ncol(dosage.ref))]
save(dosage.ref, file="./Dosages_REF.RData")

##Perform backward regression for each metabolite

cat("\n\nPerforming regression analysis for each metabolite with all metabolite specific conditionally independent variants in the model...\n\n")
index <- 0
var.conditional.ref <- NULL

                #tmp.fname <- list.files('./metabolon_metabolomics/interval/gwas/conditional_analysis/snptest_out_all_locus', full.names=T)[grep(paste0('slurm.*_', locus.no, '.out'), list.files('./metabolon_metabolomics/interval/gwas/conditional_analysis/snptest_out_all_locus', full.names=T))]
                #tmp.dat <- fread(paste0("egrep 'threshold used for conditional analysis' ", tmp.fname), header=F, sep=":")
                #p.threshold <- as.numeric(tmp.dat$V2[1])


unlink('metabolites_gws_forced_in.txt')

for(metabolite.ref in colnames(pheno.ref)[15:ncol(pheno.ref)]) {
index <- index+1

	setwd(metabolite.ref)
	cat("\r",paste0("\tProcessing metabolite ",index, "/", n.metabolite))
	var.condition <- gsub(':','_',fread('all_forward_independent.snps', header=F, sep="\n")[[1]])

	formula.string <- as.formula(paste0(metabolite.ref, " ~ ", paste(c(var.condition, paste0("EPICN_PC_", 1:5), paste0("INTERVAL_PC_",1:5), "study"), collapse=" + ")))
	lm.dat <- lm(formula.string, data=dosage.ref)

	summary.dat <- cbind(data.table(markername=names(summary(lm.dat)$coefficients[,4])), data.table(summary(lm.dat)$coefficients))
	summary.dat <- summary.dat[-1,]
	summary.dat <- summary.dat[-grep("EPICN_PC|INTERVAL_PC|study",markername),]
	setnames(summary.dat, c('markername', 'beta', 'se', 'tvalue', 'p'))
	summary.dat[,log10p:=-pchisq((beta/se)^2, df=1, lower.tail=F, log.p=T)/log(10)]
	summary.dat[, log10p:=as.numeric(round(log10p, digits=2), scientific=F)]

	if(nrow(summary.dat)==1 & nrow(summary.dat[log10p>logp.conditional,])==0) {

		var.conditional.ref <- rbind(var.conditional.ref, data.table(metabolite=metabolite.ref, variants=paste(summary.dat[,markername], collapse=",")))
	
	} else {

		var.conditional.ref <- rbind(var.conditional.ref, data.table(metabolite=metabolite.ref, variants=paste(summary.dat[log10p>logp.conditional,markername], collapse=",")))

	}

	setwd('..')

	write(paste(locus.no, metabolite.ref, summary.dat[,paste(markername,log10p, sep="\t")], sep="\t"), 'metabolites_gws_forced_in.txt', append=T)

}

cat("\n\nCreating final variant x metabolite association matrix...\n\n")
##Create final variant x metabolite association matrix
matrix.final <- data.table(markername=unique(unlist(strsplit(paste(unlist(var.conditional.ref[,.(variants)]),collapse=","),split=","))))
index <- 0
for(metabolite.ref in colnames(pheno.ref)[15:ncol(pheno.ref)]) {
index <- index+1

        cat("\r",paste0("\tProcessing metabolite ",index, "/", n.metabolite))

        ###Couldn't add study.. add study here...
        formula.string <- as.formula(paste0(metabolite.ref, " ~ ", paste(c(unlist(strsplit(var.conditional.ref[metabolite==metabolite.ref,variants], split=",")), paste0("EPICN_PC_", 1:5), paste0("INTERVAL_PC_",1:5), "study"), collapse=" + ")))
        lm.dat <- lm(formula.string, data=dosage.ref)
	summary.dat <- cbind(data.table(markername=names(summary(lm.dat)$coefficients[,4])), data.table(summary(lm.dat)$coefficients))
	summary.dat <- summary.dat[-1,]
	summary.dat <- summary.dat[-grep("EPICN_PC|INTERVAL_PC|study",markername),]

	setnames(summary.dat, c('markername', 'beta', 'se', 'tvalue', 'p'))

	summary.dat <- summary.dat[,log10p:=-pchisq((beta/se)^2, df=1, lower.tail=F, log.p=T)/log(10)]	
	summary.dat[, log10p:=as.numeric(round(log10p, digits=2), scientific=F)]

	matrix.final[, c(metabolite.ref):=summary.dat[match(matrix.final[,markername], markername), log10p]]

}

if(nrow(matrix.final[markername=="",])>0) matrix.final <- matrix.final[-which(markername==""),]
save(matrix.final, file="MATRIX_sub.RData")

if(nrow(matrix.final)>10000000000) {

cat("\nCreating correlation matrix...")
##Create a correlation matrix

	cor.matrix <- setnames(cbind(data.table(variant=matrix.final$markername), data.table(diag(1, nrow(matrix.final)))), c('variant', matrix.final$markername))
	for(var.one in cor.matrix$variant) {
		for(var.two in cor.matrix$variant) cor.matrix[match(var.one, variant), c(var.two) := list(cor(dosage.ref[,get(var.one)], dosage.ref[,get(var.two)]))]
	}

save(cor.matrix, file='Correlation_Matrix_AllVariants.RData')
###Create final LD pruned matrix

cat("\n\nPruning final matrix...")

matrix.prune <- NULL
for(index in 1:nrow(matrix.final)) matrix.prune <- rbind(matrix.prune, data.table(variant=matrix.final[index, markername], pvalue=unlist(matrix.final[index, 2:ncol(matrix.final), with=F]), metabolite=colnames(matrix.final)[2:ncol(matrix.final)]))


	metabolite.picked <- NULL
	create.matrix <- TRUE
	matrix.out <- NULL
	variants.ref <- NULL

	while(1) {

		matrix.prune <- matrix.prune[order(pvalue, decreasing=T), ]
		metabolite.ref <- matrix.prune[1, metabolite]
		metabolite.picked <- c(metabolite.picked, metabolite.ref)

			if(create.matrix) {
			
				matrix.out <- matrix.final[,c(1, match(metabolite.ref, colnames(matrix.final))),with=F]
				#matrix.out <- matrix.out[complete.cases(matrix.out),]
				matrix.out <- matrix.out[order(get(metabolite.picked), decreasing=T),]
				variants.ref <- matrix.out[complete.cases(matrix.out),markername]

			} else {
			
				variant.condition <- data.table(variant.ref=matrix.final[complete.cases(get(metabolite.ref)), markername])
				variant.condition[variant.ref%in%variants.ref==T, variant.replace:=variant.ref]

				if(nrow(variant.condition[is.na(variant.replace)==T,])>0) {
					
					for(variant.check in variant.condition[is.na(variant.replace)==T,variant.ref]) {

						cor.tmp <- data.table(cor(dosage.ref[,c(variant.check),with=F], dosage.ref[,c(variants.ref),with=F]))
						var.replace <- colnames(cor.tmp)[which(abs(cor.tmp[1,])>0.90)]

						if(length(var.replace)>0) variant.condition[variant.ref==variant.check, variant.replace:=colnames(cor.tmp)[which(cor.tmp[1,]==max(cor.tmp[1,]))]]
						if(length(var.replace)==0) variant.condition[variant.ref==variant.check, variant.replace:=variant.check]

					}
					
				}

				        formula.string <- as.formula(paste0(metabolite.ref, " ~ ", paste(c(variant.condition[,unique(c(variant.replace, variant.ref))], paste0("EPICN_PC_", 1:5), paste0("INTERVAL_PC_",1:5), "study"), collapse=" + ")))
					lm.dat <- lm(formula.string, data=dosage.ref)
					summary.dat <- cbind(data.table(markername=names(summary(lm.dat)$coefficients[,4])), data.table(summary(lm.dat)$coefficients))
					summary.dat <- summary.dat[-1,]
					summary.dat <- summary.dat[-grep("study|PC",markername),]

					setnames(summary.dat, c('markername', 'beta', 'se', 'tvalue', 'p'))

					summary.dat <- summary.dat[,log10p:=-pchisq((beta/se)^2, df=1, lower.tail=F, log.p=T)/log(10)]
					summary.dat[, log10p:=as.numeric(round(log10p, digits=2), scientific=F)]

						variant.condition[,variant.ref.beta:=summary.dat[match(variant.condition$variant.ref, markername), beta]]
						variant.condition[,variant.ref.log10p:=summary.dat[match(variant.condition$variant.ref, markername), log10p]]
						variant.condition[,variant.replace.beta:=summary.dat[match(variant.condition$variant.replace, markername), beta]]
						variant.condition[,variant.replace.log10p:=summary.dat[match(variant.condition$variant.replace, markername), log10p]]

						variant.condition[,index:=1:nrow(variant.condition)]
						setkey(variant.condition, index)
						variant.condition[,cor.val:=cor(dosage.ref[,variant.ref,with=F], dosage.ref[,variant.replace,with=F]), by=index]

					matrix.out[match(summary.dat$markername, markername), c(metabolite.ref):=summary.dat[,log10p]]
					matrix.out <- matrix.out[order(get(metabolite.picked), decreasing=T),]

					matrix.out[,na_count:=mapply(function(x) sum(is.na(matrix.out[x,])), 1:nrow(matrix.out))]
					variants.ref <- matrix.out[na_count!=(ncol(matrix.out)-2), markername]
					matrix.out[,na_count:=NULL]


			}
			matrix.prune <- matrix.prune[metabolite%in%metabolite.picked==F,]
			if(nrow(matrix.prune)==0) break;
			create.matrix <- FALSE

	}

	matrix.out[,na_count:=mapply(function(x) sum(is.na(matrix.out[x,])), 1:nrow(matrix.out))]
	matrix.out <- matrix.out[na_count!=(ncol(matrix.out)-2), ]
	matrix.out[,na_count:=NULL]	

matrix.pruned <- NULL
matrix.pruned <- matrix.out
save(matrix.pruned, file='MATRIX_pruned.RData')

}
##
reformat.matrix <- function(x.dat) {

	x.dat[,c(1:ncol(x.dat)):=lapply(.SD, as.character), .SDcols=c(1:ncol(x.dat))]
	for(metabolite.ref in colnames(x.dat)[grepl('^M', colnames(x.dat))]) {
	
		formula.string <- as.formula(paste0(metabolite.ref, " ~ ", paste(c(x.dat[is.na(get(metabolite.ref))==F, markername], paste0("EPICN_PC_", 1:5), paste0("INTERVAL_PC_",1:5), "study"), collapse=" + ")))
		lm.dat <- lm(formula.string, data=dosage.ref)
		summary.dat <- cbind(data.table(markername=names(summary(lm.dat)$coefficients[,4])), data.table(summary(lm.dat)$coefficients))
		summary.dat <- summary.dat[-1,]
		summary.dat <- summary.dat[-grep("EPICN_PC|INTERVAL_PC|study",markername),]

		setnames(summary.dat, c('markername', 'beta', 'se', 'tvalue', 'p'))

		summary.dat <- summary.dat[,log10p:=-pchisq((beta/se)^2, df=1, lower.tail=F, log.p=T)/log(10)]
		summary.dat[, log10p:=as.numeric(round(log10p, digits=2), scientific=F)]		

		for(marker.name in x.dat[,markername]) if(nrow(summary.dat[markername%in%marker.name==T,])==1) x.dat[markername==marker.name, c(metabolite.ref):=paste0(get(metabolite.ref), " (", summary.dat[markername==marker.name, round(beta, digits=4)], ", ", summary.dat[markername==marker.name, round(se, digits=4)], ")")]

	}	

	annotation <- fread('./metabolon_metabolomics/interval/gwas/reference_files/Metabolite_pathway_annotation.csv',header=T,sep=",")
	setkey(annotation, c_ln_comp_id)
	annotation[,comp_id:=mapply(function(x) paste0("M",paste(rep("0",(5-nchar(x))),collapse=""),x), c_ln_comp_id)]
	setnames(x.dat, annotation[match(colnames(x.dat), comp_id), BIOCHEMICAL])
	x.dat[,c('chr','pos'):=tstrsplit(markername, split="_", fixed=T)[1:2]]
	x.dat[,chr:=gsub('chr','',chr)]
	write(x.dat[,paste(chr,pos,sep=":")], 'snps_include.txt')
	sys.cmd <- paste0(qctool.exe, " -g interval_epicn.bgen -incl-positions snps_include.txt -snp-stats variants_pruned.stats 2> /dev/null")
	system(sys.cmd)

	snp.stats <- fread('variants_pruned.stats',header=T)

	write(x.dat[,paste(paste0("chr",chr),pos,sep="\t")], 'snps_include.txt')

		sys.cmd <- "/scratch/sk752/dbsnp_files/UCSC_b147/test_fetch/fetch_rsIDs.sh -i snps_include.txt -t snp147_hg19 -o output.tsv"
		dev.null <- system(sys.cmd, intern=T)
		rsid.info <- snp.stats[,.(chrom=paste0("chr",chromosome), chromEnd=position, name="-")]
		
		if(file.exists('output.txt')) rsid.info <- fread('output.tsv',header=T,sep="\t",showProgress=F)
		unlink('snps_include.txt')
		unlink('output.tsv')

	x.dat[,rsid:=rsid.info[match(x.dat[,paste(chr,pos,sep=":")], paste(gsub("chr","",chrom), chromEnd, sep=":")), name]]
	x.dat[,c('maf','info'):=list(snp.stats[match(x.dat$markername, gsub(":", "_", SNPID)), MAF], snp.stats[match(x.dat$markername, gsub(":", "_", SNPID)), information])]
	x.dat <- x.dat[,c(match(c('rsid','markername','chr','pos','maf','info'), colnames(x.dat)), 2:(ncol(x.dat)-5)), with=F]
	x.dat[,c(3:6):=lapply(.SD, as.numeric), .SDcols=c(3:6)]

	setnames(x.dat, c('rsID','MarkerName','Chromosome','Position','MAF','INFO', colnames(x.dat)[7:ncol(x.dat)])) 
	return(x.dat)

}

matrix.final.annotated <- NULL
matrix.final.annotated <- reformat.matrix(matrix.final)
save(matrix.final.annotated, file="MATRIX_sub_annotated.RData")

if(FALSE) {
matrix.pruned.annotated <- NULL
matrix.pruned.annotated <- reformat.matrix(matrix.pruned)
save(matrix.pruned.annotated, file="MATRIX_pruned_annotated.RData")

file.out <- paste0('Locus_',args[1],'_Summary.xlsx')
rm(matrix.final.annotated)

	load('MATRIX_sub_annotated.RData')
		matrix.final.annotated[ is.na(matrix.final.annotated)==T ] <- "-"
		write.xlsx(matrix.final.annotated, file.out, row.names=F, sheetName='Matrix_BeforePruning')

rm(matrix.pruned.annotated)
	load('MATRIX_pruned_annotated.RData')
		matrix.pruned.annotated[ is.na(matrix.pruned.annotated)==T ] <- "-"
		write.xlsx(matrix.pruned.annotated, file.out, row.names=F, sheetName='Matrix_AfterPruning', append=T)

rm(cor.matrix)
	load('Correlation_Matrix_AllVariants.RData')
		write.xlsx(cor.matrix, file.out, row.names=F, sheetName='Correlation_Matrix', append=T)

	cat(paste0("\nFinal matrices written to ", getwd(), "/", file.out))
	cat("\n\n")
}

cat("\n\nDone")
end.time <- Sys.time()
cat("\nEnded on ", format(Sys.time(), "%a %b %d %X %Y"))
cat(paste0("\nConditional analysis completed on ",Sys.info()[["nodename"]]," in ",round(as.numeric(end.time-start.time, units="mins"),digits=2), " minutes\n"))
rm(list=ls())














