#!/bin/R

rm(list=ls())
options(warn=-1)

closeAllConnections()
setwd("./files_to_share_final")

start.time <- Sys.time()

##Functions

getAssoc <- function(metabolite.name) {

        metabolite.info <- fread("./reference_files/metabolites.list")
        system('clear')

	interval.path <- paste0(interval.bolt.path,"/",metabolite.name,"_GWAS_BOLT")

        cat(paste0("\n\nMetabolite to be processed - ",metabolite.name,"\n"))
        fname.list <- list.files(interval.path,full.names=T)[grep(paste0(metabolite.name,'_BOLT-LMM_imputed.*.assoc$'),list.files(interval.path,full.names=T))]
	n.pheno <- 2*nrow(pheno.ref[complete.cases(get(metabolite.name)),])

        assoc <- NULL
        cat(paste0("\nTotal number of INTERVAL files to merge : ",length(fname.list)))
        cat("\nMerging...\n")
        index <- 1

        for(fname in fname.list) {

                cat("\r",paste0("File ",index))
                index <- index+1

                temp.infile <- suppressWarnings(fread(fname,showProgress=F))

		temp.infile[,MAF:=A1FREQ]
		temp.infile[which(A1FREQ>0.5), MAF:=1-A1FREQ]

		temp.infile[,MAC:=MAF*n.pheno]

                assoc <- rbind(assoc,temp.infile[complete.cases(SE),.(CHR,BP,ALLELE1,ALLELE0,N=n.pheno/2,A1FREQ,MAC,BETA,SE,P_BOLT_LMM_INF,INFO)])
                temp.infile <- NULL

        }

	cat("\nMerging imputed results complete\n")
        return(assoc)

}

###Get commandline arguments

	args <- commandArgs(trailing=TRUE)

	metabolite.info <- fread("./reference_files/metabolites.list")
	pheno.ref <- fread("./bolt_final/reference_files/BOLT.pheno",showProgress=F,header=T,colClasses=c('character','character',rep('numeric',995)))

	metabolite.name <- metabolite.info[which(id==args[1]),metabolite]
	fcount <- fread("./bolt_final/file_report/metabolite_filecount.txt",header=T)

	if(fcount[metabolite.name==metabolite,n_files]!=23) stop("File count not equal to 23.. please check")

	#args <- c('metabolite=M00054')
	#if(length(args)!=1) stop("\nPlease enter a metabolite name in the format\n\tmetabolite=[metabolite-name]\n")
	#if(length(grep('metabolite=',args))!=1) stop("\nPlease enter a metabolite name in the format\n\tmetabolite=[metabolite-name]\n")
		#metabolite.name <- gsub("metabolite=","",args[grep('metabolite=',args)])

	#Get combined association summary for metabolite from INTERVAL (BOLT-LMM results)
	assoc <- getAssoc(metabolite.name)
		
		cat("\nFormatting INTERVAL\n")
		##Create a marker name column
		setnames(assoc,c("chromosome","position","alleleA","alleleB","N","frequency_A","MAC","beta","se","pvalue","info"))
		assoc <- assoc[-which(info<0.3 | abs(beta)>10),]
		cat("Processing imputed data...\n")
		##Create a marker name column
                snps.index <- which(nchar(assoc[,alleleA])==1 & nchar(assoc[,alleleB])==1)
                indels.index <- which(nchar(assoc[,alleleA])>1 | nchar(assoc[,alleleB])>1)
                assoc[snps.index,name:=paste(chromosome,position,alleleA,alleleB,sep=":")]
                assoc[intersect(snps.index,which(alleleA>alleleB)),name:=paste(chromosome,position,alleleB,alleleA,sep=":")]
                assoc[indels.index,name:=paste(chromosome,position,alleleA,alleleB,sep=":")]
                assoc[intersect(indels.index,which(nchar(alleleB)>nchar(alleleA))),name:=paste(chromosome,position,alleleB,alleleA,sep=":")]
                assoc[,name:=paste0("chr",name)]

		imputed.assoc <- assoc
		assoc <- NULL

		cat("Processing genotyped data...\n")
		temp.infile <- fread(paste0(interval.bolt.path,"/",metabolite.name,"_GWAS_BOLT/",metabolite.name,"_BOLT-LMM_genotyped.assoc"),header=T,showProgress=F)
		n.pheno <- 2*nrow(pheno.ref[complete.cases(get(metabolite.name)),])

		temp.infile[,MAF:=A1FREQ]
                temp.infile[which(A1FREQ>0.5), MAF:=1-A1FREQ]
                temp.infile[,MAC:=MAF*n.pheno]

		assoc <- rbind(assoc,temp.infile[complete.cases(SE),.(CHR,BP,ALLELE1,ALLELE0,N=n.pheno/2,A1FREQ,MAC,BETA,SE,P_BOLT_LMM_INF)])
		setnames(assoc,c("chromosome","position","alleleA","alleleB","N","frequency_A","MAC","beta","se","pvalue"))
		rm(temp.infile)
		
		##Create a marker name column
                snps.index <- which(nchar(assoc[,alleleA])==1 & nchar(assoc[,alleleB])==1)
                indels.index <- which(nchar(assoc[,alleleA])>1 | nchar(assoc[,alleleB])>1)
                assoc[snps.index,name:=paste(chromosome,position,alleleA,alleleB,sep=":")]
                assoc[intersect(snps.index,which(alleleA>alleleB)),name:=paste(chromosome,position,alleleB,alleleA,sep=":")]
                assoc[indels.index,name:=paste(chromosome,position,alleleA,alleleB,sep=":")]
                assoc[intersect(indels.index,which(nchar(alleleB)>nchar(alleleA))),name:=paste(chromosome,position,alleleB,alleleA,sep=":")]
                assoc[,name:=paste0("chr",name)]

		genotyped.assoc <- assoc
		assoc <- NULL

		cat("Creating final dataset based on the rules described in the analysis plan\n")	
		##Create final data table to share with MRC-Epi
		res.out <- NULL

			#Get list of variants only in the imputed data
			var.pick <- setdiff(imputed.assoc[,name],genotyped.assoc[,name])
			res.out <- rbind(res.out,imputed.assoc[match(var.pick,name),])

			#Remove all variants unique to imputed dataset
			imputed.assoc <- imputed.assoc[-match(var.pick,name),]

			##Add well imputed genotyped variants to final data to share
			res.out <- rbind(res.out,imputed.assoc[info>0.6,])
			var.remove.from.genotyped <- intersect(imputed.assoc[info>0.6,name],genotyped.assoc[,name])
			genotyped.assoc <- genotyped.assoc[-match(var.remove.from.genotyped,name),]

			#Add genotyped variants with non "-" alleles to the final data to share
			var.pick <- genotyped.assoc[which(alleleA!="-" & alleleB!="-"),name]
			res.out <- rbind(res.out,data.table(genotyped.assoc[match(var.pick,name),.(chromosome,position,alleleA,alleleB,N,frequency_A,MAC,beta,se,pvalue)],info="-",genotyped.assoc[match(var.pick,name),.(name)]))

			genotyped.assoc <- genotyped.assoc[-match(var.pick,genotyped.assoc[,name]),]

			##Get the right alleles from Affy manifest that Mihir updated
			affy.manifest <- fread("./files_to_share_final/Modified_List.txt",header=T)
			common.var <- intersect(genotyped.assoc[,paste(chromosome,position,sep=":")],affy.manifest[,paste(chromosome,position,sep=":")])

			affy.manifest <- affy.manifest[match(common.var,paste(chromosome,position,sep=":")),]
			genotyped.assoc <- genotyped.assoc[match(common.var,paste(chromosome,position,sep=":")),]

			assoc <- NULL
			assoc <- cbind(genotyped.assoc,affy.manifest[,.(Mod_Ref,Mod_Alt)])

			##Create updated name for the final list of genotyped INDELS
			indels.index <- which(nchar(assoc[,Mod_Ref])>1 | nchar(assoc[,Mod_Alt])>1)
			assoc[indels.index,name_mod:=paste(chromosome,position,Mod_Ref,Mod_Alt,sep=":")]
			assoc[intersect(indels.index,which(nchar(Mod_Alt)>nchar(Mod_Ref))),name_mod:=paste(chromosome,position,Mod_Alt,Mod_Ref,sep=":")]
			assoc[,name_mod:=paste0("chr",name_mod)]

			genotyped.assoc <- assoc[,.(chromosome,position,alleleA,alleleB,N,frequency_A,MAC,beta,se,pvalue,name=name_mod)]
			if(length(intersect(genotyped.assoc[,name],res.out[,name]))>0) {
				genotyped.assoc <- genotyped.assoc[-which(name%in%intersect(genotyped.assoc[,name],res.out[,name])==TRUE),]
			}
			res.out <- rbind(res.out,data.table(genotyped.assoc[,.(chromosome,position,alleleA,alleleB,N,frequency_A,MAC,beta,se,pvalue)],info="-",genotyped.assoc[,.(name)]))

			rm(assoc,genotyped.assoc,imputed.assoc)
			##Finally change insertions and deletions to D/I description
			res.out.indels <- res.out[which(nchar(alleleA)>1 | nchar(alleleB)>1 | alleleA=="-" | alleleB=="-"),]
			res.out <- res.out[-which(nchar(alleleA)>1 | nchar(alleleB)>1 | alleleA=="-" | alleleB=="-"),]

			index.check <- which(nchar(res.out.indels[,alleleA])>nchar(res.out.indels[,alleleB]))
			res.out.indels[index.check,alleleA:="I"]
			res.out.indels[index.check,alleleB:="D"]

			index.check <- which(nchar(res.out.indels[,alleleA])<nchar(res.out.indels[,alleleB]))
			res.out.indels[index.check,alleleA:="D"]
			res.out.indels[index.check,alleleB:="I"]

			index.check <- which(res.out.indels[,alleleA]=="-")
			res.out.indels[index.check,alleleA:="D"]
                        res.out.indels[index.check,alleleB:="I"]

			index.check <- which(res.out.indels[,alleleB]=="-")
                        res.out.indels[index.check,alleleA:="I"]
                        res.out.indels[index.check,alleleB:="D"]

			index.d <- which(res.out.indels[,alleleA]=="D")
			res.out.indels[index.d,frequency_A:=1-frequency_A]
			res.out.indels[index.d,beta:=beta*-1]
			
			res.out.indels[index.d,alleleA:="I"]
			res.out.indels[index.d,alleleB:="D"]

			res.out.indels <- res.out.indels[-which(nchar(alleleA)>1 | nchar(alleleB)>1),]
			res.out <- rbind(res.out,res.out.indels)

			res.final <- res.out[,.(name,alleleA,alleleB,N,frequency_A,MAC,beta,se,pvalue)]
			res.out <- NULL

			setnames(res.final, c("MARKERNAME","EFFECT_ALLELE","NON_EFFECT_ALLELE","N","EAF","MAC","BETA","SE","PVAL"))
			cat("\nWriting the final output to file..")
			write.table(res.final,paste0("INTERVAL_",metabolite.name,"_formattedForMeta.txt"),col.names=T,row.names=F,quote=F,sep="\t")
			cat("\nCompressing the file...")			
			system(paste0("gzip ",paste0("INTERVAL_",metabolite.name,"_formattedForMeta.txt")))


		cat("Done\nMetal input files for metabolite",metabolite.name,"are in the directory ",getwd(),"\n")

###Done. Write runtime stats
end.time <- Sys.time()
cat("\nProcessed",metabolite.name,"on",Sys.info()[["nodename"]],"in",round(as.numeric(end.time-start.time, units="mins"),digits=2),"\n")
rm(list=ls())
		


	




	
