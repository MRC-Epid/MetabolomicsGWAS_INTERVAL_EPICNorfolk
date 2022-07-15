#!/bin/R

rm(list=ls())
options(warn=-1)

closeAllConnections()
setwd(getwd())
source('/home/ps629/scripts/metabolomics/METABOLOMICS.info.R')
mrc.files <- "./mrc_files_for_meta"

args <- commandArgs(trailing=T)
metabolite.info <- fread("./reference_files/metabolites_reference.list")
metabolite.name <- metabolite.info[which(id==args[1]),metabolite]
fname.interval.list <- list.files(metal.input.path,full.names=T,recursive=T)[grep('.txt.gz',list.files(metal.input.path,full.names=T,recursive=T))]
fname.mrc.list <- list.files(mrc.files,full.names=T,recursive=T)[grep('MAformat.gz',list.files(mrc.files,full.names=T,recursive=T))]

unlink(metabolite.name, recursive=T)
dir.create(metabolite.name,showWarnings=F)
setwd(metabolite.name)

start.time <- Sys.time()

##Create metal script for the input files just created
metal.script.path <- getwd()
metal.outfile <- paste0(metal.script.path,"/",metabolite.name,"_METAL.sh")

unlink(metal.outfile)

	#filter.tab <- paste0("#ANALYSIS ",metabolite.name,"\n\nSEPARATOR TAB\nCOLUMNCOUNTING STRICT\nMARKER MARKERNAME\nALLELE EFFECT_ALLELE NON_EFFECT_ALLELE\nEFFECT BETA\nFREQLABEL EAF\nSTDERRLABEL SE\n#WEIGHTLABEL N\nPVALUE PVALUE\nSCHEME STDERR\nAVERAGEFREQ ON\nMINMAXFREQ ON\n#CUSTOMVARIABLE N\nCUSTOMVARIABLE MAC\nADDFILTER MAC > 10")
	filter.tab <- paste0("#ANALYSIS ",metabolite.name,"\n\nSEPARATOR TAB\nCOLUMNCOUNTING STRICT\nSCHEME STDERR\nAVERAGEFREQ ON\nMINMAXFREQ ON\nCUSTOMVARIABLE N\nCUSTOMVARIABLE MAC\nADDFILTER MAC > 10\nLABEL N AS N\nGENOMICCONTROL OFF\nMARKER MARKERNAME\nALLELELABELS EFFECT_ALLELE NON_EFFECT_ALLELE\nEFFECTLABEL BETA\nFREQLABEL EAF\nSTDERRLABEL SE\n")
	write(filter.tab,metal.outfile,append=T)
	fname.list <- c(fname.interval.list[grep(paste0("_",metabolite.name,"_"),fname.interval.list)],fname.mrc.list[grep(paste0("_",metabolite.name,"_"),fname.mrc.list)])
				
	for(fname in fname.list) {
		write(paste0("PROCESSFILE ",fname), metal.outfile,append=T)
	}
	write(paste0("\nOUTFILE ",getwd(), "/INTERVAL_MRC-Epi_",metabolite.name," .tbl\nANALYZE HETEROGENEITY\n\nCLEAR\n\n"),metal.outfile,append=TRUE)

	##Run METAL on the input files
	fname <- metal.outfile
	log.fname <- paste0(fread(paste0("egrep OUTFILE ",fname),header=F)[[2]],".tbl.log")
	metal.exe <- "/usr/local/Cluster-Apps/ceuadmin/metal_updated/metal-hiprec"

	metal.cmd <- paste0(metal.exe," ",fname," > ",log.fname)
	system(metal.cmd)


#for(fname in list.files()[grep('tbl',list.files())]) {
#	file.rename(fname,gsub("1.tbl",".tbl",fname))
#}

if(FALSE) {
####Get significant results from METAL
metal.res.path <- getwd()
###This list of genes was provided by mk863 based on the co-ordinates (1K recombination) from Joseph Pickerell (provided by James Staley)
gene.list <- fread('./reference_files/joseph_pickrell_recombination_information_gene_list_from_mk863.txt',header=F)
setnames(gene.list,c('genes','chr','start','end'))
gene.list[,locus:=paste(chr,start,end,sep=":")]

add.locus <- function(x) {

	chromosome <- as.numeric(gsub("chr","",unlist(strsplit(x,split=":"))[1]))
	position <- as.numeric(unlist(strsplit(x,split=":"))[2])
	return(gene.list[which(chr==chromosome & start<position & end>position),.(genes,locus)])

}

        #Get metal results and filter the signficant ones
        assoc.fname <- list.files(metal.res.path,full.names=T)[grep(paste0("_",metabolite.name,".tbl$"),list.files(metal.res.path,full.names=T))]
        if(length(assoc.fname)!=1 | !file.exists(assoc.fname)) stop("\nMore than one input file for the given metabolite or file not present\n\n")

        assoc <- fread(assoc.fname,header=T,showProgress=F)
        setnames(assoc,gsub("P-value","Pvalue",colnames(assoc)))
        assoc.out <- assoc[which(Pvalue<5e-08 & (Direction=="++" | Direction=="--")),]

                ###Get study results for the significant variants
                input.files <- fread(paste("grep 'Input File'",paste0(assoc.fname,".info")),header=F)[[7]]
                mrc.assoc <- fread(input.files[grep('MRC',input.files)],header=T,showProgress=F)
                interval.assoc <- fread(input.files[grep('INTERVAL',input.files)],header=T,showProgress=F)
                setnames(mrc.assoc, paste0("MRCEpi_",colnames(mrc.assoc)))
                setnames(interval.assoc, paste0("INTERVAL_",colnames(interval.assoc)))

                ##Combine all results

                assoc.out <- cbind(assoc.out,
                                mrc.assoc[match(assoc.out[,MarkerName],MRCEpi_MARKER),.(MRCEpi_REF_ALLELE,MRCEpi_OTHER_ALLELE,MRCEpi_FREQ,MRCEpi_EFFECT,MRCEpi_STDERR,MRCEpi_PVALUE)],
                                interval.assoc[match(assoc.out[,MarkerName],INTERVAL_MARKER),.(INTERVAL_REF_ALLELE,INTERVAL_OTHER_ALLELE,INTERVAL_FREQ,INTERVAL_EFFECT,INTERVAL_STDERR,INTERVAL_PVALUE)])

                rm(interval.assoc)
                rm(mrc.assoc)
                assoc.out[which(toupper(Allele1)!=MRCEpi_REF_ALLELE),
                                c("MRCEpi_REF_ALLELE","MRCEpi_OTHER_ALLELE","MRCEpi_FREQ","MRCEpi_EFFECT"):=list(toupper(Allele1),toupper(Allele2),1-MRCEpi_FREQ,-1*MRCEpi_EFFECT)]
                assoc.out[which(toupper(Allele1)!=INTERVAL_REF_ALLELE),
                                c("INTERVAL_REF_ALLELE","INTERVAL_OTHER_ALLELE","INTERVAL_FREQ","INTERVAL_EFFECT"):=list(toupper(Allele1),toupper(Allele2),1-INTERVAL_FREQ,-1*INTERVAL_EFFECT)]


                ##Get final list of variants with consistent direction of effects, significant at GWS and significant at a mininum p-value threshold of .01 in both studies
                assoc.final <- assoc.out[which(Pvalue<5e-08 & MRCEpi_PVALUE<0.01 & INTERVAL_PVALUE<0.01),]
		setkey(assoc.final,MarkerName)

		assoc.final[, c("genes","locus"):=add.locus(MarkerName), by=MarkerName]

                if(nrow(assoc.final)>0) {
                        assoc.final <- assoc.final[,.(MarkerName,toupper(Allele1),toupper(Allele2),Direction,HetPVal,Freq1,INTERVAL_FREQ,MRCEpi_FREQ,Effect,INTERVAL_EFFECT,MRCEpi_EFFECT,StdErr,INTERVAL_STDERR,MRCEpi_STDERR,Pvalue,INTERVAL_PVALUE,MRCEpi_PVALUE,genes,locus)]
                        setnames(assoc.final,c("MarkerName","Allele1","Allele2",colnames(assoc.final)[4:(ncol(assoc.final)-2)],"GenesInLocus","Locus"))
                        assoc.final[,Metabolite:=metabolite.name]
                        write.table(assoc.final,paste0("INTERVAL_MRC-Epi_",metabolite.name,"_significant.txt"),col.names=T,row.names=F,quote=F,sep="\t")

		}

}
###Done. Write runtime stats
end.time <- Sys.time()
cat("\nProcessed",metabolite.name,"on",Sys.info()[["nodename"]],"in",round(as.numeric(end.time-start.time, units="mins"),digits=2))

rm(list=ls())
		


	




	
