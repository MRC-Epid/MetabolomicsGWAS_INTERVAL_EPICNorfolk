#!/usr/bin/Rscript

rm(list=ls())
closeAllConnections()
system('clear')

library(data.table)
library(here)

source("/home/ps629/scripts/metabolomics/METABOLOMICS.info.R")
str_eval <- function(x) {return(eval(parse(text=x)))}

op <- options(digits.secs = 6)
args <- commandArgs(trailing=T)

#args <- c("2-3")
array.ids <- as.numeric(unlist(sapply(gsub("-",":",unlist(strsplit(args[1],split=","))),str_eval)))
user.name <- Sys.info()[["user"]]
job.name.list <- NULL

	##GET METABOLITE NAMES (COMP_IDS)

	metabolite.info <- fread("metabolites.list")
	metabolite.names <- metabolite.info[id%in%array.ids,metabolite]

	cat("\n\nPreparing reference files to run BOLM-LMM analysis for the metabolites listed below:\n\n")

	cat(paste0("\t",paste(metabolite.names,collapse=","),"\n\n"))

	cat("Preparing reference phenotype files...\n")
	pheno.ref <- fread('phenotype_subsetted.samples',showProgress=F,header=T,colClasses=c('character','character',rep('numeric',996)))
	pheno.ref <- pheno.ref[-1,]
	col.pick <- c('id_1','id_2',metabolite.names)
	
	pheno.out <- pheno.ref[,match(col.pick,colnames(pheno.ref)),with=F]
	setnames(pheno.out,c('FID','IID',metabolite.names))
	dir.create('reference_files',showWarnings=F)
	write.table(pheno.out,'./reference_files/BOLT.pheno',col.names=T,row.names=F,quote=F,sep=" ")

	cat("Getting list of variants for model building...\n")
	##Get the set of files required to analysed genotyped variants
	
	genotyped <- gsub(".fam","",list.files(genotyped.path,full.names=T)[grep('fam$',list.files(genotyped.path,full.names=T))])
	pheno.fname <- list.files(getwd(),full.names=T,recursive=T)[grep('/BOLT.pheno$',list.files(getwd(),full.names=F,recursive=T))]
	model.pheno.file <- "/scratch/curated_genetic_data/reference_files/interval/interval_bolt_modelbuilding_2.2.16"
	write(fread(paste0(model.pheno.file,".bim"),header=F,showProgress=F)[[2]],'./reference_files/BOLT_ModelVariants.list')
	modelvar.fname <- list.files(getwd(),full.names=T,recursive=T)[grep('reference_files/BOLT_ModelVariants.list$',list.files(getwd(),full.names=F,recursive=T))]
	
	cat("Submitting jobs to analyse genotyped variants and imputed variants for the metabolites listed above\n")

	##Submite jobs to analyse each metabolite
	
	for(metabolite.name in metabolite.names) {

	  slurm.id.list <- NULL
	  temp.folder.name <- gsub(":|\\.| ","-",gsub(" BST","",Sys.time()))
	  
	  dir.create(temp.folder.name,showWarnings=FALSE)
	  start.time <- Sys.time()
	  
	  #Job name for this process pipeline
	  workdir.path <- paste0(getwd(),"/",temp.folder.name)
	  setwd(workdir.path)
	  
	  ##RUN BOLT-LMM for genotyped variants
	  
	  bolt.exe <- "/usr/local/Cluster-Apps/ceuadmin/boltlmm/2.2/bolt"	
	  job.name <- paste0(metabolite.name,'_Genotyped')
	  job.name.list <- unique(c(job.name.list,job.name))
	  slurm.outfile <- paste0(metabolite.name,'_genotyped_SLURM.out')
	  genotyped.res.fname <- paste0(metabolite.name,"_BOLT-LMM_genotyped.assoc")
	  slurm.script <- "genotyped_imputed.slurm"
	  
	  bolt.cmd <- paste(bolt.exe,
	                    " --bed=",paste0(genotyped,".bed"),
	                    " --bim=",paste0(genotyped,".bim"),
	                    " --fam=",paste0(genotyped,".fam"),
	                    " --phenoFile=",pheno.fname,
	                    " --phenoCol=",metabolite.name,
	                    #" --covarFile=",covar.fname,
	                    #" --qCovarCol=COV_{1:3}",
	                    #" --covarCol=COV_19",
	                    " --lmmForceNonInf",
	                    " --modelSnps=",modelvar.fname,
	                    " --maxMissingPerSnp=1",
	                    " --numThreads=4",
	                    " --LDscoresMatchBp",
	                    " --LDscoresFile=/usr/local/Cluster-Apps/ceuadmin/boltlmm/2.2/tables/LDSCORE.1000G_EUR.tab.gz",
	                    " --geneticMapFile=/usr/local/Cluster-Apps/ceuadmin/boltlmm/2.2/tables/genetic_map_hg19.txt.gz",
	                    " --verboseStats",
	                    " --statsFile=",genotyped.res.fname,
	                    sep="")
	  system('cp /home/ps629/scripts/metabolomics/BOLT_SLURM_TEMPLATE.slurm genotyped_imputed.slurm')
	  write(bolt.cmd,slurm.script,append=T)
	  
	  #dev.null <- system(paste0('sbatch -J ',job.name,' --workdir=',workdir.path,' --output=',slurm.outfile,' ',genotyped_imputed.slurm),intern=T)
	  
	  ##RUN analysis for imputed variants
	  
	  #Get the path to the executables and module headers
	  
	  bolt.exe <- "/usr/local/Cluster-Apps/ceuadmin/boltlmm/2.2/bolt"
	  job.name <- paste0(metabolite.name,'_Imputed')
	  job.name.list <- unique(c(job.name.list,job.name))
	  imputed.res.fname <- paste0(metabolite.name,"_BOLT-LMM_imputed_chr_$SLURM_ARRAY_TASK_ID\\.assoc")
	  
	  bolt.cmd <- paste(bolt.exe,
	                    " --bfile=",model.pheno.file,
	                    " --phenoFile=",pheno.fname,
	                    " --phenoCol=",metabolite.name,
	                    #" --covarFile=",covar.fname,
	                    #" --qCovarCol=COV_{1:4}",
	                    #" --covarCol=COV_{5:6}",
	                    " --lmmForceNonInf",
	                    " --modelSnps=",modelvar.fname,
	                    " --maxMissingPerSnp=1",
	                    " --numThreads=4",
	                    " --LDscoresFile=/usr/local/Cluster-Apps/ceuadmin/boltlmm/2.2/tables/LDSCORE.1000G_EUR.tab.gz",
	                    " --geneticMapFile=/usr/local/Cluster-Apps/ceuadmin/boltlmm/2.2/tables/genetic_map_hg19.txt.gz",
	                    " --statsFile=/dev/null",
	                    " --bgenFile=",paste0(imputed.path,"/impute_$SLURM_ARRAY_TASK_ID\\_interval.bgen"),
	                    " --sampleFile=/scratch/curated_genetic_data/genotype_data/imputed/interval/interval.samples",
	                    " --statsFileBgenSnps=",imputed.res.fname,
	                    sep="")
	  
	  system('cp /home/ps629/scripts/metabolomics/BOLT_SLURM_TEMPLATE.slurm imputed.slurm')
	  write(bolt.cmd,'imputed.slurm',append=T)
	  
	  dev.null <- system(paste0('sbatch -J ',job.name,' --array=1-22 --workdir=',workdir.path,' imputed.slurm'),intern=T)
	  
	  setwd('..')	  
	
	 }

	job.info <- fread(paste0("squeue -u ",user.name," -p all -l --format='%.7A %.8u %.82j'"),skip=1)
	dependency.string <- paste(job.info[grep(paste(unique(job.name.list),collapse="|"),NAME),JOBID],collapse=":")
	sacct.string <- paste0("sacct --format 'JobID,jobname%20,AveVMSize,MaxVMSize,elapsed,state' | egrep -v batch | egrep '",paste(job.info[grep(paste(unique(job.name.list),collapse="|"),NAME),JOBID],collapse="|"),"'")
	
	system('cp /home/ps629/scripts/metabolomics/sendmail.slurm .')
	write(sacct.string,'sendmail.slurm',append=T)
	
	##Write a bit of bash code to move various files to a single directories
	bash.code <- paste0("\nmkdir -p bolt_output;\n",
                      "for fname in `find /rds/user/ps629/rds-jmmh2-projects/metabolon_metabolomics/interval/gwas/bolt_analysis -name *.assoc`\n",
	                    "do\n\t",
	                    "IFS='/' read -r -a array <<< \"$fname\"\n",
	                    "\tmv $fname ./bolt_output/${array[((${#array[@]}-1))]}\n",
	                    "done\n")
	write(bash.code,'sendmail.slurm',append=TRUE)
	
	bash.code <- paste0("\nmkdir -p bolt_slurm_output;\n",
	                    "for fname in `find /rds/user/ps629/rds-jmmh2-projects/metabolon_metabolomics/interval/gwas/bolt_analysis -name *.out`\n",
	                    "do\n\t",
	                    "IFS='/' read -r -a array <<< \"$fname\"\n",
	                    "\tmv $fname ./bolt_slurm_output/${array[((${#array[@]}-1))]}\n",
	                    "done\n")
	write(bash.code,'sendmail.slurm',append=TRUE)
	
	bash.code <- paste0("\ndir_name=`echo \"$(date)\" | awk '{gsub(/ /,\"_\");}1'`\n",
	                    "mkdir $dir_name\n",
	                    "mv bolt_slurm_output $dir_name/\n",
	                    "mv bolt_output $dir_name/\n",
	                    "mv Jobs_RunStats.report $dir_name/\n",
	                    "mv reference_files $dir_name/\n",
	                    "rm -rf 2016* *.slurm\n")
	
	write(bash.code,'sendmail.slurm',append=TRUE)
	
	dev.null <- system(paste0("sbatch --dependency=afterany:",dependency.string," sendmail.slurm"),intern=T)
	
	cat("Jobs have been submitted. SLURM scheduler will send an e-mail once the jobs are complete and a file that contains job report - 'Jobs_RunStats.report' will be created withing the working directory\n")
	cat("EOF\n\n")
	
rm(list=ls())
	

