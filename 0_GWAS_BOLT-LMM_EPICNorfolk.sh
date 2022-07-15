#!/bin/bash

# specificy shell to interpret the job script
#$ -S /bin/bash

# Redirect output and error to this location
#$ -o /scratch2/Metabolon_GWAS/06_gwas_HRCUK10K1000G/BOLT/output
#$ -e /scratch2/Metabolon_GWAS/06_gwas_HRCUK10K1000G/BOLT/output

#Set name of job
#$ -N BOLT

#send status information to this email address.
#$ -M Isobel.Stewart@mrc-epid.cam.ac.uk

#send an email when the job is done.
#$ -m e

#specify queue - you check what queues are available with qconf -sql.  qconf -sq <queue name> will give you details of the queue
#$ -q all.q

export BOLT=/genetics/bin/BOLT-LMM_v2.2/bolt
export DIR=/genetics/data/GWA/metabolomics/metabolon_epicn/06_gwas_HRCUK10K1000G
export OUTDIR=/scratch2/Metabolon_GWAS/06_gwas_HRCUK10K1000G/BOLT/results


TRAIT="$(cat $DIR/data/phenotypes/metabolitelist.txt | awk '{print $1}' | head -n ${SGE_TASK_ID} | tail -n 1 )"


$BOLT \
--bfile=/genetics/data/GWA/metabolomics/metabolon_epicn/02_model_validation_gwas/data/final_dataset/imputed_data/EPICN_OMICS_5841extract \
--geneticMapFile=/genetics/data/GWA/bolt/tables/genetic_map_hg19.txt.gz \
--LDscoresFile=/genetics/data/GWA/bolt/tables/LDSCORE.1000G_EUR.tab.gz \
--modelSnps=/genetics/data/GWA/metabolomics/metabolon_epicn/02_model_validation_gwas/data/EPICNorfolk_SNPskeep.txt \
--phenoFile=$DIR/data/phenotypes/epicn_metabolon_Apr2016_gwascohort_residuals_5841.txt \
--phenoCol=$TRAIT \
--lmm \
--impute2FileList=$DIR/data/Imputed/impute_filelist.txt \
--impute2FidIidFile=$DIR/data/Imputed/UK10K-HRC-EPIC.chr1_met5841_formatted.sample \
--statsFile=$OUTDIR/EPICNorfolk_BOLT_genotypedSNPs_${TRAIT}_HRC1KGUK10K \
--statsFileImpute2Snps=$OUTDIR/EPICNorfolk_BOLT_imputedSNPs_${TRAIT}_HRC1KGUK10K \
--numThreads 1 > /scratch2/Metabolon_GWAS/06_gwas_HRCUK10K1000G/BOLT/output/EPICNorfolk_BOLT_${TRAIT}_HRC1KGUK10K.log


gzip $OUTDIR/EPICNorfolk_BOLT_genotypedSNPs_${TRAIT}_HRC1KGUK10K
gzip $OUTDIR/EPICNorfolk_BOLT_imputedSNPs_${TRAIT}_HRC1KGUK10K

