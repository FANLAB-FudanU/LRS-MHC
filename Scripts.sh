#!/bin/bash

##1. Reads mapping##
##1.1 LRS data##

meryl count k=15 output merylDB ${hg38_ref_genome}

meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt

winnowmap -W repetitive_k15.txt -L --cs --MD -t 16 -ax map-ont ${hg38_ref_genome} ${sample}.fq.gz > ${sample}.sam

##1.2 SRS data##
##1.2.1 Zhongyuan cohort##

bwa-mem2 mem -t 16 -B 4 -O 6 -E 1 -M -R ${tags} ${hg38_ref_genome} ${sample}_in1.fastq.gz ${sample}_in2.fastq.gz > ${sample}.sam

##1.2.2 the independent Han Chinese cohort##

bwa mem -t 16 -B 4 -O 6 -E 1 -M -R ${tags} ${hg38_ref_genome} ${sample}_in1.fastq.gz ${sample}_in2.fastq.gz > ${sample}.sam

##2. Coverage estimation##
##2.1 whole genome##

mosdepth -x -n -Q 20 -t 4 -b 5000 ${sample}_output ${sample}.bam

##2.2 MHC region##

samtools coverage -q 20 -r ${region_of_MHC/HLA_genes} ${sample}

##3. HLA typing##
##3.1 LRS data##

perl HLA-LA.pl --workingDIR ${outDIR} --sampleID ${sample} --BAM ${sample}.bam --longReads ont2d --graph PRG_MHC_GRCh38_withIMGT --maxThreads 4

python -m NomenCleaner --2field --hped ${hped_file1} --out ${output1}

##3.2 SRS data##

perl HLA-LA.pl --workingDIR ${outDIR} --sampleID ${sample} --BAM ${sample}.bam --graph PRG_MHC_GRCh38_withIMGT --maxThreads 4

python -m NomenCleaner --2field --hped ${hped_file2} --out ${output2}

##4. SNP calling##
##4.1 LRS data##

singularity exec --bind ${HOME} pepper_deepvariant_r0.7.sif run_pepper_margin_deepvariant call_variant \
-b ${sample}.bam -f ${hg38_ref_genome} -o ${OUTPUT_DIR} \
-p ${OUTPUT_PREFIX} -t 40 --ont_r9_guppy5_sup --gvcf

singularity exec --bind ${HOME} glnexus_vl.4.1.sif glnexus_cli \
--config DeepVariantWGS --thread 10 --list ${gvcf_files} >> lrs_jointcalling.bcf

whatshap genotype --chromosome ${chr} --reference ${hg38_ref_genome} \
-o ${sample}.${chr}.regenotype.vcf ${sample}.${chr}.jointcalling.vcf ${sample}.bam

##4.2 SRS data##

java -jar picard.jar MarkDuplicates \
--INPUT ${sample}.bam --OUTPUT ${sample}.dedup.bam \
--METRICS_FILE ${sample}.metric.txt \
--ASSUME_SORT_ORDER coordinate \
--VALIDATION_STRINGENCY SILENT \
--CREATE_INDEX true

gatk ApplyBQSR -R ${ref} -I ${dedup_bam} \
-bqsr ${recal_table} -O ${recal_bam}

gatk HaplotypeCaller -ERC GVCF --intervals ${chr} \
-R ${hg38_ref_genome} --dbsnp dbsnp_151.GRCh38.vcf \
-I ${recal_bam} -O ${sample}.${chr}.g.vcf.gz

gatk CombineGVCFs -R ${hg38_ref_genome} -V ${gvcfs_chr} \
-O combined.${chr}.g.vcf.gz --intervals ${chr}

gatk GenotypeGVCFs -R ${hg38_ref_genome} \
-V combined.${chr}.g.vcf.gz -O genotyped.${chr}.vcf.gz \
--dbsnp dbsnp_151.GRCh38.vcf --intervals ${chr}

##two steps: SNPs and InDels separately##

gatk ApplyVQSR -R ${hg38_ref_genome} --mode SNP\
-V ${genotyped_whole_genome}.vcf.gz \
-O ${appliedVQSR}.snps.vcf.gz \
--recal-file ${snp_recal_file} \
--truth-sensitivity-filter-level 99.0 \
--tranches-file ${snp_tranches_file} \
--seconds-between-progress-updates 30

gatk ApplyVQSR -R ${hg38_ref_genome} --mode INDEL\
-V ${appliedVQSR}.snps.vcf.gz \
-O ${appliedVQSR}.snps.indels.vcf.gz \
--recal-file ${indel_recal_file} \
--truth-sensitivity-filter-level 99.0 \
--tranches-file ${indel_tranches_file} \
--seconds-between-progress-updates 30

bcftools view -f PASS -m2 -M2 -v snps,indels -O z -o ${filtered}.vcf.gz ${appliedVQSR}.snps.indels.vcf.gz

##5. Variants annotation##

vcftools --vcf ${whole_genome_vcf} --chr 6 --from-bp 25000000 --to-bp 34000000 --recode --recode-INFO-all --out ${MHC_vcf}

vcftools --vcf ${MHC_vcf} --min-alleles 2 --max-alleles 2 --max-missing 0.99 --maf 0.01 --hwe 1e-7 --recode --recode-INFO-all --out ${vcf_out}

perl convert2annovar.pl -allsample -withfreq -format vcf4 ${input}.vcf > ${input}.avinput

perl annotate_variation.pl -out ${input} -build hg38 ${input}.avinput humandb

perl table_annovar.pl ${input}.avinput humandb -buildver hg38 -out ${input} \
-remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a -operation gx,r,f,f,f -nastring . -csvout -polish

##6. Construction of LRS-based MHC reference panel##

liftover ${oldfile} ${chain_$hg38_to_hg19_$hg19_to_hg18} ${mapped_new_file} ${unmapped_postions_file}

tcsh MakeReference.csh ${bfile_hg18} ${HLA_allele_file} ${mhc_panel} plink

##7. Imputation to the independent Han Chinese cohort##
##7.1 HLA imputation##

python -m MakeGeneticMap -i ${bfile} -hg 18 -ref ${mhc_panel} -o ${tmp}

python CookHLA.py -hg 18 -ref ${mhc_panel} -o ${output} -gm ${tmp}.mach_step.avg.clpsB \
-ae ${tmp}.aver.erate -mem 200g -bgl4 -i ${bfile}

python -m bMarkerGenerator --bfile ${bfile} --chped ${chped} --hg 18 \
--dict-AA HLA_DICTIONARY_AA.hg18.imgt3320 --dict-SNPS HLA_DICTIONARY_SNPS.hg18.imgt3320 --out ${output}.bmarker

##7.2 SNP imputation##

eagle --bfile ${target_bfile} --geneticMapFile genetic_map_hg18_withX.txt.gz --chrom 6 --outPrefix ${target_phased} --numThreads 8

minimac4 --refHaps ${mhc_panel_M3VCF_format} --haps ${target_phased}.vcf \
--chr 6 --start 25000000 --end 34000000 --format GT --allTypedSites --prefix ${target_minimac}

##8. PheWAS and conditional analysis##

plink2.0 --bfile ${bfile} --pheno ${phenotypes} --covar ${covariate} --glm --quantile-normalize --out ${output}

gcta --bfile ${target} --chr 6 --maf 0.01 --cojo-file ${summary_statistics_of_each_phenotype}.ma --cojo-slct --out ${cojo_out}

##9. Fine-mapping##

finemap --sss --in-files ${master_file} --n-threads 10 --n-causal-snps ${num_of_independent_signals_GCTA-COJO} --dataset ${1..N}
