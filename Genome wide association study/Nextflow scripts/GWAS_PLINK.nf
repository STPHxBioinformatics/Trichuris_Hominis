nextflow.enable.dsl=2
params.reference = "/scicore/home/schpie00/baer0006/Trichuris_Trichiura_CI/analysis_pipeline/results_vcf/Run_4_own_genome/"


log.info """\
	
	=============================================
	SWISS-TPH WHOLE-GENOME-PIPELINE: CAN'T RESIST 2
	=============================================

	"""
	.stripIndent()



\* Filter based on hwe once and also once without *\

process GWAS_QC {
  input:
  path(genoVCF)


  script:
  """
  plink --vcf ${genoVCF} --make-bed --out trichuris_data --double-id --allow-extra-chr --allow-no-sex
  plink --bfile trichuris_data --make-bed --out trichuris_pheno_data --pheno /scicore/home/schpie00/baer0006/Trichuris_Trichiura_CI/pairwise_stat_SNP/Main_Phenotype_File.txt --double-id --allow-extra-chr --allow-no-sex
  plink --bfile trichuris_pheno_data --make-bed --out trichuris_pheno_sex_data --update-sex /scicore/home/schpie00/baer0006/Trichuris_Trichiura_CI/pairwise_stat_SNP/Main_Phenotype_File_sex.txt --double-id --allow-extra-chr --allow-no-sex
  plink --bfile trichuris_pheno_sex_data --geno 0.05 --mind 0.1 --make-bed --out trichuris_pheno_sex_data_filt --double-id --allow-extra-chr --allow-no-sex
  plink --bfile trichuris_pheno_sex_data_filt --maf 0.05 --make-bed --out trichuris_pheno_sex_data_filt2 --double-id --allow-extra-chr --allow-no-sex
  plink --bfile trichuris_pheno_sex_data_filt2 --hwe 0.001 --make-bed --out trichuris_pheno_sex_data_filt_hwe --double-id --allow-extra-chr --allow-no-sex
  plink --bfile trichuris_pheno_sex_data_filt2 --pca --out data_pca --double-id --allow-extra-chr --allow-no-sex
  plink --bfile trichuris_pheno_sex_data_filt2 --assoc --adjust --out gw_results --double-id --allow-extra-chr --allow-no-sex
  plink --bfile trichuris_pheno_sex_data_filt_hwe --pca --out data_pca --double-id --allow-extra-chr --allow-no-sex
  plink --bfile trichuris_pheno_sex_data_filt_hwe --assoc --adjust --out gw_results --double-id --allow-extra-chr --allow-no-sex
  """
}

workflow {
	 vcf_ch = channel.fromPath(params.reference, checkIfExists: true)

	 GWAS_QC(vcf_ch)
}