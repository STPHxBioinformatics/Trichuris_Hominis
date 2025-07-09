nextflow.enable.dsl=2

params.whole_genome_assembly_freeze = "/scicore/home/schpie00/baer0006/ONT/de_novo_annotation_pipeline/trichuris_suis_PRJNA208416.fasta"
params.prothint_output = "/scicore/home/schpie00/baer0006/ONT/de_novo_annotation_pipeline/prothint_augustus.gff"
params.protein_db_metazoa = "/scicore/home/schpie00/baer0006/ONT/de_novo_annotation_pipeline/Metazoa.fa"

log.info """\

	=============================================
	SWISS-TPH WHOLE-GENOME-PIPELINE: CAN'T RESIST
	=============================================

	Babedibubbi

	"""
	.stripIndent()

workflow annotation {
  repeatmasker_input_ch = Channel.fromPath('/scicore/home/schpie00/baer0006/ONT/de_novo_annotation_pipeline/whole_genomes/*.fa')
  repeatmasker(repeatmasker_input_ch)
  repeatmasker_output_ch = repeatmasker.out

  braker2(repeatmasker_output_ch)
  braker2_ch = braker2.out

  gffread(braker2_ch)
  stats(braker2_ch)
}

workflow {
  annotation()
}

/*
Masking repeat regions with Repeatmarker
*/
process repeatmasker {
  input:
  path(genome_assembly)

  output:
  path("${genome_assembly}_masked.fasta")

  script:
  """
  mkdir run1 ;
  ml load RepeatMasker/4.1.4-foss-2021a ;
  perl  /scicore/soft/apps/RepeatMasker/4.1.4-foss-2021a/RepeatMasker -species nematoda -pa 8 -xsmall -dir ./run1 ${genome_assembly} ;
  cp run1/${genome_assembly}.masked ./${genome_assembly}_masked.fasta
  """
}

/*
Run ProtHint on assembly
*/
process prothint {
	input:
	path(repeat_masked_assembly)

	output:
	path("prothint_augustus.gff")
	
	script:
	"""
	ml purge;
	ml load Python/3.9.5-GCCcore-10.3.0;
  ml load Perl/5.32.1-GCCcore-10.3.0;
  ml load DIAMOND/2.0.15-GCC-10.3.0;
  source $HOME/venv_prothint/bin/activate;
  export PERL5LIB="$HOME/perl5/lib/perl5:$PERL5LIB";
  export PATH="$HOME/perl5/bin:$PATH";

  /scicore/home/schpie00/baer0006/packages/ProtHint/bin/prothint.py \
  --threads 16 \
  ${repeat_masked_assembly} \
  ${params.protein_db_metazoa};

  cp ProtHint/prothint_augustus.gff ./
  """
}

/*
Run Braken3 on assembly
--prot_seq=${params.protein_db_metazoa} \
  --hints=${params.prothint_output} \
*/
process braker3 {
	input:
	path(repeat_masked_assembly)

	output:
	path("braker3.gff")

	script:
	"""
	source $HOME/venv_braker/bin/activate;

	singularity exec /scicore/home/schpie00/baer0006/packages/braker3_latest.sif braker.pl \
  --threads=16 \
  --species=TrichurisSuis --useexisting \
  --genome=${repeat_masked_assembly} \
  --hints=${params.prothint_output} \
  --rnaseq_sets_ids=SRR1041650,SRR1041649,SRR1041648,SRR1041647,SRR1041646,SRR1041645 \
  --rnaseq_sets_dirs=/scicore/home/schpie00/baer0006/reference_RNA_seq_data/suis/ \
  --gff3 \
   --workingdir=. \
   ;

  cp braker/braker.gff3 ./braker3.gff
  
	"""
}

/*
Run Braken2 on assembly
  --prot_seq=${params.protein_db_metazoa} \
  	source $HOME/venv_braker/bin/activate;
  ml load AUGUSTUS/3.4.0-foss-2021a;
*/
process braker2 {
	input:
	path(repeat_masked_assembly)

	output:
	tuple path("${repeat_masked_assembly}"), path("braker2.gff")

	script:
	"""
	export TMPDIR=/tmp/;
	singularity exec /scicore/home/schpie00/baer0006/packages/braker3_latest.sif braker.pl \
  --threads=16 \
  --genome=${repeat_masked_assembly} \
  --prot_seq=${params.protein_db_metazoa} \
  --gff3  --species=trichuris  --useexisting;
  cp braker/braker.gff3 ./braker2.gff
	"""
}

process gffread {
  publishDir "/scicore/home/schpie00/baer0006/ONT/de_novo_annotation_pipeline/gff_files/", mode: 'copy'

  input:
  tuple path(repeat_masked_assembly), path(gff_file)

  output:
  path("${repeat_masked_assembly}gene_sequences.fasta")

  script:
  """
  ml load BEDTools/2.30.0-GCC-10.3.0;
  grep "gene" ${gff_file} > genes.gff;
  awk 'BEGIN {OFS="\\t"} \$3 == "gene" {split(\$9, a, ";"); split(a[1], b, "="); \$3 = b[1]"="b[2]; print}' genes.gff > output_file.gff;
  sed -E '/^>/ s/ /_/g' ${repeat_masked_assembly} > ${repeat_masked_assembly}_clean_header.fasta;
  bedtools getfasta -fi ${repeat_masked_assembly}_clean_header.fasta -bed output_file.gff -fo ${repeat_masked_assembly}gene_sequences.fasta -name
  """
}

process stats {
  publishDir "/scicore/home/schpie00/baer0006/ONT/de_novo_annotation_pipeline/stats_genome/", mode: 'copy'

  input:
  tuple path(repeat_masked_assembly), path(gff_file)

  output:
  tuple path("${repeat_masked_assembly}_report.tsv"), path("${repeat_masked_assembly}_gene_stats.tsv"), path("${repeat_masked_assembly}_busco/short_summary.specific..${repeat_masked_assembly}_busco.txt")

   script:
  """
  ml load QUAST/5.0.2-foss-2018b-Python-3.6.6;
  mkdir quast_out;
  quast.py ${repeat_masked_assembly} -o quast_out/;
  cp quast_out/report.tsv ${repeat_masked_assembly}_report.tsv

  exon_count=\$(grep -c "exon" ${gff_file});
  mean_exon_length=\$(awk '\$3 == "exon" {print \$5 - \$4 + 1}' ${gff_file} | awk '{ sum += \$1 } END { if (NR > 0) print sum / NR }');
  median_exon_length=\$(awk '\$3 == "exon" {print \$5 - \$4 + 1}' ${gff_file} | sort -n | awk '{ a[i++] = \$1; } END { if (i % 2 == 0) print (a[i / 2 - 1] + a[i / 2]) / 2; else print a[int(i / 2)]; }');

  intron_count=\$(grep -c "intron" ${gff_file});
  mean_intron_length=\$(awk '\$3 == "intron" {print \$5 - \$4 + 1}' ${gff_file} | awk '{ sum += \$1 } END { if (NR > 0) print sum / NR }');
  median_intron_length=\$(awk '\$3 == "intron" {print \$5 - \$4 + 1}' ${gff_file} | sort -n | awk '{ a[i++] = \$1; } END { if (i % 2 == 0) print (a[i / 2 - 1] + a[i / 2]) / 2; else print a[int(i / 2)]; }');


  cds_length=\$(awk '\$3 == "CDS" {sum += \$5 - \$4 + 1} END {print sum}' ${gff_file});

  echo -e "Exon_count\\t\$exon_count" >> "output_results.tsv";
  echo -e "Mean_exon_length\\t\$mean_exon_length" >> "output_results.tsv";
  echo -e "Median_exon_length\\t\$median_exon_length" >> "output_results.tsv";
  echo -e "Intron_count\\t\$intron_count" >> "output_results.tsv";
  echo -e "Mean_intron_length\\t\$mean_intron_length" >> "output_results.tsv";
  echo -e "Median_intron_length\\t\$median_intron_length" >> "output_results.tsv";
  echo -e "CDS_length\\t\$cds_length" >> "output_results.tsv";

  cp output_results.tsv ${repeat_masked_assembly}_gene_stats.tsv;

  ml purge;
  ml load BUSCO/5.1.2-foss-2018b-Python-3.6.6;

  busco -f -i ${repeat_masked_assembly} -l /scicore/home/schpie00/baer0006/ONT/de_novo_pipeline/assembly/flye/run_2/BUSCO_trial/busco_downloads/lineages/metazoa_odb10/ -o ${repeat_masked_assembly}_busco -m genome --offline
  
  """


}
