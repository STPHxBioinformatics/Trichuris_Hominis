nextflow.enable.dsl=2

params.illumina_reads = "/scicore/home/schpie00/GROUP/20230601_Raw_sequencing_Trichuris_whole_genome/20230530165811198-27395/GFB-7577_HWNGYDSX5_3_TRICHURIS777Tiagba_S147_L003_R{1,2}_001.fastq.gz"
params.nanopore_reads_dir = "/scicore/home/schpie00/baer0006/ONT/de_novo_pipeline/fastq_files/duplex_called_reads/"
params.reference_whole_genome_doyle = "/scicore/home/schpie00/baer0006/ONT/de_novo_pipeline/assembly/results/trichuris_trichiura.fa"

log.info """\

	=============================================
	SWISS-TPH WHOLE-GENOME-PIPELINE: CAN'T RESIST
	=============================================

	Babedibubbi

	"""
	.stripIndent()

workflow assembly {
    chopper()
    nanopore_reads_ch = chopper.out

    kraken2_nanopore(nanopore_reads_ch)
    nanopore_reads_cleaned_ch = kraken2_nanopore.out

    illumina_reads_ch = channel.fromFilePairs(params.illumina_reads, checkIfExists: true)
    kraken2_illumina(illumina_reads_ch)
    illumina_reads_cleaned_ch = kraken2_illumina.out

    masurca_run(illumina_reads_cleaned_ch, nanopore_reads_cleaned_ch)
    raw_masurca_assembly_ch = masurca_run.out
    
    linkage_group_MaSuRCA(raw_masurca_assembly_ch)

    ragtag_correct(masurca_run.out)
    ragtag_corrected_ch = ragtag_correct.out

    ragtag_scaffold(ragtag_corrected_ch)
    ragtag_scaffolded_ch = ragtag_scaffold.out

    linkage_group_ragtag(ragtag_scaffolded_ch)
    final_genome_ch = linkage_group_ragtag.out

    quality_control(final_genome_ch)
}

workflow {
  assembly()
}

/*
Run chopper on nanopore reads - done
*/
process chopper {
	output:
	path("nanopore_reads.fastq")
	
	script:
	"""
	ml load Python/3.9.5-GCCcore-10.3.0
  cat ${params.nanopore_reads_dir}barcode* | /scicore/home/schpie00/baer0006/packages/chopper-0.6.0/target/release/chopper -l 5000 -q 10 --headcrop 34 --maxlength 60000 --threads 8 > nanopore_reads.fastq
	"""
}

/*
Kraken2 cleaning of Illumina reads - done
*/
process kraken2_illumina {
	input:
	tuple val(pair_id), path(reads)

	output:
	tuple val(pair_id), path("${pair_id}uc_{1,2}.fastq")

	script:
	"""
	ml load Kraken2/2.1.1-foss-2018b-Perl-5.28.0
	kraken2 \
	--db /scicore/data/managed/.store/kraken_191029/ \
	--threads 16 \
	--report ${pair_id}_kraken2report \
	--paired \
	--unclassified-out ${pair_id}uc#.fastq \
	--output kraken_out \
	${reads[0]} ${reads[1]}
	"""
}

/*
Kraken2 cleaning of nanopore reads, basically find out percentage of reads then filter every read with >25% contamination. >25% seems fair. Use python to determine the contamination. Use R to execute the assignment step. Encorporate system variables.
#sed -n '1~4s/^@/>/p;2~4p' nanopore_reads_c.fastq > nanopore_reads_c.fasta ;
*/
process kraken2_nanopore {
	input:
	path(nanopore_reads)

	output: 
	path("clean_nanopore_reads.fastq")

	script:
	"""
	ml load Kraken2/2.1.1-foss-2018b-Perl-5.28.0;
	kraken2 \
	--db /scicore/data/managed/.store/kraken_191029/ \
    --threads 16 \
	--report nanopore_kraken2report \
	--classified-out nanopore_reads_c.fastq \
	--unclassified-out nanopore_reads_uc.fastq \
	--output kraken_out \
	${nanopore_reads} ;

	ml purge;
	ml load Python/3.9.5-GCCcore-10.3.0;
	python /scicore/home/schpie00/baer0006/ONT/de_novo_pipeline/ratio_of_contaminant_from_kraken.py kraken_out > contamination_ratio.txt ;


    ml purge
    ml load R/4.3.0-foss-2021a
    Rscript /scicore/home/schpie00/baer0006/ONT/de_novo_pipeline/contamination_ratio_commandline.R "contamination_ratio.txt" ;

    awk 'NR==FNR { headers[\$0]; next } FNR%4 == 1 { header = \$0 } header in headers { next } 1' filtered_data.tsv ${nanopore_reads} > clean_nanopore_reads.fastq
	"""
}

process masurca_run {
	input:
	tuple val(pair_id), path(illumina_reads)
	path(nanopore_reads)

	output:
	path("assembly.fasta")

    script:
    """
    ml load MaSuRCA/4.0.9-foss-2021a-Perl-5.32.1;
    masurca -t 32 \
    -i ${illumina_reads[0]},${illumina_reads[1]} \
    -r ${nanopore_reads};

    cp CA.mr.99.17.15.0.02/primary.genome.scf.fasta ./assembly.fasta
    """
}


process linkage_group_MaSuRCA {


  input:
	path(raw_masurca_assembly)

	output:
	path("genome_with_new_headers_sorted.fasta")

  script:
    """
    ml load minimap2/2.20-GCCcore-10.3.0;
    ml load Python/3.9.5-GCCcore-10.3.0;
    source $HOME/venv_liftoff/bin/activate;

    liftoff -db /scicore/home/schpie00/baer0006/ONT/de_novo_pipeline/assembly/results/trichuris_muris.PRJEB126.WBPS17.annotations.gff3_db \
    ${raw_masurca_assembly} /scicore/home/schpie00/baer0006/ONT/de_novo_pipeline/assembly/results/trichuris_muris.PRJEB126.WBPS17.genomic.fa \
    -o Masurca_liftoff.gff \
    -u Masurca_unmapped_liftoff;

    grep ">" ${raw_masurca_assembly} > raw_masurca_assembly_headers

    ml purge;
    ml load R/4.3.0-foss-2021a;
    Rscript /scicore/home/schpie00/baer0006/ONT/de_novo_pipeline/Assign_LG_to_scaffolds_clean.R Masurca_liftoff.gff \
    /scicore/home/schpie00/baer0006/ONT/de_novo_pipeline/trichuris_muris.PRJEB126.WBPS17.annotations.onlygenes.gff \
    raw_masurca_assembly_headers;

    awk -F'\t' 'FNR==NR {f2[\$2]=\$1;next} \$2 in f2 {\$2=f2[\$2]}1' New_Old_names3 FS='>' OFS='>' ${raw_masurca_assembly} > genome_with_new_headers.fasta

    ml purge
    conda init bash;
    eval "\$(conda shell.bash hook)";
    conda activate seqkit;
    seqkit sort -nN genome_with_new_headers.fasta -o genome_with_new_headers_sorted.fasta
    """

}


process ragtag_correct {

    input:
	path(raw_masurca_assembly)

	output:
	path("ragtag_corrected.fasta")


    script:
    """
    ml load minimap2/2.20-GCCcore-10.3.0;
    ml load Python/3.9.5-GCCcore-10.3.0;
    source $HOME/venv_ragtag/bin/activate;
    ragtag.py correct \
    ${params.reference_whole_genome_doyle} \
    ${raw_masurca_assembly} \
    -o ragtag_corrected;

    cp ragtag_corrected/ragtag.correct.fasta ragtag_corrected.fasta
    """

}

process ragtag_scaffold {
    input:
	path(masurca_corrected_assembly)

	output:
	path("ragtag_scaffolded2.fasta")


    script:
    """
    ml load minimap2/2.20-GCCcore-10.3.0;
    ml load Python/3.9.5-GCCcore-10.3.0;
    source $HOME/venv_ragtag/bin/activate;
    ragtag.py scaffold -u \
    ${params.reference_whole_genome_doyle} \
    ${masurca_corrected_assembly} \
    -o ragtag_scaffolded;

    cp ragtag_scaffolded/ragtag.scaffold.fasta ragtag_scaffolded.fasta;

    ragtag.py scaffold -u \
    ${params.reference_whole_genome_doyle} \
    ragtag_scaffolded.fasta \
    -o ragtag_scaffolded2;

    cp ragtag_scaffolded2/ragtag.scaffold.fasta ragtag_scaffolded2.fasta
    """

}

process linkage_group_ragtag {

  input:
	path(masurca_assembly_final)

	output:
	path("Trichuris_cote_divoire_freeze.fasta")

    script:
    """
    ml load minimap2/2.20-GCCcore-10.3.0;
    ml load Python/3.9.5-GCCcore-10.3.0;
    source $HOME/venv_liftoff/bin/activate;

    liftoff -db /scicore/home/schpie00/baer0006/ONT/de_novo_pipeline/assembly/results/trichuris_muris.PRJEB126.WBPS17.annotations.gff3_db \
    ${masurca_assembly_final} /scicore/home/schpie00/baer0006/ONT/de_novo_pipeline/assembly/results/trichuris_muris.PRJEB126.WBPS17.genomic.fa \
    -o Masurca_liftoff.gff \
    -u Masurca_unmapped_liftoff;

    grep ">" ${masurca_assembly_final} > raw_masurca_assembly_headers

    ml purge;
    ml load R/4.3.0-foss-2021a;
    Rscript /scicore/home/schpie00/baer0006/ONT/de_novo_pipeline/Assign_LG_to_scaffolds_clean.R Masurca_liftoff.gff \
    /scicore/home/schpie00/baer0006/ONT/de_novo_pipeline/trichuris_muris.PRJEB126.WBPS17.annotations.onlygenes.gff \
    raw_masurca_assembly_headers;

    awk -F'\t' 'FNR==NR {f2[\$2]=\$1;next} \$2 in f2 {\$2=f2[\$2]}1' New_Old_names3 FS='>' OFS='>' ${masurca_assembly_final} > genome_with_new_headers.fasta

    ml purge
    conda init bash;
    eval "\$(conda shell.bash hook)";
    conda activate seqkit;
    seqkit sort -nN genome_with_new_headers.fasta -o Trichuris_cote_divoire_freeze.fasta
    """
}


process quality_control {
  input:
	path(final_assembly)

	output:
	path("busco_out_1")


    script:
    """
    ml purge;
    ml load BUSCO/5.1.2-foss-2018b-Python-3.6.6
    busco -i ${final_assembly} \
    -l /scicore/home/schpie00/baer0006/ONT/de_novo_pipeline/assembly/flye/run_2/BUSCO_trial/busco_downloads/lineages/metazoa_odb10/ \
    -o busco_out \
    -m genome \
    -c 32 \
    --offline;
    cp busco_out/short_summary.specific..busco_out.txt busco_out_1

    ml purge;
    ml load QUAST/5.0.2-foss-2018b-Python-3.6.6;
    quast ${final_assembly}
    """
}

