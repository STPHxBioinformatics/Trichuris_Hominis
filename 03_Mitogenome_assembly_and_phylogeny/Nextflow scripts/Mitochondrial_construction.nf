nextflow.enable.dsl=2
params.reads = "/scicore/home/schpie00/GROUP/20230601_Raw_sequencing_Trichuris_whole_genome/2023*/GFB-*TRI*_R{1,2}_001.fastq.gz"
params.trimmomatic = "/scicore/home/schpie00/baer0006/IluminiPilot/Nextflow/packages/Trimmomatic-0.39/trimmomatic-0.39.jar"
params.path_getorganelle = "/scicore/home/schpie00/baer0006/Trichuris_Trichiura_CI/Packages/GetOrganelle/get_organelle_from_reads.py"
params.mito_reference = "/scicore/home/schpie00/baer0006/Trichuris_Trichiura_CI/analysis_pipeline/Mitochondrial_references/*.fasta"

log.info """\
	
	=============================================
	SWISS-TPH WHOLE-GENOME-PIPELINE: CAN'T RESIST
	=============================================

	reads directory: 	${params.reads}

	"""
	.stripIndent()

/*
Trim raw reads
*/
process trimming {
	maxForks 100
	
	input:
	tuple val(pair_id), path(reads)
	
	output:
	tuple val(pair_id), path("${pair_id}_trimmed_R{1,2}.fastq")
	
	script:
	"""
	java -jar $baseDir/packages/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33\
		${reads[0]}\
		${reads[1]}\
		${pair_id}_trimmed_R1.fastq\
		${pair_id}_trimmed_R1_unpaired.fastq\
		${pair_id}_trimmed_R2.fastq\
		${pair_id}_trimmed_R2_unpaired.fastq\
		ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50
	"""
}

/*
Use getorganelle to assemble mitochondrial genomes. Only circular mitochondrial genomes which were assembled to completeness are processed further.
*/
process getorganelle {
   	maxForks 60

	input: 
	tuple val(pair_id), path(reads)
	
	output: 
	tuple val(pair_id), path("${pair_id}.fasta"), optional: true
	
	script:
	"""
	ml load Python/3.9.5-GCCcore-10.3.0
	${params.path_getorganelle} \
	-1 ${reads[0]}\
	-2 ${reads[1]}\
	-o ${pair_id}_GetOrganelle\
	-F animal_mt -t 4 --overwrite --reverse-lsc;
	
	if [ -f ${pair_id}_GetOrganelle/animal_mt.K115.complete.graph1.1.path_sequence.fasta ]; then

	cp ${pair_id}_GetOrganelle/animal_mt.K115.complete.graph1.1.path_sequence.fasta ./${pair_id}.fasta

	fi
	"""
}

/*
As a result of the software the reverse complement of the mitochondrial genomes is generated. This little python script flips it around.
*/
process reverse_complement_mito_fasta {
	input:
    tuple val(pair_id), path(mito_fasta)

	output:
	tuple val(pair_id), path("${pair_id}_reversed.fasta")
	
	script:
	"""
	ml load Python/3.9.5-GCCcore-10.3.0-bare
	
	python /scicore/home/schpie00/baer0006/Trichuris_Trichiura_CI/Mitobim_test/rearrange_mito_reverse_complement.py ${mito_fasta}
	"""
}

/*
The mitochondrial genome is annotated using mitoZ to generate a gbf file. There is no mito_fasta output.
*/
process mitoz_assembly {
  maxForks 60

	input:
  tuple val(pair_id), path(mito_fasta)
	
	output:
	tuple val(pair_id), path("${pair_id}_reversed.fasta"), path("MitoZ_${pair_id}_tmp_header.fasta_mitoscaf.fa.gbf")

	script:
	def sampleID = pair_id.substring(30, pair_id.length() - 15).findAll( /\d+/ )*.toInteger()

	"""
	ml load GD/2.66-goolf-1.7.20-Perl-5.22.2;
	ml load SAMtools/1.7-goolf-1.7.20;
	ml load BWA/0.7.17-goolf-1.7.20;
	ml load Perl/5.22.2-goolf-1.7.20;


	sed -e "s/^>.*/>${sampleID} topology=circular/" ${mito_fasta} > ${pair_id}_tmp_header.fasta;

	/scicore/home/schpie00/baer0006/packages/mitoz/MitoZ_v3.6.sif mitoz annotate --clade Nematoda --fastafiles ${pair_id}_tmp_header.fasta --outprefix MitoZ;


	if grep -B 1 'COX1' MitoZ.${pair_id}_tmp_header.fasta.result/MitoZ_${pair_id}_tmp_header.fasta_mitoscaf.fa.gbf | grep -q 'complement'; then \
	
	ml purge; ml load Python/3.9.5-GCCcore-10.3.0;
	
	python /scicore/home/schpie00/baer0006/Trichuris_Trichiura_CI/Mitobim_test/rearrange_mito_reverse_complement.py ${pair_id}_tmp_header.fasta;
	
	ml purge;

	ml load GD/2.66-goolf-1.7.20-Perl-5.22.2;
	ml load SAMtools/1.7-goolf-1.7.20;
	ml load BWA/0.7.17-goolf-1.7.20;
	ml load Perl/5.22.2-goolf-1.7.20;

	/scicore/home/schpie00/baer0006/packages/mitoz/MitoZ_v3.6.sif mitoz annotate --clade Nematoda --fastafiles ${pair_id}_tmp_header_reversed.fasta --outprefix MitoZ;
	cp MitoZ.${pair_id}_tmp_header_reversed.fasta.result/MitoZ_${pair_id}_tmp_header_reversed.fasta_mitoscaf.fa.gbf ./MitoZ_${pair_id}_tmp_header.fasta_mitoscaf.fa.gbf;
	mv ${pair_id}_tmp_header_reversed.fasta ${pair_id}_reversed.fasta;
	else mv ${pair_id}_tmp_header.fasta ${pair_id}_reversed.fasta;
	cp MitoZ.${pair_id}_tmp_header.fasta.result/MitoZ_${pair_id}_tmp_header.fasta_mitoscaf.fa.gbf ./MitoZ_${pair_id}_tmp_header.fasta_mitoscaf.fa.gbf; fi

	"""
}

/*
The starting point of the mitochondrial fasta file is adjusted to the start of the COX1 gene using the GBF file from mitoZ.
*/
process adjust_mito_fasta {
	input:
    	tuple val(pair_id), path(mito_fasta), path(mito_gbf)

	output:
	tuple val(pair_id), path("${pair_id}_reversed_adjusted.fasta")
	
	script:
	"""
	ml load Python/3.9.5-GCCcore-10.3.0-bare
	
	starting_postition="\$(grep -B 1 'COX1' ${mito_gbf} | grep ' gene' | grep -oP '\\d+(?=\\.\\.)' | tail -1)";
	
	python /scicore/home/schpie00/baer0006/Trichuris_Trichiura_CI/Nextflow/rearrange_mito.py ${mito_fasta} "\$(grep -B 1 'COX1' ${mito_gbf} | grep ' gene' | grep -oP '\\d+(?=\\.\\.)' | tail -1)"
	"""
}

/*
The newly adjusted fasta file is assembled again to generate the gbf file with COX1 as a starting position.
*/
process mitoz_assembly_adjusted {
	maxForks 60

	input:
    tuple val(pair_id), path(mito_fasta)
	
	output:
	tuple val(pair_id), path("MitoZ.${pair_id}_tmp_header.fasta.result/MitoZ_${pair_id}_tmp_header.fasta_mitoscaf.fa.gbf"), path("${pair_id}.fasta")

	script:
	def sampleID = pair_id.substring(30, pair_id.length() - 15).findAll( /\d+/ )*.toInteger()

	"""
	ml load GD/2.66-goolf-1.7.20-Perl-5.22.2;
	ml load SAMtools/1.7-goolf-1.7.20;
	ml load BWA/0.7.17-goolf-1.7.20;
	ml load Perl/5.22.2-goolf-1.7.20;
	

	sed -e "s/^>.*/>${sampleID} topology=circular/" ${mito_fasta} > ${pair_id}_tmp_header.fasta;

	/scicore/home/schpie00/baer0006/packages/mitoz/MitoZ_v3.6.sif mitoz annotate --clade Nematoda --fastafiles ${pair_id}_tmp_header.fasta --outprefix MitoZ;
	
	cp ${pair_id}_tmp_header.fasta ${pair_id}.fasta
	"""
}

/*
The mitochondrial genome is visualised using MitoZs visualise option that uses circos and the coverages are also included in the plot.
*/
process mitoz_visualisation {
	maxForks 60
    stageInMode 'copy'
	publishDir '/scicore/home/schpie00/baer0006/Trichuris_Trichiura_CI/mitochondrial_baiting_and_analysis_pipeline/Results_samples/CDS_and_non_CDS/', mode: 'copy'
	
	input:
   	tuple val(pair_id), path(mito_gbf), path(final_fasta), path(reads)
	
	output:
	tuple val(pair_id), path("${pair_id}_mito_circos.png"), path("${final_fasta}"), path("${mito_gbf}")

	script:
	"""
	ml load GD/2.66-goolf-1.7.20-Perl-5.22.2;
	ml load SAMtools/1.7-goolf-1.7.20;
	ml load BWA/0.7.17-goolf-1.7.20;
	ml load Perl/5.22.2-goolf-1.7.20;
   	singularity run /scicore/home/schpie00/baer0006/packages/mitoz/MitoZ_v3.6.sif mitoz visualize \
	--circos /scicore/home/schpie00/baer0006/Trichuris_Trichiura_CI/Nextflow/test_mitoz/circos-0.69-9/bin/circos \
	--gb ${mito_gbf} --outdir ./results --run_map yes \
	--fq1 ${reads[0]} \
	--fq2 ${reads[1]};

	cp ./results/circos.png ./${pair_id}_mito_circos.png
	"""
}



/*
gene_selection
*/
process mafft_new_sequences {
	input:
  path(mito_gbf)

	script:
	"""
	ls | grep "TRICHURISP-" | xargs rm -f;
  ml load Python/3.9.5-GCCcore-10.3.0-bare
  ml load MAFFT/7.490-GCC-10.3.0-with-extensions
  cp /scicore/home/schpie00/baer0006/Trichuris_Trichiura_CI/mitochondrial_baiting_and_analysis_pipeline/results_nucleic_acid_phylogeny/generate_nucleicacid_fastas.py ./
  python generate_nucleicacid_fastas.py;
  for file in *.fasta; do
  sed -i 's/^>/>s/' "\$file";
  sed -i 's/\\[//' "\$file";
  sed -i 's/\\]//' "\$file";
  mafft --auto "\$file" > "\$file"aligned.fasta
  done 


	"""
}










/* 
Extract CDS regions only 
*/
process extract_cds_mito_fasta {
	publishDir '/scicore/home/schpie00/baer0006/Trichuris_Trichiura_CI/mitochondrial_baiting_and_analysis_pipeline/Results_samples/CDS_only/', mode: 'copy'
	
	input:
    	tuple val(pair_id), path(png), path(mito_fasta), path(mito_gbf)

	output:
	tuple val(pair_id), path("${pair_id}_CDS.fasta")
	
	script:
	"""
	ml load Python/3.9.5-GCCcore-10.3.0-bare
	
	python /scicore/home/schpie00/baer0006/Trichuris_Trichiura_CI/Nextflow/extract_coding_sequence.py ${mito_gbf} ${mito_fasta}
	"""
}




/* 
Check if all genes are present visually 
*/
process mitoz_assembly_cds {
	input:
  tuple val(pair_id), path(mito_fasta)
	
	output:
	path("${pair_id}_CDS.png")

	script:

	"""
	ml load GD/2.66-goolf-1.7.20-Perl-5.22.2;
	ml load SAMtools/1.7-goolf-1.7.20;
	ml load BWA/0.7.17-goolf-1.7.20;
	ml load Perl/5.22.2-goolf-1.7.20;
	
	/scicore/home/schpie00/baer0006/packages/mitoz/MitoZ_v3.6.sif mitoz annotate --clade Nematoda --fastafiles ${mito_fasta} --outprefix MitoZ;
	cp MitoZ.${pair_id}_CDS.fasta.result/circos.png ./${pair_id}_CDS.png
	"""
}

/*
The reference is treated in the same way as the samples to avoid bias from the software. In a first step the reference is annotated to generate a gbf file.
*/
process mito_ref_prep {
	maxForks 60

	input:
	path(ref_mito_fasta)
	
	output:
	tuple path("${ref_mito_fasta.baseName.substring(0, 6)}.fasta"), path("MitoZ.${ref_mito_fasta.baseName.substring(0, 6)}.fasta.result/MitoZ_${ref_mito_fasta.baseName.substring(0, 6)}.fasta_mitoscaf.fa.gbf")

	script:
	def refID = ref_mito_fasta.baseName.substring(0, 6)

	"""
	ml load GD/2.66-goolf-1.7.20-Perl-5.22.2;
	ml load SAMtools/1.7-goolf-1.7.20;
	ml load BWA/0.7.17-goolf-1.7.20;
	ml load Perl/5.22.2-goolf-1.7.20;
	

	sed -e "s/^>.*/>${refID} topology=circular/" ${ref_mito_fasta} > ${refID}.fasta;

	/scicore/home/schpie00/baer0006/packages/mitoz/MitoZ_v3.6.sif mitoz annotate --clade Nematoda --fastafiles ${refID}.fasta --outprefix MitoZ;

	"""
}

/*
The starting position of the reference is adjusted to the starting position of the COX1 gene.
*/
process adjust_mito_ref_fasta {
	input:
    tuple path(mito_ref_fasta), path(mito_ref_gbf)

	output:
	path("${mito_ref_fasta.baseName.substring(0, 6)}_adjusted.fasta")
	
	script:
	def refID = mito_ref_fasta.baseName.substring(0, 6)
	"""
	ml load Python/3.9.5-GCCcore-10.3.0-bare
	
	python /scicore/home/schpie00/baer0006/Trichuris_Trichiura_CI/Nextflow/rearrange_mito.py ${mito_ref_fasta} "\$(grep -B 1 'COX1' ${mito_ref_gbf} | grep ' gene' | grep -oP '\\d+(?=\\.\\.)' | head -1)"
	"""
}

/*
The reference is re-annotated with the adjusted starting position and the results put into the publishDir.
*/
process mitoz_ref_assembly_adjusted {
	maxForks 60

	publishDir '/scicore/home/schpie00/baer0006/Trichuris_Trichiura_CI/mitochondrial_baiting_and_analysis_pipeline/Results_references/', mode: 'copy'
	
	input:
    path(mito_ref_fasta)
	
	output:
	tuple path("${mito_ref_fasta}"), path("MitoZ.${mito_ref_fasta}.result/MitoZ_${mito_ref_fasta}_mitoscaf.fa.gbf"), path("${mito_ref_fasta.baseName.substring(0, 6)}.png")

	script:
	def refID = mito_ref_fasta.baseName.substring(0, 6)

	"""
	ml load GD/2.66-goolf-1.7.20-Perl-5.22.2;
	ml load SAMtools/1.7-goolf-1.7.20;
	ml load BWA/0.7.17-goolf-1.7.20;
	ml load Perl/5.22.2-goolf-1.7.20;
	
	/scicore/home/schpie00/baer0006/packages/mitoz/MitoZ_v3.6.sif mitoz annotate --clade Nematoda --fastafiles ${mito_ref_fasta} --outprefix MitoZ;
	cp MitoZ.${mito_ref_fasta}.result/circos.png ./${refID}.png
	"""
}

/*
Extract CDS of reference. The final file still has to be checked and manually curated. Had some issues with the first numbers, didnt sprinkle on much syntax sugar yet
*/
process extract_cds_ref_mito_fasta {
	publishDir '/scicore/home/schpie00/baer0006/Trichuris_Trichiura_CI/mitochondrial_baiting_and_analysis_pipeline/Results_references/CDS_fastas/', mode: 'copy'
	
	input:
    	tuple path(mito_ref_fasta), path(mito_ref_gbf), path(png_dummy)

	output:
	path("${mito_ref_fasta.baseName.substring(0, 6)}_adjusted_CDS.fasta")
	
	script:
	def refID = mito_ref_fasta.baseName.substring(0, 6)

	"""
	ml load Python/3.9.5-GCCcore-10.3.0-bare
	
	python /scicore/home/schpie00/baer0006/Trichuris_Trichiura_CI/Nextflow/extract_coding_sequence.py ${mito_ref_gbf} ${mito_ref_fasta}
	"""
}

workflow {
	read_pairs_ch = channel.fromFilePairs(params.reads, checkIfExists: true)
	reference_mito_ch = channel.fromPath(params.mito_reference)

	getorganelle(read_pairs_ch)
	assembled_mito_ch = getorganelle.out
	
	reverse_complement_mito_fasta(assembled_mito_ch)
	reverse_complemented_mito_ch = reverse_complement_mito_fasta.out
	
	mitoz_assembly(reverse_complemented_mito_ch)
	mito_gbf_file_ch = mitoz_assembly.out

	adjust_mito_fasta(mito_gbf_file_ch)
	adjusted_mito_ch = adjust_mito_fasta.out
	
	mitoz_assembly_adjusted(adjusted_mito_ch)
	mito_adjusted_gbf_file_ch = mitoz_assembly_adjusted.out	
	
	mitoz_visualisation_ch = mito_adjusted_gbf_file_ch.join(read_pairs_ch)
	
	mitoz_visualisation(mitoz_visualisation_ch)
	mitoz_visualitation_new_ch = mitoz_visualisation.out

	mafft_new_sequences(mitoz_visualitation_new_ch.map { it[3] }.collect())

	extract_cds_mito_fasta(mitoz_visualitation_new_ch)
	cds_adjusted_mito_ch = extract_cds_mito_fasta.out

	mitoz_assembly_cds(cds_adjusted_mito_ch)

	mito_ref_prep(reference_mito_ch)
	mito_ref_gbf_ch = mito_ref_prep.out

	adjust_mito_ref_fasta(mito_ref_gbf_ch)
	adjusted_mito_ref_ch = adjust_mito_ref_fasta.out

	mitoz_ref_assembly_adjusted(adjusted_mito_ref_ch)
	mitoz_ref_assembly_adjusted_ch = mitoz_ref_assembly_adjusted.out

	extract_cds_ref_mito_fasta(mitoz_ref_assembly_adjusted_ch)
	adjusted_mito_ref_cds_ch = extract_cds_ref_mito_fasta.out
}