# Genome wide association study 

Author: Max Bär, max.baer[at]swisstph.ch

## Contents:

This section presents the code used for the genome wide association study of drug sensitive and non-sensitive *T. hominibus* populations.
 It is presented by the individual nextflow processes which are combined in a final workflow or individual sections for further processing outside
of Nextflow. The scripts were run on the SciCORE computing cluster at the University of Basel with most modules pre-installed.

Table of contents

1. [Variant calling pipeline](#variants)
   1. [Nextflow settings](#settings)
   2. [Indexing assembled reference - process index](#index)
   3. [Aligning to reference genome - process alignment](#align)
   4. [Sorting and filtering bam file - process sorting](#sort)
   5. [QC flagstat, Kraken - process flagstat](#flagstat)
   6. [Mark duplicates - process mark_duplicates_spark](#duplicates)
   7. [Classify gender using coverage - process classify_gender](#gender)
   8. [Initial basecalling for BQSR - process initial_basecalling_BQSR](#initial_BQSR)
   9. [Initial hard filtering for BQSR - process initial_SNP_BQSR](#initial)
   10. [BQSR - process BQSR](#BQSR)
   11. [GVCF generation and merging - process GVCF](#GVCF)
   12. [Genotyping and filtering](#geno)
2. [Nextflow workflow](#NF1)
3. [GWAS in PLINK](#plink)
___
## Variant calling pipeline <a name="variants"></a>
GATK was used to call variants.

### Nextflow environment variables <a name="settings"></a>
DSL 2 option was used with nextflow and directories of genomes and databases used were defined outside of processes.

      nextflow.enable.dsl=2
      params.reads = "/*/GROUP/20230601_Raw_sequencing_Trichuris_whole_genome/2023*/GFB-*TRI*_R{1,2}_001.fastq.gz"
      params.outdir = "/*/Trichuris_Hominibus_CI/Nextflow/"
      params.reference = "/*/analysis_pipeline/Reference/Trichuris_cote_divoire_freeze.fasta"
      params.reference_path = "/*/analysis_pipeline/Reference/Trichuris_cote_divoire_freeze"
      params.path_kraken2_db = "/*/data/managed/.store/kraken_191029/"
      params.flagstat_output = "/*/analysis_pipeline/flagstat/"


### Indexing assembled reference - process index <a name="index"></a>
BWA, version 0.7.17, and SAMtools, version 1.14, was used to index the assembled reference genome.
      
      process index {
          
          input:
          path(reference)
      
          output:
          val(reference.baseName)
      
          script:
          """
          ml load BWA/0.7.17-GCC-10.3.0
          ml load SAMtools/1.14-GCC-10.3.0
          bwa index ${reference};
          
          samtools faidx ${reference};
          samtools dict ${reference} > ${params.reference_path}.dict
          """
      }

### Aligning to reference genome - process alignment <a name="align"></a>
BWA, version 0.7.17, and SAMtools, version 1.14, was used to align trimmed illumina reads to the assembled reference genome.

      process alignment {
          maxForks 60
      
          input:
          val(reference)
          tuple val(pair_id), path(reads)
          
          output:
          tuple val(pair_id), path("${pair_id}.bam")
          
          script:
          readGroup = "@RG\\tID:${pair_id}\\tLB:${pair_id}\\tPL:ILLUMINA\\tPM:NOVASEQ\\tSM:${pair_id}"
          """
          ml load BWA/0.7.17-GCC-10.3.0
          ml load SAMtools/1.14-GCC-10.3.0
          bwa mem \
          -t 4 $params.reference \
          -R \"${readGroup}\" \
          ${reads[0]} \
          ${reads[1]} | samtools view --threads 4 -b - \
          > ${pair_id}.bam
          """
      }

### Sorting and filtering bam file - process sorting <a name="sort"></a>
Bam file was sorted using SAMtools, version 1.14.

      process sorting {
          maxForks 20
      
          input:
          tuple val(pair_id), path(reads)
      
          output:
          tuple val(pair_id), path("${pair_id}_filtered.bam")
      
          script:
          """
          ml load SAMtools/1.14-GCC-10.3.0
          samtools sort -o ${pair_id}_sorted.bam $reads;
          samtools view -F 4 -b -o ${pair_id}_filtered.bam ${pair_id}_sorted.bam;
          samtools index ${pair_id}_filtered.bam ${pair_id}_filtered.bam.bai
          """
      }

### QC flagstat, Kraken - process flagstat <a name="flagstat"></a>
SAMtools version 1.14, MultiQC version 1.11, and Kraken2 version 2.1.1 was used to generate QC files. QC files are located in
the QC files folder.

      process flagstat {
          maxForks 20
          
          input:
          tuple val(pair_id), path(reads)
      
          output: 
          file("${pair_id}.flagstat")
      
          script:
          """
          ml load SAMtools/1.14-GCC-10.3.0
          samtools flagstat ${reads} > ${pair_id}.flagstat
          """
      }
      
      process multiqc_flagstat {
          publishDir "$params.flagstat_outputrun2", mode: 'copy'
      
          input:
          file(flagstat)
      
          output:
          file "multiqc_report_flagstat.html"
          
          script:
          """
          ml load MultiQC/1.11-foss-2018b-Python-3.6.6
          multiqc .;
          mv multiqc_report.html multiqc_report_flagstat.html
          """
      }
      
      process kraken2_impurities {
          maxForks 2
          
          input:
          tuple val(pair_id), path(reads)
          
          output:
          file("${pair_id}_kraken2report")
          
          script:
          """
          ml load Kraken2/2.1.1-foss-2018b-Perl-5.28.0
          kraken2 \
          --db $params.path_kraken2_db \
          --report ${pair_id}_kraken2report \
          --paired ${reads[0]} ${reads[1]}
          """
      }
      
      process multiqc_kraken2 {
          publishDir "$params.flagstat_output*", mode: 'copy'
      
          input:
          file(kraken2)
      
          output:
          file "multiqc_report_kraken2.html"
          
          script:
          """
          ml load MultiQC/1.11-foss-2018b-Python-3.6.6
          multiqc .;
          mv multiqc_report.html multiqc_report_kraken2.html
          """
      }

### Mark duplicates - process mark_duplicates_spark <a name="duplicates"></a>
Duplicated reads are marked for basecalling.

      process mark_duplicates_spark {
          maxForks 5
          
          input:
          tuple val(pair_id), path(reads)
      
          output:
          tuple val(pair_id), path("${pair_id}_sorted_dedup.bam")
          
          script:
          """
          ml purge
          ml load Java/8u212b03
          ml load GATK/4.2.6.1-foss-2018b-Java-1.8
          gatk MarkDuplicatesSpark \
          -I ${reads} \
          -M ${pair_id}_dedup_metrics.txt \
          -O ${pair_id}_sorted_dedup.bam 
          """
      }


### Classify gender using coverage - process classify_gender <a name="gender"></a>
Gender was classified using the coverage ratio of the largest X-chromosome scaffold. To the largest Scaffold on Chromosome 2.
Gender data is provided in the QC folder.

      process classify_gender {
          input:
          tuple val(pair_id), path(bam_file)
      
          output:
          path("${pair_id}.txt")
      
          script:
          """
          ml load SAMtools/1.16.1-GCC-10.3.0
          samtools index ${bam_file};
          cov_1="\$(samtools coverage -r TTRE_chr1_scaffold1 ${bam_file} | awk 'NR==2 {print \$7}')";
          cov_2="\$(samtools coverage -r TTRE_chr2_scaffold39 ${bam_file} | awk 'NR==2 {print \$7}')";
          ratio=\$(echo "\$cov_1 / \$cov_2" | bc -l)
          echo "\$cov_1\n\$cov_2\n\$ratio" >> ${pair_id}.txt
          """
      
      }
      
      process merge_gender_data {
          input:
          path(gender_files)
      
          output:
          path("merged_coverage.txt")
      
          script:
          """
          ml load Python/3.9.5-GCCcore-10.3.0
          python /scicore/home/schpie00/baer0006/Trichuris_Trichiura_CI/analysis_pipeline/merge_coverage.py
          """
      
      }



### Initial basecalling for BQSR - process initial_basecalling_BQSR <a name="initial_BQSR"></a>
Basecalling was done in a first round to obtain a subset of very likely SNP's for base quality score recalibration.
Plots were generated for QUAL, DP, QD, FS, MQ, MQRankSum, SQR, ReadPosRankSum. This code is largely adapted from the 
Doyle et al. paper on ancient and modern Trichuris
![Plot_Stats](plot_nuclear_variant_summaries.png)
![Plot_Stats](table_nuclear_variant_quantiles.png)

      process initial_basecalling_BQSR {
          maxForks 56
          
          input:
          tuple val(pair_id), path(reads)
      
          output:
          tuple val(pair_id), path("${pair_id}_first_raw_variants.vcf")
          
          script:
          """
          ml load  GATK/4.2.4.1-GCCcore-10.3.0-Java-11 
          gatk HaplotypeCaller \
          -R ${params.reference} \
          -I ${reads} \
          --heterozygosity 0.015 \
              --indel-heterozygosity 0.01 \
          -O ${pair_id}_first_raw_variants.vcf
          """
      }


### Extract coding sequences - process extract_cds_mito_fasta <a name="MitoZ6"></a>
A custom python script was used to generally extract only coding sequences, not by gene.

      process extract_cds_mito_fasta {
          publishDir '/*/Results_samples/CDS_only/', mode: 'copy'
          
          input:
              tuple val(pair_id), path(png), path(mito_fasta), path(mito_gbf)
      
          output:
          tuple val(pair_id), path("${pair_id}_CDS.fasta")
          
          script:
          """
          ml load Python/3.9.5-GCCcore-10.3.0-bare
          
          python /*/extract_coding_sequence.py ${mito_gbf} ${mito_fasta}
          """
      }

python script

      import argparse
      from Bio import SeqIO
      
      # Parse GBF file and extract coding regions
      def extract_cds(gbf_file, fasta_file):
          coding_regions = []
          with open(gbf_file, "r") as gbf:
              for line in gbf:
                  if "CDS" in line:
                      parts = line.split()
                      parts_new = parts[1].strip('complement()<').split(sep='..')
                      start = int(parts_new[0])
                      end = str(parts_new[1])
                      end2 = int(''.join(c for c in end if c.isdigit()))
                      print(start, end2)
                      coding_regions.append((start, end2))
      
          #Get header
          with open(fasta_file, "r") as file:
              lines = file.readlines()
              header = lines[0].strip()
      
      
          # Read mitochondrial sequence FASTA file and extract coding regions
          coding_sequence = ""
          for record in SeqIO.parse(fasta_file, "fasta"):
              sequence = str(record.seq)
              for start, end in coding_regions:
                  coding_sequence += sequence[start-1:end]
      
          # Create a new FASTA file with coding regions
          adjusted_file = fasta_file.replace(".fasta", "_CDS.fasta")
          with open(adjusted_file, "w") as file:
              file.write(f"{header}\n{coding_sequence}")
      
      parser = argparse.ArgumentParser()
      parser.add_argument("gbf_file", help="Path to the gbf file")
      parser.add_argument("fasta_file", help="Path to the fasta_file file")
      args = parser.parse_args()
      
      cds_adjusted_fasta_file = extract_cds(args.gbf_file, args.fasta_file)
      print("Adjusted FASTA file:", cds_adjusted_fasta_file)

### Annotation and visualization of coding sequences as a QC - process mitoz_assembly_cds <a name="MitoZ7"></a>
The concatenated coding sequences from the mitogenome were visualized again using MitoZ

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

### Reference mitochondrial genome annotation - process adjust_mito_ref_fasta, mitoz_ref_assembly_adjusted, extract_cds_ref_mito_fasta <a name="MitoZ8"></a>
Slightly repetitive to the annotation, starting point adjustment and visualisation of the assembled mitogenomes is the 
annotation of the reference genomes. 

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
      
          /*/mitoz/MitoZ_v3.6.sif mitoz annotate --clade Nematoda --fastafiles ${refID}.fasta --outprefix MitoZ;
      
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
          
          python /*/rearrange_mito.py ${mito_ref_fasta} "\$(grep -B 1 'COX1' ${mito_ref_gbf} | grep ' gene' | grep -oP '\\d+(?=\\.\\.)' | head -1)"
          """
      }
      
      /*
      The reference is re-annotated with the adjusted starting position and the results put into the publishDir.
      */

      process mitoz_ref_assembly_adjusted {
          maxForks 60
      
          publishDir '/*/Results_references/', mode: 'copy'
          
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
          
          /*/mitoz/MitoZ_v3.6.sif mitoz annotate --clade Nematoda --fastafiles ${mito_ref_fasta} --outprefix MitoZ;
          cp MitoZ.${mito_ref_fasta}.result/circos.png ./${refID}.png
          """
      }
      
      /*
      Extract CDS of reference. The final file still has to be checked and manually curated. Had some issues with the first numbers, didnt sprinkle on much syntax sugar yet
      */

      process extract_cds_ref_mito_fasta {
          publishDir '/*/Results_references/CDS_fastas/', mode: 'copy'
          
          input:
          tuple path(mito_ref_fasta), path(mito_ref_gbf), path(png_dummy)
      
          output:
          path("${mito_ref_fasta.baseName.substring(0, 6)}_adjusted_CDS.fasta")
          
          script:
          def refID = mito_ref_fasta.baseName.substring(0, 6)
      
          """
          ml load Python/3.9.5-GCCcore-10.3.0-bare
          
          python /*/extract_coding_sequence.py ${mito_ref_gbf} ${mito_ref_fasta}
          """
      }

## Nextflow workflow <a name="NF1"></a>
All processes were stitched together in one workflow

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

## Phylogenetic inference <a name="pyhlo"></a>
BEAST2 (version 2.7.5) was used for all phylogenetic inference with a pure birth process (Yule process). 
Aligned concatenated nucleic acid sequences from all assembled and reference mitochondrial genomes were 
treated as homochromous. The JC69 substitution model with gamma-distributed rate heterogeneity was used 
(JC69 + Γ4) and a strict molecular clock was assumed. Trees and clock models were linked and all model 
parameters were estimated jointly. A Markov Chain Monte Carlo was run for each analysis. Tracer 
(version 1.7.2) was used to assess convergence and effective sample size (should be over 200).
The percentage of samples discarded as burn-in was at least 10%. For the amino acid inferred species tree,
2 assembled sequences of each subclade clade with a posterior distribution over 0.95 were selected at random.
The mito REV substitution model was used and all genes were partitioned separately, allowing for different 
rates for each gene. Effective sample size (ESS) was at least 900 for each of the inferred parameters. The 
intraspecies phylogeny was inferred using a TN96 + Γ4 model. Each gene was partitioned separately and codon 
positions 1 and 2 were partitioned separately to 3. The mean mutation rates of the genes and there codon positions 
are provided in the supporting information and show a faster mutation rate for codon position 3 in all cases.

Snippet to run BEAST2:

      ml load beagle-lib/3.1.2-GCCcore-10.3.0  
      /*/beast/bin/beast -overwrite -threads 32 concatenated_JC69_model_for_rate.xml

Snippet to run Treeannotator, using an older BEAST version:

      ml load Beast/2.6.6-foss-2018b
      treeannotator -lowMem concatenated_JC69_model_for_rate-aligned_merged_mitochondria.trees concatenated_JC69_model_for_rate.tree

XML and log-files are provided in corresponding folders.

