# Mitochondrial assembly and phylogenetic inference

Author: Max Bär, max.baer[at]swisstph.ch

## Contents:

This section presents the code used for the assembly of *T. incognita* mitogenomes and phylogenetic inference. It is presented by
the individual nextflow processes which are combined in a final workflow or individual sections for further processing outside
of Nextflow. The scripts were run on the SciCORE computing cluster at the University of Basel with most modules pre-installed.

Table of contents

1. [Mitochondrial assembly and annotation](#annotation)
   1. [Nextflow settings](#settings)
   2. [Mitochondrial assembly - process getorganelle](#assembly)
   3. [Identify and flip reverse complement sequences - process reverse_complement_mito_fasta](#flip)
   4. [Mitochondrial annotation using MitoZ - process mitoz_assembly](#MitoZ1)
   5. [Startingpoint of mitochondrial genome is set to the COX1 gene - process adjust_mito_fasta](#MitoZ2)
   6. [Re-annotation of newly adjusted starting position - process mitoz_assembly_adjusted](#MitoZ3)
   7. [Visualization of mitochondrial genomes using MitoZ - process mitoz_visualisation](#MitoZ4)
   8. [Extraction of gene sequences and MAFFT alignment of individual genes - process mafft_new_sequences](#Mafft)
   9. [Extract coding sequences - process extract_cds_mito_fasta](#MitoZ6)
   10. [Annotation and visualization of coding sequences as a QC - process mitoz_assembly_cds](#MitoZ7)
   11. [Reference mitochondrial genome annotation - process adjust_mito_ref_fasta, mitoz_ref_assembly_adjusted, extract_cds_ref_mito_fasta](#MitoZ8)
2. [Nextflow workflow](#NF1)
3. [Phylogenetic inference](#pyhlo)
___
## Mitochondrial assembly and annotation <a name="annotation"></a>
GetOrganelle was used to bait and reconstruct mitochondrial sequences from whole genome data.

### Nextflow environment variables <a name="settings"></a>
DSL 2 option was used with nextflow and directories of genomes and databases used were defined outside of processes.

      nextflow.enable.dsl=2
      params.reads = "/*/2023*/GFB-*TRI*_R{1,2}_001.fastq.gz"
      params.path_getorganelle = "/*/GetOrganelle/get_organelle_from_reads.py"
      params.mito_reference = "/*/Mitochondrial_references/*.fasta"

### Mitochondrial assembly - process getorganelle <a name="assembly"></a>
Getorganelle, version 1.7.7.0 was used to assemble mitochondrial genomes. Only circular mitochondrial genomes which were 
assembled to completeness are processed further.
      
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

### Identify and flip reverse complement sequences - process reverse_complement_mito_fasta <a name="flip"></a>
As a result of the software the reverse complement of the mitochondrial genomes is generated. A custom python script 
identifies reverse complement sequences and flips them around.

      process reverse_complement_mito_fasta {
          input:
          tuple val(pair_id), path(mito_fasta)
      
          output:
          tuple val(pair_id), path("${pair_id}_reversed.fasta")
          
          script:
          """
          ml load Python/3.9.5-GCCcore-10.3.0-bare
          
          python /*/Mitobim_test/rearrange_mito_reverse_complement.py ${mito_fasta}
          """
      }

Python script

      import argparse
      from Bio.Seq import Seq
      
      
      def modify_fasta_file(fasta_file):
          with open(fasta_file, "r") as file:
              lines = file.readlines()
      
          header = lines[0].strip()
          seq = Seq("".join(lines[1:]).replace("\n", ""))
      
          # Reverse complement sequence
          reversed_complement_sequence = str(seq.reverse_complement())
      
          # Write the adjusted sequence to a new file
          reversed_file = fasta_file.replace(".fasta", "_reversed.fasta")
          with open(reversed_file, "w") as file:
              file.write(f"{header}\n{reversed_complement_sequence}")
      
          return reversed_file
      
      # Provide the path to your GetOrganelle-generated mitochondrial genome in FASTA format
      
      parser = argparse.ArgumentParser()
      parser.add_argument("fasta_file", help="Path to the FASTA file")
      args = parser.parse_args()
      
      reversed_fasta_file = modify_fasta_file(args.fasta_file)
      print("Adjusted FASTA file:", reversed_fasta_file)

### Mitochondrial annotation using MitoZ - process mitoz_assembly <a name="MitoZ1"></a>
MitoZ, version 3.6, was used to annotate the mitochondrial genomes. In cases where the reversed complement was present the
mitochondrial genome was flipped again and re-annotated.

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
      
         /*/mitoz/MitoZ_v3.6.sif mitoz annotate --clade Nematoda --fastafiles ${pair_id}_tmp_header.fasta --outprefix MitoZ;
      
      
         if grep -B 1 'COX1' MitoZ.${pair_id}_tmp_header.fasta.result/MitoZ_${pair_id}_tmp_header.fasta_mitoscaf.fa.gbf | grep -q 'complement'; then \
          
         ml purge; ml load Python/3.9.5-GCCcore-10.3.0;
          
         python /*/rearrange_mito_reverse_complement.py ${pair_id}_tmp_header.fasta;
          
         ml purge;
      
         ml load GD/2.66-goolf-1.7.20-Perl-5.22.2;
         ml load SAMtools/1.7-goolf-1.7.20;
         ml load BWA/0.7.17-goolf-1.7.20;
         ml load Perl/5.22.2-goolf-1.7.20;
      
         /*/mitoz/MitoZ_v3.6.sif mitoz annotate --clade Nematoda --fastafiles ${pair_id}_tmp_header_reversed.fasta --outprefix MitoZ;
         cp MitoZ.${pair_id}_tmp_header_reversed.fasta.result/MitoZ_${pair_id}_tmp_header_reversed.fasta_mitoscaf.fa.gbf ./MitoZ_${pair_id}_tmp_header.fasta_mitoscaf.fa.gbf;
         mv ${pair_id}_tmp_header_reversed.fasta ${pair_id}_reversed.fasta;
         else mv ${pair_id}_tmp_header.fasta ${pair_id}_reversed.fasta;
         cp MitoZ.${pair_id}_tmp_header.fasta.result/MitoZ_${pair_id}_tmp_header.fasta_mitoscaf.fa.gbf ./MitoZ_${pair_id}_tmp_header.fasta_mitoscaf.fa.gbf; fi
         """
      }

### Startingpoint of mitochondrial genome is set to the COX1 gene - process adjust_mito_fasta <a name="MitoZ2"></a>
The starting point of the mitochondrial genome was set to the COX1 gene using a custom python script.

      process adjust_mito_fasta {
         input:
         tuple val(pair_id), path(mito_fasta), path(mito_gbf)
      
         output:
         tuple val(pair_id), path("${pair_id}_reversed_adjusted.fasta")
          
         script:
         """
         ml load Python/3.9.5-GCCcore-10.3.0-bare
          
         starting_postition="\$(grep -B 1 'COX1' ${mito_gbf} | grep ' gene' | grep -oP '\\d+(?=\\.\\.)' | tail -1)";
          
         python /*/rearrange_mito.py ${mito_fasta} "\$(grep -B 1 'COX1' ${mito_gbf} | grep ' gene' | grep -oP '\\d+(?=\\.\\.)' | tail -1)"
         """
      }
 
python script

      import argparse
      
      
      def modify_fasta_file(fasta_file, cox1_start_position):
          with open(fasta_file, "r") as file:
              lines = file.readlines()
      
          header = lines[0].strip()
          sequence = "".join(lines[1:]).replace("\n", "")
      
          # Circularize the mitochondrial sequence with the COX1 gene at the desired starting position
          adjusted_sequence = sequence[cox1_start_position:] + sequence[:cox1_start_position]
      
          # Write the adjusted sequence to a new file
          adjusted_file = fasta_file.replace(".fasta", "_adjusted.fasta")
          with open(adjusted_file, "w") as file:
              file.write(f"{header}\n{adjusted_sequence}")
      
          return adjusted_file
      
      # Provide the path to your GetOrganelle-generated mitochondrial genome in FASTA format
      fasta_file = "Mitobim_test_original.fasta"
      parser = argparse.ArgumentParser()
      parser.add_argument("fasta_file", help="Path to the FASTA file")
      parser.add_argument("starting_position", type=int, help="Starting position of the COX1 gene")
      args = parser.parse_args()
      
      adjusted_fasta_file = modify_fasta_file(args.fasta_file, args.starting_position)
      print("Adjusted FASTA file:", adjusted_fasta_file) 

### Re-annotation of newly adjusted starting position - process mitoz_assembly_adjusted <a name="MitoZ3"></a>
Mitochondrial genome is re-annotated after adjustment to COX1 as starting position.

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
      
          /*/mitoz/MitoZ_v3.6.sif mitoz annotate --clade Nematoda --fastafiles ${pair_id}_tmp_header.fasta --outprefix MitoZ;
          
          cp ${pair_id}_tmp_header.fasta ${pair_id}.fasta
          """
      }

### Visualization of mitochondrial genomes using MitoZ - process mitoz_visualisation <a name="MitoZ4"></a>
Mitochondrial genomes were visualized using MitoZ, including coverage plots in circos with the raw reads.

      process mitoz_visualisation {
          maxForks 60
          stageInMode 'copy'
          publishDir '/*/mitochondrial_baiting_and_analysis_pipeline/Results_samples/CDS_and_non_CDS/', mode: 'copy'
          
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
          singularity run /*/mitoz/MitoZ_v3.6.sif mitoz visualize \
          --circos /*/circos-0.69-9/bin/circos \
          --gb ${mito_gbf} --outdir ./results --run_map yes \
          --fq1 ${reads[0]} \
          --fq2 ${reads[1]};
      
          cp ./results/circos.png ./${pair_id}_mito_circos.png
          """
      }


### Extraction of gene sequences and MAFFT alignment of individual genes - process mafft_new_sequences <a name="Mafft"></a>
Coding sequences were extracted from the generated gbf files using a custom python script and saved by gene in the
 fasta format, before alignment with MAFFT.
These sequences are only from human infecting *Trichuris* species without references.

      process mafft_new_sequences {
         input:
         path(mito_gbf)
      
         script:
         """
         ls | grep "TRICHURISP-" | xargs rm -f;
         ml load Python/3.9.5-GCCcore-10.3.0-bare
         ml load MAFFT/7.490-GCC-10.3.0-with-extensions
         cp /results_nucleic_acid_phylogeny/generate_nucleicacid_fastas.py ./
         python generate_nucleicacid_fastas.py;
         for file in *.fasta; do
         sed -i 's/^>/>s/' "\$file";
         sed -i 's/\\[//' "\$file";
         sed -i 's/\\]//' "\$file";
         mafft --auto "\$file" > "\$file"aligned.fasta
         done 
      
      
          """
      }

python script

      from Bio import SeqIO
      import os
      
      # Define the folder containing your GBF files
      folder_path = "./"
      
      # Dictionary to store gene sequences
      gene_sequences = {}
      
      # Process all GBF files in the folder
      for file_name in os.listdir(folder_path):
          if file_name.endswith(".gbf"):
              file_path = os.path.join(folder_path, file_name)
              for record in SeqIO.parse(file_path, "genbank"):
                  for feature in record.features:
                      if feature.type == "CDS" and "gene" in feature.qualifiers:
                          gene_name = feature.qualifiers["gene"][0]
                          gene_sequence = str(feature.location.extract(record).seq)
                          if gene_name not in gene_sequences:
                              gene_sequences[gene_name] = {}
                          gene_sequences[gene_name][record.id] = gene_sequence
      
      # Create FASTA files for each gene
      for gene, sequences in gene_sequences.items():
          fasta_file_name = os.path.join(folder_path, f"{gene}.fasta")
          with open(fasta_file_name, 'w') as fasta_file:
              for taxon, sequence in sequences.items():
                  fasta_file.write(f">{taxon}\n{sequence}\n")
      
      print("FASTA files have been created for each gene.")


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

