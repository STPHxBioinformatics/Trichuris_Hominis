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
