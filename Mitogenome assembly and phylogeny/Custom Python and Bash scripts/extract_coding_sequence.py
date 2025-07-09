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

