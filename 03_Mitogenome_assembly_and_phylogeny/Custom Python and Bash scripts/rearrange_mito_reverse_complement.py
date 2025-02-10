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
