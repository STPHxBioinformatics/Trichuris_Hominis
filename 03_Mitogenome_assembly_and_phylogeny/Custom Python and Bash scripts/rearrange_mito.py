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
