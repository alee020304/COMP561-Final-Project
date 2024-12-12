def parse_fasta(file_path):
    """
    Parse a FASTA file into a generator of (header, sequence) tuples.
    This handles large files efficiently.
    """
    current_header = None
    current_sequence = []
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if current_header:
                    yield current_header, ''.join(current_sequence)
                current_header = line
                current_sequence = []
            else:
                current_sequence.append(line)
        if current_header:
            yield current_header, ''.join(current_sequence)


def match_genes_to_genomes(genome_file, gene_file, output_file):
    """
    Parse the genome and genes, match genes to all genomes, and save results to a new FASTA file.
    Prevent duplicates by ensuring unique start positions for each match.
    """
    with open(output_file, 'w') as output:
        # Parse all genomes
        genome_iterator = parse_fasta(genome_file)

        # Parse all genes into a list for reuse
        gene_list = list(parse_fasta(gene_file))

        for genome_header, genome_sequence in genome_iterator:
            genome_id = genome_header.lstrip('>').strip()  # Remove '>' and extract genome ID

            # Track used positions in the current genome
            used_positions = set()

            # Match each gene to the current genome
            for gene_header, gene_sequence in gene_list:
                # Extract the gene's genome ID from the header
                if "|" in gene_header:
                    gene_genome_id = gene_header.split('|')[1].strip()  # Extract genome ID from gene header
                else:
                    gene_genome_id = gene_header.lstrip('>').strip()

                # Append ".mid" to match genome ID format
                if not gene_genome_id.endswith('.mid'):
                    gene_genome_id += '.mid'

                # Skip genes that don't match the current genome
                if gene_genome_id != genome_id:
                    continue

                # Search for the gene sequence in the genome
                start_pos = -1
                while True:
                    # Find the next occurrence of the gene sequence
                    start_pos = genome_sequence.find(gene_sequence.upper(), start_pos + 1)

                    # Break if no further matches are found
                    if start_pos == -1:
                        print(f"[WARNING] Gene {gene_header} not found in genome {genome_id}.")
                        break

                    # Check if this position is already used
                    if start_pos not in used_positions:
                        # Mark position as used and write the match to the output file
                        used_positions.add(start_pos)
                        output_header = f"{gene_header} | Genome: {genome_id} | Start Position: {start_pos}"
                        output.write(f">{output_header}\n{gene_sequence}\n")
                        break  # Move to the next gene

    print(f"[INFO] Matches saved to {output_file}")


# File paths
genome_file = '../data/all_midi.fa'  # Replace with your genome file path
gene_file = '../data/all_notes.fa'     # Replace with your gene file path
output_file = 'output_matches.fasta'

# Run the matching process
match_genes_to_genomes(genome_file, gene_file, output_file)
