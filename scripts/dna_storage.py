def parse_fasta(file_path):
    """
    Parse a FASTA file into a dictionary with headers as keys and sequences as values.
    """
    sequences = {}
    current_header = None
    headers = []
    current_sequence = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                if current_header:
                    sequences[current_header] = ''.join(current_sequence)
                current_header = line.strip()
                current_sequence = []
            else:
                current_sequence.append(line.strip())
        if current_header:
            sequences[current_header] = ''.join(current_sequence)
    return sequences


def preprocess_sequences(sequences):
    """
    Create a dictionary mapping genome IDs (including .mid) to their sequences for faster lookup.
    """
    genome_map = {}
    for header, sequence in sequences.items():
        try:
            if "|" in header:
                genome_id = header.split('|')[1].strip()  # Extract genome ID
            else:
                genome_id = header.lstrip('>').strip()  # Use the entire header without '>'
            genome_map[genome_id] = sequence
        except IndexError:
            print(f"[ERROR] Could not extract genome ID from header: {header}")
    print(f"[DEBUG] Preprocessed {len(genome_map)} genomes.")
    return genome_map


def match_genes_to_sequences(genes, genome_map):
    """
    Match genes to their position in the specified genome based on IDs in the headers,
    while ensuring sequence order is respected.
    """
    matches = []
    current_positions = {}  # Track the current position for each genome

    for idx, (gene_header, gene_seq) in enumerate(genes.items()):
        # Extract the genome ID from the gene header
        try:
            if "|" in gene_header:
                header_parts = gene_header.split('|')
                genome_id = header_parts[1].strip()  # Extract genome ID after '|'
            else:
                genome_id = gene_header.lstrip('>').strip()  # Use the entire header without '>'

            # Append ".mid" to the genome ID if not already present
            if not genome_id.endswith('.mid'):
                genome_id += '.mid'
        except (IndexError, ValueError):
            print(f"[ERROR] Could not extract genome ID from gene header: {gene_header}")
            continue

        # Find the corresponding genome sequence
        genome_sequence = genome_map.get(genome_id, None)
        if genome_sequence:
            # Normalize case
            gene_seq = gene_seq.strip().upper()
            genome_sequence = genome_sequence.strip().upper()

            # Start search from the last position in this genome
            start_pos = current_positions.get(genome_id, 0)

            # Perform the search
            match_pos = genome_sequence.find(gene_seq, start_pos)

            if match_pos != -1:
                print(f"[DEBUG] Found gene {gene_header} in genome {genome_id} at position {match_pos}.")
                matches.append({
                    "gene_header": gene_header,
                    "gene_sequence": gene_seq,
                    "genome_id": genome_id,
                    "position": match_pos
                })
                # Update current position for this genome
                current_positions[genome_id] = match_pos + len(gene_seq)
            else:
                # Move the position forward to avoid infinite loops
                next_pos = start_pos + 1
                current_positions[genome_id] = next_pos
                print(f"[DEBUG] Gene {gene_seq} not found in genome {genome_sequence} starting from position {start_pos}. Moving to {next_pos}.")
                matches.append({
                    "gene_header": gene_header,
                    "gene_sequence": gene_seq,
                    "genome_id": genome_id,
                    "position": None
                })
        else:
            print(f"[ERROR] Genome ID {genome_id} not found in genome map.")
            matches.append({
                "gene_header": gene_header,
                "gene_sequence": gene_seq,
                "genome_id": None,
                "position": None
            })

        # Progress logging
        if idx % 100 == 0:  # Log every 100 genes
            print(f"[PROGRESS] Processed {idx + 1}/{len(genes)} genes.")

    return matches


def write_matches_to_fasta(matches, output_file):
    """
    Write the matched genes to a new FASTA file with start position and genome ID.
    """
    with open(output_file, 'w') as file:
        for match in matches:
            if match['genome_id'] is not None and match['position'] is not None:
                header = f"{match['gene_header']} | Midi File ID: {match['genome_id']} | Start Position: {match['position']}"
                file.write(f">{header}\n{match['gene_sequence']}\n")
            else:
                header = f"{match['gene_header']} | NO MATCH FOUND"
                file.write(f">{header}\n{match['gene_sequence']}\n")


# File paths
genes_file = '../data/all_notes.fa'  # Genes with IDs in headers
entire_sequences_file = '../data/all_midi.fa'  # Genomes (full sequences)
output_file = 'output_matched_genes.fasta'

# Parse FASTA files
genes = parse_fasta(genes_file)
sequences = parse_fasta(entire_sequences_file)

# Preprocess genome sequences for fast lookup
genome_map = preprocess_sequences(sequences)

# Perform matching
matches = match_genes_to_sequences(genes, genome_map)

# Write matches to a new FASTA file
write_matches_to_fasta(matches, output_file)

print(f"[INFO] Matches have been written to {output_file}")
