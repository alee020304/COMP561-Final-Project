import itertools
import json
import math
from Bio import SeqIO


# Constants
STATES = ["header", "time_delta", "on", "note", "velocity", "intergenic"]

# Initialize Transition Matrix
transition_matrix = {
    state: {target_state: 0 for target_state in STATES}
    for state in STATES
}
# Define valid transitions
transition_matrix["time_delta"]["on"] = 1
transition_matrix["on"]["note"] = 1
transition_matrix["note"]["velocity"] = 1
transition_matrix["velocity"]["intergenic"] = 1

header_lengths = []  # To store lengths of headers
intergenic_lengths = []  # To store lengths of intergenic regions

# Initialize nucleotide counts and codon frequency dictionary
header_nuc_count = {"A": 0, "T": 0, "G": 0, "C": 0}
intergenic_nuc_count = {"A": 0, "T": 0, "G": 0, "C": 0}
state_codon_freq = {
    "time": {''.join(codon): 0 for codon in itertools.product("ATGC", repeat=4)},
    "on": {''.join(codon): 0 for codon in itertools.product("ATGC", repeat=4)},
    "note": {''.join(codon): 0 for codon in itertools.product("ATGC", repeat=4)},
    "velocity": {''.join(codon): 0 for codon in itertools.product("ATGC", repeat=4)},
}
delta_codon_freq = {''.join(codon): 0 for codon in itertools.product("ATGC", repeat=4)}
on_codon_freq = {''.join(codon): 0 for codon in itertools.product("ATGC", repeat=4)}
note_codon_freq = {''.join(codon): 0 for codon in itertools.product("ATGC", repeat=4)}
velocity_codon_freq = {''.join(codon): 0 for codon in itertools.product("ATGC", repeat=4)}


# Initialize Emission Matrix
emission_matrix = {
    state: {
        **{nucleotide: 0 for nucleotide in "AGCT"},  # Single nucleotides for "header"
        **{''.join(codon): 0 for codon in itertools.product("ATGC", repeat=4)}  # Tetramers for other states
    }
    for state in STATES
}

def parse_genome(file_path):
    """
    Parse a genome FASTA file and return a dictionary with IDs (excluding .mid) as keys,
    and a dictionary containing sequences and their lengths as values.
    """
    genome_sequences = {}

    with open(file_path, 'r') as fasta_file:
        current_id = None
        current_sequence = []

        for line in fasta_file:
            line = line.strip()
            if line.startswith(">"):
                # Save the previous record
                if current_id:
                    sequence = ''.join(current_sequence)
                    genome_sequences[current_id] = {
                        "sequence": sequence,
                        "length": len(sequence)
                    }

                # Parse the new ID (exclude ".mid")
                current_id = line[1:].replace(".mid", "").strip()
                current_sequence = []
            else:
                # Append sequence lines
                current_sequence.append(line)

        # Save the last record
        if current_id:
            sequence = ''.join(current_sequence)
            genome_sequences[current_id] = {
                "sequence": sequence,
                "length": len(sequence)
            }

    return genome_sequences

def parse_midi(file_path):
    """
    Parse a FASTA file and extract start and stop positions for each ID.
    For each entry, assume the stop position is start + 12.
    Returns a dictionary where each key is an ID, and its value is a list of dictionaries with 'start' and 'stop'.
    """
    results = {}

    with open(file_path, 'r') as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if line.startswith(">"):
                # Parse current_id
                current_id_parts = line[1:].split("|")
                
                # Extract the genome identifier (e.g., "mz_570_1")
                if len(current_id_parts) > 1:
                    genome_id = current_id_parts[1].strip()
                else:
                    print(f"Warning: Unable to parse genome ID from line: {line}")
                    continue
                
                # Extract the start position
                start_position = None
                if "Start Position:" in current_id_parts[-1]:
                    try:
                        start_position = int(current_id_parts[-1].split(":")[1].strip())
                    except ValueError:
                        print(f"Warning: Unable to parse start position from line: {line}")
                        continue
                
                # Calculate stop position
                if start_position is not None:
                    stop_position = start_position + 12
                
                    # Store start and stop in results
                    if genome_id not in results:
                        results[genome_id] = []
                    results[genome_id].append({"start": start_position, "stop": stop_position})
                else:
                    print(f"Warning: Start position missing in line: {line}")

    return results





import itertools

def calculate_frequencies(genome, matches):
    """
    Calculate nucleotide counts for headers, intergenic regions, and codon frequencies for genes,
    with codon frequencies separated into different states: time, on, note, velocity.
    Also calculates and prints the average lengths of headers and intergenic regions.
    
    Args:
    - genome: Path to the genome FASTA file.
    - matches: Path to the gene matches file.
    """
    # Initialize totals for lengths
    total_header_length = 0
    total_header_count = 0
    total_intergenic_length = 0
    total_intergenic_count = 0

    # Parse genome and matches
    genome_dict = parse_genome(genome)
    gene_dict = parse_midi(matches)

    for genome_id, genome_data in genome_dict.items():
        sequence = genome_data["sequence"]
        genes = gene_dict.get(genome_id, [])
        # Sort genes by start position
        genes = sorted(genes, key=lambda x: x["start"])
        
        # Process header (everything before the first "time" region)
        if genes:
            first_time_start = genes[0]["start"] - 4  # "time" starts 4 nucleotides before the first gene
            header_sequence = sequence[:first_time_start]
            total_header_length += len(header_sequence)
            total_header_count += 1

            for nucleotide in header_sequence:
                if nucleotide in header_nuc_count:
                    header_nuc_count[nucleotide] += 1

        # Debug: Print genome and first few genes
        print(f"Genome ID: {genome_id}")
        for gene in genes[:3]:  # Show only the first 3 genes for brevity
            print(f" Here is the gene positions: Start: {gene['start']}, Stop: {gene['stop']}")

        # Process intergenic and gene regions
        prev_stop = 0
        for gene in genes:
            # Intergenic region: everything between the previous stop and the current "time" start
            time_start = gene["start"] - 4
            intergenic_sequence = sequence[prev_stop:time_start]
            total_intergenic_length += len(intergenic_sequence)
            total_intergenic_count += 1

            for nucleotide in intergenic_sequence:
                if nucleotide in intergenic_nuc_count:
                    intergenic_nuc_count[nucleotide] += 1
            
            # Gene region: separate processing for "time", "on", "note", and "velocity"
            time_sequence = sequence[time_start:time_start + 4]
            on_sequence = sequence[gene["start"]:gene["start"] + 4]
            note_sequence = sequence[gene["start"] + 4:gene["start"] + 8]
            velocity_sequence = sequence[gene["start"] + 8:gene["start"] + 12]

            # Update state codon frequencies
            state_codon_freq["time"][time_sequence] = state_codon_freq["time"].get(time_sequence, 0) + 1
            state_codon_freq["on"][on_sequence] = state_codon_freq["on"].get(on_sequence, 0) + 1
            state_codon_freq["note"][note_sequence] = state_codon_freq["note"].get(note_sequence, 0) + 1
            state_codon_freq["velocity"][velocity_sequence] = state_codon_freq["velocity"].get(velocity_sequence, 0) + 1

            prev_stop = gene["stop"]  # Update the previous stop position

        # Process final intergenic region (after the last gene)
        final_intergenic_sequence = sequence[prev_stop:]
        total_intergenic_length += len(final_intergenic_sequence)
        total_intergenic_count += 1
        for nucleotide in final_intergenic_sequence:
            if nucleotide in intergenic_nuc_count:
                intergenic_nuc_count[nucleotide] += 1

    # Calculate and print average lengths
    avg_header_length = total_header_length / total_header_count if total_header_count > 0 else 0
    avg_intergenic_length = total_intergenic_length / total_intergenic_count if total_intergenic_count > 0 else 0

    print(f"\nAverage Header Length: {avg_header_length}")
    print(f"Average Intergenic Length: {avg_intergenic_length}")

def calculate_emission_matrix():
    """
    Calculate the emission matrix for each state, maintaining the structure:
    - Header and intergenic states: Include nucleotide frequencies and empty codon frequencies.
    - Other states: Include codon frequencies and empty nucleotide frequencies.
    """
    emission_matrix = {}

    # Header and Intergenic: Use nucleotide counts with empty codon frequencies
    for state in ["header", "intergenic"]:
        freq = header_nuc_count if state == "header" else intergenic_nuc_count
        total = sum(freq.values())
        if total > 0:
            nucleotide_probs = {nucleotide: count / total for nucleotide, count in freq.items()}
        else:
            nucleotide_probs = {nucleotide: 0 for nucleotide in "AGCT"}

        codon_probs = {''.join(codon): 0 for codon in itertools.product("ATGC", repeat=4)}
        emission_matrix[state] = {**nucleotide_probs, **codon_probs}

    # Other states: Include codon frequencies with empty nucleotide frequencies
    for state in ["time", "on", "note", "velocity"]:
        codon_freq = state_codon_freq[state]
        total_codons = sum(codon_freq.values())
        if total_codons > 0:
            codon_probs = {codon: count / total_codons for codon, count in codon_freq.items()}
        else:
            codon_probs = {''.join(codon): 0 for codon in itertools.product("ATGC", repeat=4)}

        nucleotide_probs = {nucleotide: 0 for nucleotide in "AGCT"}
        emission_matrix[state] = {**nucleotide_probs, **codon_probs}

    return emission_matrix



# Save matrices to JSON files
def save_matrix_to_json(matrix, file_name):
    with open(file_name, 'w') as file:
        json.dump(matrix, file, indent=4)

# Main function
def main():
    # Define file paths
    genome_file = "../data/all_midi.fa"  # Replace with your genome file
    matches_file = "output_matches.fasta"  # Replace with your matches file

    # Process the genome and matches to populate matrices
    calculate_frequencies(genome_file, matches_file)

    # Save matrices to JSON files
    emission_matrix = calculate_emission_matrix()
    save_matrix_to_json(emission_matrix, "emission_matrix.json")


# Execute main
if __name__ == "__main__":
    main()
