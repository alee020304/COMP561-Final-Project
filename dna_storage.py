from algos import Storage_Algorithm, SA_Naive
import os
import random

# Illumina sequencing
class Mutation_Simulator:
	def __init__(self, mutation_rates):
		self.mutation_rates = mutation_rates  # Dictionary of mutation types and rates

	# TODO: implement more realistic mutation simulation (find research paper that has probabilities and such?)
	def apply_mutations(self, sequence):
		mutated_sequence = list(sequence)
		mutation_count = 0
		for i in range(len(mutated_sequence)):
			mutation_type = random.choices(list(self.mutation_rates.keys()), weights=list(self.mutation_rates.values()))[0]
			if random.random() < self.mutation_rates[mutation_type]:
				if mutation_type == "insertion":
					mutated_sequence.insert(i, random.choice("ACGT"))
				elif mutation_type == "deletion":
					if i < len(mutated_sequence) -1:
						mutated_sequence.pop(i)
				elif mutation_type == "substitution":
					mutated_sequence[i] = random.choice("ACGT")
				mutation_count += 1
		return "".join(mutated_sequence), mutation_count


class DNA_Storage:
	def __init__(self, algo, mutation_sim=None):
		self.algo = algo
		self.mutation_sim = mutation_sim
		self.seq = ""
		if not isinstance(self.algo, Storage_Algorithm):
			raise Exception("Invalid algorithm provided")

	def simulate_mutations(self):
		if not self.seq:
			print("No DNA sequence generated, please encode first")
			return -1

		if not self.mutation_sim:
			print("No mutation simulator set")
			return -1

		print("Simulating mutations...")
		self.seq, mutation_count = self.mutation_sim.apply_mutations(self.seq)
		print(f"{mutation_count} ({(mutation_count/len(self.seq))*100:.2f}% mutated) mutations applied successfully.")


	# Set mutation simulator
	def set_mutation_simulator(self, mutation_simulator):
		if not isinstance(mutation_simulator, Mutation_Simulator):
			raise Exception("Invalid algorithm provided")
		self.mutation_simulator = mutation_simulator

	# Update algorithm
	def set_algo(self, algo):
		if not isinstance(algo, Storage_Algorithm):
			raise Exception("Invalid algorithm provided")
		self.algo = algo

	# Encode using the provided algorithm
	def encode(self, filename):
		print(f"Encoding sequence of {filename}")
		return self.algo.encode(self, filename)

	# Decode using the provided algorithm
	def decode(self, output):
		if not self.seq:
			print("No DNA sequence generated, please encode first")
			return -1

		print(f"Decoding sequence...")
		decoded = self.algo.decode(self)

		if not decoded:
			print("Failed to decode")
			return

		with open(output, "wb") as decoded_f:
			decoded_f.write(decoded)

		print(f"Saved decoded sequence to {output}")

	def export(self, output):
		if not self.seq:
			print("No DNA sequence generated, please encode first")
			return -1

		# TODO: do I trim filename stored in fasta (ie. samples/minion1.jpg -> minion1.jpg)
		with open(output, "w") as fastfa_f:
			fastfa_f.write(f">{self.filename} | {self.algo.name}\n")
			fastfa_f.write(self.seq)

	# TODO: add more information and also allow users to import saved DNA sequence if wanted
	def performance(self, debug=False):
		if not self.seq:
			print("No DNA sequence generated, please encode first")
			return -1

		# shouldn't need code below since self.seq will set a filename property
		# filename = self.filename if hasattr(self, "filename") else filename
		# if not filename:
		# 	print("Missing filename from either encoder or manually entered")

		orig_len, seq_len = os.path.getsize(self.filename) * 8, len(self.seq)
		compression_rate = seq_len / orig_len
		gc_ratio = (self.seq.count("G")+self.seq.count("C"))/seq_len
		if debug:
			print(f"Before: {orig_len // 4} bytes")
			print(f"After: {seq_len // 4} bytes")
			print(f"Compression rate of {1/compression_rate} bits/nt ({compression_rate * 100:.2f}% of original file)")
			print("__Extra__")
			print(f"GC content: {gc_ratio*100:.3f}% (should be within 40-60%)")
		return compression_rate


if __name__ == "__main__":
	naive_algo = SA_Naive()
	mutation_sim = Mutation_Simulator({
		'insertion': 0.01,
		'deletion': 0.01,
		'substitution': 0.01
	})
	a = DNA_Storage(naive_algo, mutation_sim=mutation_sim)
	a.encode("samples/minion1.jpg")
	a.simulate_mutations()
	# a.export("output/minion1.fa")
	# a.decode("output/minion1_decoded.jpg")
	# a.performance(debug=True)



'''

Optimally storing MP3 files in DNA

-


TODO:
- look into

'''
# Reed-Solomon


'''
interesting ideas to explore
- take advantage of dna structure
	- using expected structure for data integrity?
	- avoiding secondary structures that can interfere with sequencing


ideas:
- HMM based error correction for mutated sequence (somewhat more common, not sure what we can improve here)
- focus on specific file format and use information about the file type to optimize encoding for just that
- really investigate encoding optimization/factor in dna structure it would take on if sequenced/use that as factor (ie. avoiding secondary structures that can interfere with sequencing)

Module 2: Sequence Evolution and Alignment
Connection to DNA Storage:

When storing data in DNA, mutations (substitutions, insertions, deletions) may occur during synthesis or sequencing. These are analogous to sequence evolution processes.
Use sequence alignment techniques to align the original and mutated sequences to assess the integrity of the stored data and identify error patterns.
Potential Project Ideas:

Develop an error-correction algorithm inspired by global/local sequence alignment techniques.
Simulate the evolution of encoded DNA sequences over time and measure how effectively they retain the original data.


Module 3: Genome Sequencing and Assembly
Connection to DNA Storage:

Retrieving data stored in DNA involves sequencing and, in some cases, assembling fragments (similar to de novo genome assembly).
DNA storage systems must handle problems like read mapping, overlap detection, and reconstruction of the original sequence.
Potential Project Ideas:

Implement a simplified genome assembly algorithm for reconstructing stored data from fragmented DNA reads.
Explore how sequencing depth and read length impact the accuracy and efficiency of data recovery.


Module 6: Hidden Markov Models (HMMs)
Connection to DNA Storage:

HMMs are highly relevant for error detection and correction in DNA sequences. They can model the probability of mutations and sequencing errors based on context.
Profile HMMs, often used in gene prediction, can be adapted for detecting specific patterns in encoded data.
Potential Project Ideas:

Design an HMM-based error-correction algorithm that models mutations in DNA storage systems.
Use the Baum-Welch algorithm to estimate parameters for your mutation simulation model and evaluate its accuracy against real-world DNA data.
'''

'''

GC balanced and no homopolymers
- base rotation encoding scheme? (Goldman et al.)
- Randomization with a pseudorandom binary sequence is a popular strategy utilized to avoid homopolymers by many works 9,
- there isn't an universally accepted standard 


- DNA degradation


modulation decoding is highly efficient and extremely robust for the detection of insertions and deletions, which can correct up to ~40% errors.


Vladimir Levenshtein (1935-2017) proved an interesting fact about deletion and insertion errors already
in the 1960’s: if a code corrects deletion errors, it can
also correct an equal number of combination of deletions
and insertions (However, an efficient encoding/decoding
algorithm for correcting deletions does not necessarily
imply an efficient algorithm for correcting deletion and
insertion errors). Moreover, a substitution error can be
regarded as a deletion error followed by an insertion.
Therefore, it is reasonable to focus on deletion errors, as
correcting deletion errors implies correcting a combination of the three types of errors.

'''

'''
Dynamic Programming (Modules 2 and 6):

Implement efficient algorithms for encoding and decoding data using dynamic programming approaches,
similar to those used in sequence alignment or HMMs.



Markov Processes (Module 6):

Use Markov models to simulate long-term stability and the likelihood of mutation patterns, integrating real-world data on DNA stability.
'''

'''
ideas:
- simulate wetlab process of dna synthesis, storage and sequencing
- for common data types we can take advantage of fact of their file structure (ie jpg, png, mp4), otherwise we stick with 50% reduction by converting 2-bit -> 1 nucleotide
- huffman encoding


Encoding arbitrary data
To encode arbitrary data in DNA, the data is typically first converted into ternary (base 3) data rather than binary (base 2) data. Each digit (or "trit") is then converted to a nucleotide using a lookup table.
To prevent homopolymers (repeating nucleotides), which can cause problems with accurate sequencing, the result of the lookup also depends on the preceding nucleotide.
Using the example lookup table below, if the previous nucleotide in the sequence is T (thymine), and the trit is 2, the next nucleotide will be G (guanine).[7][8]

Trits to nucleotides (example)
Previous	0	1	2
T			A	C	G
G			T	A	C
C			G	T	A
A			C	G	T
Various systems may be incorporated to partition and address the data, as well as to protect it from errors.

One approach to error correction is to regularly intersperse synchronization nucleotides between the information-encoding nucleotides. These synchronization nucleotides can act as scaffolds when reconstructing the sequence from multiple overlapping strands.[8]

'''