from .base_algorithm import Storage_Algorithm

'''
SA_Naive
- encode: loop over 2 bits at a time and translate the dibit to a corresponding nucleotide
- decode: loop over 4 nucleotides at a time and convert to a byte -> repeat for all nucleotides and return array of decoded bytes
'''
class SA_Naive(Storage_Algorithm):
	NUCLEOTIDES_ARR = ["A", "T", "C", "G"]

	def dibit_generator(self, filename):
		with open(filename, "rb") as f:
			byte = f.read(1)
			while byte:
				int_value = int.from_bytes(byte, byteorder="big")
				for idx in range(0, 8, 2):
					dibit = (int_value >> (6 - idx)) & 0b11
					yield dibit
				byte = f.read(1)

	def encode(self, dna_storage, filename):
		dna_storage.filename = filename
		dna_storage.seq = ""
		for dibit in self.dibit_generator(filename):
			dna_storage.seq += self.NUCLEOTIDES_ARR[dibit]
		return dna_storage.seq

	def decode(self, dna_storage):
		decoded = []
		for idx in range(0, len(dna_storage.seq), 4):
			byte_worth_of_nucleotides = map(lambda nucleotide: self.NUCLEOTIDES_ARR.index(nucleotide), dna_storage.seq[idx:idx+4])
			result = 0
			for idx, val in enumerate(byte_worth_of_nucleotides):
				result |= (val << (6 - idx*2))
			decoded.append(result)
		return bytearray(decoded)

