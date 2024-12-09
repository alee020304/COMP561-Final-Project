class Storage_Algorithm:

	# default name of a given algorithm will be its class name (can override in subclass if it seems unclear)
	@property
	def name(self):
		return self.__class__.__name__

	def encode(self, dna_storage, filename):
		raise NotImplementedError("Encode method not implemented")

	def decode(self, dna_storage):
		raise NotImplementedError("Decode method not implemented")
