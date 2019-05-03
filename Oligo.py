class Oligo:
	def __init__(self, sequence):
		self.sequence = sequence
		self.length = len(sequence)

		def overlap_by_one(another_oligo):
			if self.sequence[0:-3] == another_oligo.sequence[0:-3] or
			   self.sequence[0:-3] == another_oligo.sequence[2:] or
			   self.sequence[2:] == another_oligo.sequence[0:-3] or
			   self.sequence[2:] == another_oligo.sequence[2:]:
				return True
			else:
				return False
			
