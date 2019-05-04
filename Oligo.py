class Oligo:
	def __init__(self, sequence):
		self.sequence = sequence
		self.length = len(sequence)

	def overlap(self, another_oligo, max_overlap):
		offset = 0

		while offset < self.length - 1 and offset / 2 <= max_overlap:
			if self.sequence[offset:] == another_oligo.sequence[0:-offset]:
				return offset / 2
			offset += 2

		return False
