from Oligo import Oligo
from collections import defaultdict

class Sequencer:
	def __init__(self, s1, s2, length, start):
		self.s1 = s1
		self.s2 = s2
		self.length = length
		self.start = start

	def construct_graph(self):
		graph = defaultdict(list)

		for oligo in self.s1:
			for another_oligo in self.s1:
				if oligo.overlap(another_oligo, 1) and oligo != another_oligo:
					graph[oligo].append(another_oligo)

		return graph


	def run():
		graph = self.construct_graph()

		odd_path = [graph[]]
		even_path = []

	def get_start_seq(parity="odd"):
		node = []

		if parity == "odd":
			x, y = 0, 1		
		else:
			x, y = 1, 0

		for i in range(x, self.start - y):
			if i % 2 == y:
				node.append("x")
			else:
				node.append(self.start[i])

		return "".join(node)



