from Oligo import Oligo
from collections import defaultdict

class Sequencer:
	def __init__(self, S1, S2, n, k, start):
		self.S1 = S1
		self.S2 = S2
		self.n = n
		self.k = k
		self.start = start


	def run(self):
		graph = self.construct_graph()
		curr_path = [self.find_oligo(graph, self.get_start_seq("odd"))] # Initialized as odd path
		last_path = [self.find_unfinished_oligo(graph, self.get_start_seq("even"))] # Initialized as even path
		DNA = []

		self.show_params()

		max_steps = self.n - (2*self.k - 1) - 1

		while len(DNA) < max_steps:
			curr_oligo = curr_path[-1]
			last_oligo = last_path[-1]

			print("DNA", DNA)

			candidates = self.getCandidates(graph, curr_oligo)

			for candidate in candidates:
				if self.isValid(last_oligo, curr_oligo):
					self.appendOligo(DNA, curr_path, candidate)
					curr_path, last_path = last_path , curr_path
					break;

		return "".join(DNA)


	def getCandidates(self, graph, oligo):
		return [candidate for candidate in graph[oligo] if candidate.visited == False]


	def isValid(self, last_o, cand_o):
		return (last_o.sequence[2:] + cand_o.sequence[-1]) in self.S2


	def appendOligo(self, result_seq, path, oligo):
		result_seq.append(oligo.sequence[-1])
		current_path.append(oligo)
		oligo.visited = True


	def construct_graph(self):
		graph = defaultdict(list)

		for oligo in self.S1:
			for another_oligo in self.S1:
				if oligo.overlap(another_oligo, 1) and oligo != another_oligo:
					graph[oligo].append(another_oligo)

		return graph


	def get_start_seq(self, parity="odd"):
		if parity == "odd":
			node = list(self.start[0:])
		else:
			node = list(self.start[1:-1])
			
		for i in range(0, len(node)):
			if i % 2 == 1:
				node[i] = "X"
					
		return "".join(node)


	def find_oligo(self, graph, sequence):
		for oligo in graph:
			if oligo.sequence == sequence:
				return oligo;

		return None


	def find_unfinished_oligo(self, graph, sequence):
		nucleobases = ['A','C','T','G']
		node = None

		for n in nucleobases:
			node = self.find_oligo(graph, sequence + 'X' + n)
			if node != None:
				break

		return node


	def show_params(self):
		print("DNA length:", self.n)
		print("Oligo length:", self.k)
		print("Starting sequence:", self.start, " length:", len(self.start))