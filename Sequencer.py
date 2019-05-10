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
		self.show_params()

		graph = self.construct_graph()
		print(len(graph.items()))
		[print(key.sequence, [v.sequence for v in value]) for key, value in graph.items()]
		print(self.get_start_seq("odd"), self.get_start_seq("even"))

		curr_path = [self.find_oligo(graph, self.get_start_seq("odd"))] # Initialized as odd path
		last_path = [self.find_oligo(graph, self.get_start_seq("even"))] # Initialized as even path
		DNA = []

		# print("curr_path", curr_path)
		# print("last_path", last_path)

		while len(DNA) <= self.n:
			candidates = [oligo for oligo in graph[curr_path[-1]] if oligo.visited == False]

			# print("Candidates:", candidates)

			for candidate in candidates:
				if last_path[-1].sequence[1:] + candidate.sequence[-1] in S2:
					# print("Validation succes:", last_path[-1][1:] + candidate.sequence[-1])
					curr_path[-1].visited = True
					DNA.append(candidate.sequence)
					curr_path.append(candidate)
					curr_path, last_path = curr_path, last_path 
					break;

		return DNA


	def construct_graph(self):
		graph = defaultdict(list)

		for oligo in self.S1:
			for another_oligo in self.S1:
				if oligo.overlap(another_oligo, 1) and oligo != another_oligo:
					graph[oligo].append(another_oligo)

		return graph


	def get_start_seq(self, parity="odd"):
		if parity == "odd":
			node = list(self.start[0:-1])
		else:
			node = list(self.start[1:])
			
		for i in range(0, 2*self.k - 1):
			if i % 2 == 1:
				node[i] = "x"
					
		return "".join(node)


	def find_oligo(self, graph, sequence):
		for oligo in graph:
			if oligo.sequence == sequence:
				return oligo;

		return None;


	def show_params(self):
		print("DNA length:", self.n)
		print("Oligo length:", self.k)
		print("Starting sequence:", self.start)
