from Oligo import Oligo


class Annealing:
    def __init__(self, oli_set, length, start):
        self.start = start
        self.length = length
        self.olis = list(olis)
        self.dist = compute_distances(olis)
        self.T = 100
        self.max_dist = length//2+1

    def mutate():
        pass


    def eval_state(self, state):
        return abs(len(state) - self.length) +\
                abs(state.num_oligos



def compute_distances(olis, length):
    ret = {}    #dict of dicts: {oli_n: {k: [o1, o2, ...]}}
    max_dist = length//2+1
    for o1 in olis:
        d = {k:[] for k in range(max_dist+1)}
        for i, o2 in enumerate(olis):
            d[dist(o1, o2)].append(i)
        ret[o1] = d
    return ret 


def dist(o1, o1):
    return sum([a!=b for a,b in zip(o1, o2)])

