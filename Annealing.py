import random
from Oligo import Oligo


class Annealing:
    def __init__(self, oli_set, oli_length, DNA_length, odd, even):
        self.odd = odd
        self.even = even
        self.current_path = 'odd'
        
        self.used_olis = self.compute_olis_dict(oli_set)                        #{oli: bool}
        self.free = self.compute_free()                             #lista wolnych oligosów
        self.dist = self.compute_distances(oli_set)              #słownik dystansów
    
        self.LENGTH = DNA_length        #docelowa długość łańucha DNA
        self.NUM_OLIGO = DNA_length - oli_length + 1              #docelowa liczba użytych oligonukleotydów
        self.T = 100                #temperatura
        self.max_dist = oli_length//2+1 #maksymalna odległość między oligonukleotydami

    def compute_olis_dict(self):
        used_olis = {oli: False for oli in oli_set}
        for oli in odd:
            used_olis[oli] = True
        for oli in even:
            used_olis[oli] = True
        return used_olis

    def compute_free(self):
        return [oli for (oli, used) in self.used_olis if !used]

    def compute_distances(self, olis):
        ret = {}    #dict of dicts: {oli_n: {k: [o1, o2, ...]}}
        for o1 in olis:
            d = {k:[] for k in range(1, self.max_dist+1)}
            for i, o2 in enumerate(olis):
                if o1 != o2: d[dist(o1, o2)].append(i)
            ret[o1] = d
        return ret 




    def mutate(self):
        action = random.choice(self.delete, self.add, self.swap)
        temp_state = action()
        temp_val = self.eval_state(temp_state)
        if self.decide(temp_val, self.current_val):
            self.current_val = temp_val
            self.odd, self.even = temp_state

    def delete(self):
        odd, even = self.new_temp_state()
        index, to_delete = self.pick_random()
        [del odd[index], del even[index]][self.current_path == 'even']
        self.used_olis[to_delete] = False
        self.free.append(to_delete)
        return (odd, even)

    def add(self):
        odd, even = self.new_temp_state()
        index, _ = self.pick_random()
        chosen = self.pop_random_from_free()
        [odd.insert(index, chosen), even.insert(index,chosen)][self.current_path == 'even']
        return (odd, even)

    def swap(self):
        odd, even = self.new_temp_state()
        permitted_dist = self.compute_permitted_dist(self.T)
        index, to_swap = self.pick_random()
        chosen = self.choose_neighbour(to_swap, permitted_dist)
        if chosen:
            self.used_olis[to_swap] = False; self.used_olis[chosen] = True
            self.free.append(to_swap); self.free.remove(chosen)
            [odd[index] = chosen, even[index] = chosen][self.current_path == 'even']
        return (odd, even)
    
    def new_temp_state(self):
        if self.current_path == 'odd':
            odd = self.odd[:]
            even = self.even
        else:
            odd = self.odd
            even = self.even[:]
        return (odd, even)

    def pick_random(self):
        if self.current_path == 'odd':
            i = random.randrange(len(self.odd))
            return i, odd[i]
        i = random.randrange(len(self.even))
        return i, even[i]

    def choose_neighbour(self, oli, max_d):
        possible = set([*self.dist[oli][k] for k in range(1, max_d+1)]) & self.free
        if possible == []:
            return None
        return self.used_olis.keys()[random.choice(possible)]   #TODO: to moze być unsafe jeśli pozycje w dict sie zmieniają
                                                                # lepiej użyć listy olis

    def pop_random_from_free(self):
        i = random.randrange(len(self.free))
        oli = self.free[i]
        self.used_olis[oli] = True
        return self.free.pop(i)

    def eval_state(self, state):
        odd, even = state
        return abs(self.compute_curr_length() - self.LENGTH) +\
                abs(len(odd) + len(even) - self.NUM_OLIGO)





def dist(o1, o1):
    return sum([a!=b for a,b in zip(o1, o2)])

