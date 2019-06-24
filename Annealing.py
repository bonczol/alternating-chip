import random
import matplotlib.pyplot as plt
import math
import time
from typing import List, Dict, Set, Tuple, Optional

class Annealing:
    def __init__(self, oli_set, S2, oli_length, DNA_length, odd, even,
        max_time, cooling_factor, exponential_coef, clustering_factor):
        self.START_TIME = time.time()

        self.S1: List[str] = list(oli_set)                       #list
        self.S2: List[str] = list(S2)                            #zbiór oligosów potwierdzających
        self.oli2ind = {oli: i for i, oli in enumerate(self.S1)}
        self.LENGTH: int = DNA_length                            #docelowa długość łańucha DNA
        self.OLI_LENGTH: int = 2*oli_length-1                        #długość oligonukleotydu ze zbioru S1
        self.odd: List[str] = odd 
        self.even: List[str] = even

        self.NUM_OLIGO: int = self.LENGTH - self.OLI_LENGTH + 1        #docelowa liczba użytych oligonukleotydów
        self.max_dist: int = self.OLI_LENGTH//2+1                     #maksymalna odległość między oligonukleotydami

        self.used_olis = self.compute_olis_dict()         #{ind: bool}
        self.free: Set[int] = self.compute_free()                #lista wolnych oligosów
        self.dist = self.compute_distances()                     #słownik dystansów
        self.clustered_olis: List[int] = []
    
        self.current_path: List[str] = self.odd
        self.MAX_T = 10
        self.T: float = self.MAX_T                                      #temperatura

        self.current_val = self.eval_state((self.odd, self.even))
        self.best_found = [self.odd[:], self.even[:]]
        self.best_val = self.current_val

        self.EXPONENTIAL_COEF = exponential_coef
        self.CLUSTERING_FACTOR = clustering_factor
        self.COOLING_FACTOR = cooling_factor
        self.MAX_TIME = max_time               #seconds

    #-----------------------------------------------------------------------------
    #--------------------------Init computations----------------------------------
    #-----------------------------------------------------------------------------
    def compute_olis_dict(self) -> Dict[int, bool]:
        used_olis = {index: False for index in range(len(self.S1))}
        for oli in self.odd:
            used_olis[self.oli2ind[oli]] = True
        for oli in self.even:
            used_olis[self.oli2ind[oli]] = True
        return used_olis

    def compute_free(self) -> Set[int]:
        return set([oli_i for oli_i, used in self.used_olis.items() if not used])

    def compute_distances(self) -> Dict[str, Dict[int, List[int]]]:
        ret = {}    #dict of dicts: {oli_n: {k: [o1, o2, ...]}}
        for o1 in self.S1:
            d = {k:[] for k in range(1, self.max_dist+1)}
            for i, o2 in enumerate(self.S1):
                if o1 != o2: d[dist(o1, o2)].append(i)
            ret[o1] = d
        return ret 


    #-----------------------------------------------------------------------------
    #-------------------------------Main loop-------------------------------------
    #-----------------------------------------------------------------------------
    def run(self):
        vars, xs, bests = [], [], []
        i = 0
        while time.time() - self.START_TIME < self.MAX_TIME:
            self.mutate()
            for _ in range(self.CLUSTERING_FACTOR):
                self.try_to_cluster()
            if (1 - self.T/self.MAX_T) > random.random(): self.fix_current_path()
            self.current_path = [self.odd, self.even][i%2]
            self.cool_down()
            i += 1
        odd, even = self.best_found
        c = 100*sum([1/max(overlap(odd[i], odd[i+1]), 1)**2 for i in range(len(odd) -1)])/len(odd)
        d = 100*sum([1/max(overlap(even[i], even[i+1]), 1)**2 for i in range(len(even) -1)])/len(even)
        final_solution = self.construct_final_solution()
        return final_solution, odd, even, c, d, self.best_val

    #-----------------------------------------------------------------------------
    #--------------------------Temperature schedule and decision------------------
    #-----------------------------------------------------------------------------
    def cool_down(self):
        self.T *= self.COOLING_FACTOR
        
    def decide(self, val_old: float, val_new: float) -> bool:
        x = (val_old - val_new)/self.T
        try:
            return 1/(1 + math.exp(-self.EXPONENTIAL_COEF*x)) > random.random()
        except OverflowError:
            return True

    #-----------------------------------------------------------------------------
    #--------------------------Mutating-------------------------------------------
    #-----------------------------------------------------------------------------
    def mutate(self, option:str = None) -> None:
        '''Mutate current path; decide whether accept this change'''
        action = random.choice([self.delete, self.add, self.swap])
        if option is not None:
            action = {'del': self.delete, 'add': self.add}[option]
        changes, unchanges = action()
        if changes is not None:
            changes()
            temp_val = self.eval_state((self.odd, self.even))
            if temp_val <= self.best_val:
                self.update_best_state(temp_val)
            unchanges()
            if self.decide(temp_val, self.current_val):
                changes(confirm = True)
                self.current_val = temp_val

    def delete(self):
        '''Delete random oli from current_path'''
        if len(self.current_path) == 0: return None, None
        found = False
        index, to_delete = None, None
        for _ in range(10):
            index, to_delete = self.pick_random()
            flag1 = index != 0
            flag2 = index == len(self.current_path) - 1 or\
                self.isValid(self.current_path[index-1], self.current_path[index+1])
            if flag1 and flag2:
                found = True
                break
        if not found:
            return None, None

        def changes(confirm = False):
            if confirm:
                self.add_to_free(to_delete)
            self.current_path.pop(index)
        def unchanges():
            self.current_path.insert(index, to_delete)
        return changes, unchanges

    def add(self):
        '''Insert random free oli into random place of current_path'''
        if len(self.free) == 0: return None, None
        found = False
        index, chosen = None, None
        for i in range(10):
            if found: break
            index, _ = self.pick_random()
            chosen = self.choose_from_free()
            flag_previous = index == 0 or self.isValid(self.current_path[index], chosen)
            flag_next = index == len(self.current_path) - 1 or self.isValid(chosen, self.current_path[index+1])
            if flag_previous and flag_next:
                found = True
        if not found: return None, None

        def changes(confirm = False):
            if confirm:
                self.pop_from_free(chosen)
            self.current_path.insert(index+1, chosen)
        def unchanges():
            self.current_path.pop(index+1)
        return changes, unchanges

    def swap(self):
        '''Swap random element from current path with random unused neighbour'''
        if len(self.free) == 0: return None, None
        found = False
        index, to_swap, chosen = None, None, None
        before, after = None, None
        for _ in range(10):
            index, to_swap = self.pick_random()
            chosen = self.choose_neighbour(to_swap)
            if chosen is None: continue
            flag_one = index > 0 and self.isValid(self.current_path[index-1], chosen)
            flag_two = index == len(self.current_path) -1 or self.isValid(chosen, self.current_path[index+1])
            if flag_one and flag_two:
                found = True
                break
        if not found:
            return None, None
        def changes(confirm = False):
            self.current_path[index] = chosen
            if confirm:
                self.add_to_free(to_swap)
                self.pop_from_free(chosen)
        def unchanges():
            self.current_path[index] = to_swap
        return changes, unchanges


    def eval_state(self, state: tuple) -> float:
        '''Compute state's value'''
        odd, even = state
        a = abs(compute_length(odd) - self.LENGTH//2 - self.LENGTH%2)
        b = abs(compute_length(even) - self.LENGTH//2)
        c = 100*sum([1/max(overlap(odd[i], odd[i+1]), 1)**2 for i in range(len(odd) -1)])/len(odd)
        d = 100*sum([1/max(overlap(even[i], even[i+1]), 1)**2 for i in range(len(even) -1)])/len(even)
        #e = 10**6*len(odd) / sum([len(oli) for oli in odd])**2
        #f = 10**6*len(even) / sum([len(oli) for oli in even])**2
        return 10*(a+b+c+d)/4/self.LENGTH

    def valid_state(self, odd, even) -> Tuple[bool, bool]:
        a = compute_length(odd) == self.LENGTH//2 + self.LENGTH%2
        b = compute_length(even) == self.LENGTH//2
        return (a, b)

    def update_best_state(self, temp_val):
        up_odd, up_even = self.valid_state(self.odd, self.even)
        if up_odd:
            self.best_found[0] = self.odd[:]
        if up_even:
            self.best_found[1] = self.even[:]
        if up_odd or up_even:
            self.best_val = temp_val

    def fix_current_path(self):
        l = compute_length(self.current_path)
        i = 0
        if l < self.LENGTH//2:
            while l < self.LENGTH//2:
                self.mutate('add')
                i+=1
                if i == 100: break
                l = compute_length(self.current_path)
            return
        while l > self.LENGTH//2:
            self.mutate('del')
            i += 1
            if i == 100: break
            l = compute_length(self.current_path)


    #-----------------------------------------------------------------------------
    #--------------------------Helper functions for mutating----------------------
    #-----------------------------------------------------------------------------
    def pick_random(self) -> Tuple[int, str]:
        '''Get index and element from current path'''
        i = random.randrange(len(self.current_path))
        return i, self.current_path[i]

    def choose_neighbour(self, oli: str) -> Optional[str]:
        '''Choose unused neighbour within max_d distance from oli'''
        max_dist = 15#max(self.max_dist//2, math.ceil(self.T/self.MAX_T*max(self.dist[oli])))
        ks = [k for k in range(1, max_dist+1) if k in self.dist[oli]]
        possible = set([ind for k in ks for ind in self.dist[oli][k]]) & self.free
        if len(possible) == 0:
            return None
        chosen = self.S1[random.choice(list(possible))]   
        if chosen not in self.clustered_olis:
            return chosen
        return None

    def pop_from_free(self, oli: str = None) -> Optional[str]:
        '''Get random/chosen unused oli; mark it as used and remove it from free'''
        if oli is None:
            #random
            oli_i = self.free.pop()
            self.used_olis[oli_i] = True
            return self.S1[oli_i]
        #oli
        oli_i = self.oli2ind[oli]
        self.free.remove(oli_i)
        self.used_olis[oli_i] = True
        return None

    def choose_from_free(self) -> str:
        ret = self.free.pop()
        self.free.add(ret)
        return self.S1[ret]

    def add_to_free(self, oli: str) -> None:
        '''Mark oli as free, append i to free list'''
        oli_i = self.oli2ind[oli]
        self.free.add(oli_i)
        self.used_olis[oli_i] = False

    def isValid(self, o1: str, o2: str) -> bool:
        '''Check whether o2 can be placed after o1'''
        l1 = len(o1)
        i1, i2 = l1 - self.OLI_LENGTH + 2, self.OLI_LENGTH-1
        return o1[i1:] + o2[i2] in self.S2


    #-----------------------------------------------------------------------------
    #--------------------------------Clustering-----------------------------------
    #-----------------------------------------------------------------------------
    def try_to_cluster(self):
        index, mid = self.pick_random()
        if index == 0:
            return
        before = self.current_path[index-1]
        if self.can_cluster(before, mid):
            self.cluster(before, mid, index)

    def can_cluster(self, before, oli):
        '''Checks if clustering with 'before' oli and 'after' oli can be done'''
        OL = self.OLI_LENGTH
        b1, o1 = before[-OL:], oli[:OL]
        ov = overlap(b1, o1)
        return ov == OL//2

    def cluster(self, before, oli, i_oli):
        #i_oli - index of oli in current_path
        new_oli = before[:-self.OLI_LENGTH+2] + oli
        new_oli_i = len(self.S1)
        self.dist[new_oli] = self.compute_dist_single(new_oli, new_oli_i)
        self.used_olis[new_oli_i] = True
        self.S1.append(new_oli)
        self.oli2ind[new_oli] = new_oli_i

        for oligo in [before, oli]:
            self.dist.pop(oligo)
            oligo_i = self.oli2ind[oligo]
            self.used_olis.pop(oligo_i)
            self.clustered_olis.append(oligo)
            try: self.free.remove(oligo_i)
            except KeyError: pass
        self.current_path[i_oli] = new_oli
        self.current_path.remove(before)

    def compute_dist_single(self, oli: str, oli_i: int) -> Dict[int, List[int]]:
        d = {k:[] for k in range(1, self.max_dist+1)}
        for i, o1 in enumerate(self.S1):
            if oli != o1 and o1 not in self.clustered_olis:
                distance = dist(oli, o1)
                if distance not in d: d[distance] = []
                d[distance].append(i)
                if distance not in self.dist[o1]: self.dist[o1][distance] = []
                self.dist[o1][distance].append(oli_i)
        return d

    def construct_final_solution(self):
        odd, even = self.best_found
        odd_str, even_str = odd[0], even[0]
        for i in range(1, len(odd)):
            ov = max(2*overlap(odd[i-1], odd[i])-1, 0)
            if ov == 0: odd_str += 'X'
            odd_str += odd[i][ov:]
        for i in range(1, len(even)):
            ov = max(2*overlap(even[i-1], even[i])-1, 0)
            if ov == 0: even_str += 'X'
            even_str += even[i][ov:]
        result = ''
        for i, (o, e) in enumerate(zip(odd_str, even_str)):
            result += (o+e, '')[i%2]
        return result





def dist(o1: str, o2: str) -> int:
    l1, l2 = len(o1), len(o2)
    if l1 == l2:
        return sum([a!=b for a,b in zip(o1, o2)])
    if l2 < l1:
        o1, o2 = o2, o1     #make o1 shorter
        l1, l2 = l2, l1
    return min([sum([a!=b for a,b in zip(o1, o2[i:i+l1])]) for i in range(0, l2-l1+1, 2) ]) + l2 - l1

def compute_length(path):
    letters = len(path[0])//2 + 1
    for i in range(len(path) - 1):
        letters += len(path[i+1])//2 + 1 - overlap(path[i], path[i+1])
    return letters

def overlap(oli1, oli2):
    l1 = len(oli1)
    l2 = len(oli2)
    start = max(l1 - l2, 0)
    for i in range(start, l1, 2):
        if oli1[i:] == oli2[:l1-i]:
            a = (l1 - i)//2+1
            return a
    return 0
