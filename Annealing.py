import random
import matplotlib.pyplot as plt
import math
from typing import List, Dict, Set, Tuple, Optional

class Annealing:
    def __init__(self, oli_set, S2, oli_length, DNA_length, odd, even):
        self.S1: List[str] = list(oli_set)                       #list
        self.S2: List[str] = list(S2)                            #zbiór oligosów potwierdzających
        self.LENGTH: int = DNA_length                            #docelowa długość łańucha DNA
        self.OLI_LENGTH: int = 2*oli_length-1                        #długość oligonukleotydu ze zbioru S1
        self.odd: List[str] = odd 
        self.even: List[str] = even

        self.NUM_OLIGO: int = self.LENGTH - self.OLI_LENGTH + 1        #docelowa liczba użytych oligonukleotydów
        self.max_dist: int = self.OLI_LENGTH//2+1                     #maksymalna odległość między oligonukleotydami

        self.used_olis = self.compute_olis_dict(self.S1)         #{oli: bool}
        self.free: Set[str] = self.compute_free()                #lista wolnych oligosów
        self.dist = self.compute_distances()                     #słownik dystansów
    
        self.current_path: List[str] = self.odd
        self.current_path_str: str = 'odd'
        self.T: float = 100                                      #temperatura
        self.ITER = 50000                                       #liczba iteracji

        self.current_val = self.eval_state((self.odd, self.even))

    #-----------------------------------------------------------------------------
    #--------------------------Init computations----------------------------------
    #-----------------------------------------------------------------------------
    def compute_olis_dict(self, oli_set: List[str]) -> Dict[str, bool]:
        used_olis = {oli: False for oli in oli_set}
        for oli in self.odd:
            used_olis[oli] = True
        for oli in self.even:
            used_olis[oli] = True
        return used_olis

    def compute_free(self) -> Set[str]:
        return set([oli for oli, used in self.used_olis.items() if not used])

    def compute_distances(self) -> Dict[str, Dict[int, List[int]]]:
        ret = {}    #dict of dicts: {oli_n: {k: [o1, o2, ...]}}
        for o1 in self.S1:
            d = {k:[] for k in range(1, self.max_dist+1)}
            for i, o2 in enumerate(self.S1):
                if o1 != o2: d[dist(o1, o2)].append(i)
            ret[o1] = d
        return ret 



    def run(self):
        xs, ys, y_len, y_oli = [], [], [], []
        for i in range(self.ITER):
            self.mutate()
            self.current_path_str = ('odd', 'even')[self.current_path_str == 'odd']
            self.current_path = [self.odd, self.even][self.current_path[0] == self.odd[0]]
            self.cool_down()
            xs.append(i)
            ys.append(self.current_val)
            y_len.append(5*abs(self.compute_state_length(self.odd, self.even) - self.LENGTH))
            y_oli.append(10*abs(len(self.odd) + len(self.even) - self.NUM_OLIGO))
            #print(i, self.current_val)
            if i % 100 == 0:
                print(self.compute_state_length(self.odd, self.even),
                    len(self.odd), len(self.even))

        plt.plot(xs, ys)
        plt.plot(xs, y_len)
        plt.plot(xs, y_oli)
        plt.show()

    #-----------------------------------------------------------------------------
    #--------------------------Temperature schedule and decision------------------
    #-----------------------------------------------------------------------------
    def cool_down(self):
        self.T -= 100/self.ITER
        
    def decide(self, val_old: int, val_new: int) -> bool:
        x = (val_new - val_old)/self.T
        try:
            return 1/(1 + math.exp(-0.4*x)) > random.random()
        except OverflowError:
            return True

    def compute_permitted_dist(self) -> int:
        return math.ceil(self.T/100*self.max_dist)
        

    #-----------------------------------------------------------------------------
    #--------------------------Mutating-------------------------------------------
    #-----------------------------------------------------------------------------
    def mutate(self) -> None:
        '''Mutate current path; decide whether accept this change'''
        action = random.choice([self.delete, self.add, self.swap])
        changes, unchanges = action()
        if changes is not None:
            changes()
            temp_val = self.eval_state((self.odd, self.even))
            unchanges()
            if self.decide(temp_val, self.current_val):
                changes(confirm = True)
                self.current_val = temp_val

    def delete(self) -> Tuple[List[str], List[str]]:
        '''Delete random oli from current_path'''
        found = False
        index, to_delete = None, None
        for _ in range(100):
            index, to_delete = self.pick_random()
            flag1 = index != 0
            flag2 = True#index == len(self.current_path) - 1 or\
                #self.isValid(self.current_path[index-1], self.current_path[index+1])
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

    def add(self) -> Tuple[List[str], List[str]]:
        '''Insert random free oli into random place of current_path'''
        found = False
        index, chosen = None, None
        for i in range(100):
            if found: break
            index, _ = self.pick_random()
            chosen = self.choose_from_free()
            flag_previous = index > 0   # and self.isValid(self.current_path[index-1], chosen) 
            flag_next = True   #len(self.current_path) - 1 or self.isValid(chosen, self.current_path[index+1])
            if flag_previous and flag_next:
                found = True
        if not found: return None, None

        def changes(confirm = False):
            if confirm:
                self.pop_from_free(chosen),
            self.current_path.insert(index, chosen)
        def unchanges():
            self.current_path.pop(index)

        return changes, unchanges

    def swap(self) -> Tuple[List[str], List[str]]:
        '''Swap random element from current path with random unused neighbour'''
        found = False
        index, to_swap, chosen = None, None, None
        before, after = None, None
        permitted_dist = self.compute_permitted_dist()
        for _ in range(100):
            index, to_swap = self.pick_random()
            chosen = self.choose_neighbour(to_swap, permitted_dist)
            if chosen is None: continue
            before = self.current_path[index-1]
            after = self.current_path[index+1]
            flag_one = index != 0# and self.isValid(before, chosen)
            flag_two = True#index == len(self.current_path) -1 or self.isValid(chosen, after)
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
                if index != len(self.current_path) -1 and self.can_cluster(before, chosen, after):
                    changes.append(lambda: self.cluster(before, chosen, after, index))
        def unchanges():
            self.current_path[index] = to_swap
        return changes, unchanges
    

    def eval_state(self, state: tuple) -> float:
        '''Compute state's value'''
        odd, even = state
        return 0.3*abs(self.compute_state_length(odd, even) - self.LENGTH) +\
                5*abs(len(odd) + len(even) - self.NUM_OLIGO)
                #3*(100 - self.T) *(len(odd) != len(even))

    def compute_state_length(self, odd, even):
        odd_letters = len(odd[0])//2 + 1
        even_letters = len(even[0])//2 + 1
        for i in range(len(odd) - 1):
            odd_letters += len(odd[i+1])//2 + 1 - self.overlap(odd[i], odd[i+1])
        for i in range(len(even) - 1):
            even_letters += len(even[i+1])//2 + 1 - self.overlap(even[i], even[i+1])
        return max(2*odd_letters-1, 2*even_letters) 

    def overlap(self, oli1, oli2):
        ret = 0
        l1 = len(oli1)
        l2 = len(oli2)
        start = max(l1 - l2, 0)
        for i in range(start, l1, 2):
            if oli1[i:] == oli2[:l1-i]:
                return (max(l1, l2) - i)//2+1
        return ret


    #-----------------------------------------------------------------------------
    #--------------------------Helper functions for mutating----------------------
    #-----------------------------------------------------------------------------

    def new_temp_state(self) -> Tuple[List[str], List[str]]:
        '''Get (odd, even) pair and set current_path'''
        if self.current_path_str == 'odd':
            odd = self.odd[:]
            even = self.even
        else:
            odd = self.odd
            even = self.even[:]
        return (odd, even)

    def pick_random(self) -> Tuple[int, str]:
        '''Get index and element from current path'''
        i = random.randrange(len(self.current_path))
        return i, self.current_path[i]

    def choose_neighbour(self, oli: str, max_d: int) -> Optional[str]:
        '''Choose unused neighbour within max_d distance from oli'''
        possible = set([ind for k in range(1, max_d+1) for ind in self.dist[oli][k] ]) & self.free
        if len(possible) == 0:
            return None
        return self.S1[random.choice(list(possible))]   
                                                    

    def pop_from_free(self, oli = None) -> Optional[str]:
        '''Get random/chosen unused oli; mark it as used and remove it from free'''
        if oli is None:
            #random
            oli = self.free.pop() 
            self.used_olis[oli] = True
            return oli
        #oli
        self.free.remove(oli)
        self.used_olis[oli] = True
        return None

    def choose_from_free(self):
        ret = self.free.pop()
        self.free.add(ret)
        return ret
        

    def add_to_free(self, oli: str) -> None:
        '''Mark oli as free, append i to free list'''
        self.free.add(oli)
        self.used_olis[oli] = False


    def isValid(self, o1: str, o2: str) -> bool:
        '''Check whether o2 can be placed after o1'''
        #print(o1, o2, o1[2:] + o2[-1], self.S2[0])
        if len(o1) == len(o2) ==  self.OLI_LENGTH:
            return o1[2:] + o2[-1] in self.S2
        return True


    #-----------------------------------------------------------------------------
    #--------------------------------Clustering-----------------------------------
    #-----------------------------------------------------------------------------
    def can_cluster(before, oli, after):
        '''Checks if clustering with 'before' oli and 'after' oli can be done'''
        ret = True
        lb, la = len(before), len(after)
        if lb > len(oli) or la > len(oli):
            return False
        for l1, l2 in zip(before[2:], oli[:lb-2]):
            ret = ret and l1 == l2
        for l1, l2 in zip(oli[-(la-2):], la[:-2]):
            ret = ret and l1 == l2
        return ret

    def cluster(self, before, oli, after, i_oli):
        #i_oli - index of oli in current_path
        #TODO: clustrowanie zmienia indeksy w S1 i wtedy S1 losowanie w swap robi sie useless
        new_oli = before[:2] + oli + after[-2:]
        self.dist[new_oli] = self.compute_dist_single(new_oli)
        self.used.olis[new_oli] = False
        self.free.add(new_oli)
        self.S1.append(new_oli)

        to_del = [self.S1.find(o) for o in [before, oli, after]]
        for oli, d in self.dist:
            for k, ind_list in d:
                removals = []
                for i in range(len(ind_list)):
                    if ind_list[i] in to_del:
                        removals.append(lambda: ind_list.pop(i))
                for rem in removals: rem()
                        
        for oligo in [before, oli, after]: delete_totally(oligo)

        def delete_totally(oligo: str):
            self.dist.pop(oligo)
            self.used_olis.pop(oligo)
            self.free.remove(oligo)
            self.S1.remove(oligo)

        self.curr_path.remove(before)
        self.curr_path[i_oli] = new_oli
        self.curr_path.remove(after)

    def compute_dist_single(oli: str) -> Dict[int, List[int]]:
        d = {k:[] for k in range(1, self.max_dist+1)}
        for i, o1 in enumerate(self.S1):
            if oli != o1:
                dist = dist(oli, o1)
                d[dist].append(i)
                self.dist[o1][dist].append(oli)
        return d

def dist(o1: str, o2: str) -> int:
    l1, l2 = len(o1), len(o2)
    if l1 == l2:
        return sum([a!=b for a,b in zip(o1, o2)])
    if l2 < l1:
        o1, o2 = o2, o1     #make o1 shorter
        l1, l2 = l2, l1
    return min([sum([a!=b for a,b in zip(o1, o2[i:i+l1])]) for i in range(0, l2-l1+1, 2) ]) + l2 - l1

