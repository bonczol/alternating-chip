import random

class Randomer:
    def __init__(self):
        pass

    def rand(self, S1, start, length):
        '''
        S1 - set of oligonucleotides
        start - start od DNA sequence
        length - length of DNA sequence
        '''
        start_odd = ''.join([(start[i], 'X')[i%2] for i in range(len(start))])
        start_even = ''.join([('X', start[i])[i%2] for i in range(1, len(start))])
        
        odd, even = [start_odd], []
        visited = [start_odd]
        current = odd
        for oli in S1:
            s = oli.sequence
            if s[:-1] == start_even:
                even.append(s)
                visited.append(s)
        while len(visited) < length:     #length//2 is arbitrary
            temp = random.choice(S1)
            if temp not in visited:
                visited.append(temp)
                current.append(temp)
                current = (odd, even)[current[0] == odd[0]]
        return (odd, even)



