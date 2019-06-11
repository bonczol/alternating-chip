#!/usr/bin/python3
import sys
from Sequencer import Sequencer
from DataLoader import DataLoader
from Randomer import Randomer
from Annealing import Annealing


def main():
    URL = 'http://piotr.e.wawrzyniak.doctorate.put.poznan.pl/bio.php?n='+ sys.argv[1] + '&k=' + sys.argv[2] + '&mode=alternating&intensity=0&position=0&sqpe=' + sys.argv[3] + '&sqne=0&pose=0'
    print("URL:", URL)
    loader = DataLoader(URL)
    data = loader.getData()
    sequencer = Sequencer(data["set1"], data["set2"], data["len_DNA"], data["len_oligo"], data["start_seq"])
    res, odd, even = sequencer.run()
    odd, even = Randomer().rand(data['set1'], data['start_seq'], data['len_DNA'])
    #print(rand)
    def oli2seq(ls):
        ret = []
        for oli in ls:
            if type(oli) == str:
                ret.append(oli)
            else:
                ret.append(oli.sequence)
        return ret
    ann = Annealing(oli2seq(data['set1']),
                    data['set2'],
                    data['len_oligo'],
                    data['len_DNA'],
                    oli2seq(odd),
                    oli2seq(even))
    print("VWEVE")
    ann.run()


if __name__ == '__main__':
    main()

