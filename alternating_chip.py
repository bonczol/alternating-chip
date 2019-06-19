#!/usr/bin/python3
import sys
from Sequencer import Sequencer
from DataLoader import DataLoader
from Randomer import Randomer
from Annealing import Annealing
import os



def main():
    URL = 'http://piotr.e.wawrzyniak.doctorate.put.poznan.pl/bio.php?n='+ sys.argv[1] + '&k=' + sys.argv[2] + '&mode=alternating&intensity=0&position=0&sqpe=' + sys.argv[3] + '&sqne=0&pose=0'
    print("URL:", URL)
    loader = DataLoader(URL)
    data = loader.getData()
    print("DNA length:" + sys.argv[1])
    print("oligo length:" + sys.argv[2])
    print("Starting sequence:" + data['start_seq'])

    if sys.argv[4] in 'GB':
        sequencer = Sequencer(data["set1"], data["set2"], data["len_DNA"], data["len_oligo"], data["start_seq"])
        result, _, _ = sequencer.run()
        print(result)
    if sys.argv[4] in 'HB':
        max_time = 10
        cooling_factor = 0.9999
        exponential_coef = 1.0
        clustering_factor = 20

        odd, even = Randomer().rand(data['set1'], data['start_seq'], data['len_DNA'])
        ann = Annealing(oli2seq(data['set1']),
                        data['set2'],
                        data['len_oligo'],
                        data['len_DNA'],
                        oli2seq(odd),
                        oli2seq(even),
                        max_time,
                        cooling_factor,
                        exponential_coef,
                        clustering_factor)
        res_odd, res_even, cluster_odd, cluster_even, best_val = ann.run()
        print("Ścieżka nieparzysta:")
        print(res_odd)
        print("Ścieżka parzysta:")
        print(res_even)

def oli2seq(ls):
    ret = []
    for oli in ls:
        if type(oli) == str:
            ret.append(oli)
        else:
            ret.append(oli.sequence)
    return ret


def test_quality_over_time(loader):
    cooling_factor = 0.9999
    exponential_coef = 1.0
    clustering_factor = 20
    for m_time in [10, 15, 20, 25, 30, 35, 40]:
        sum_odd, sum_even, sum_val = 0, 0, 0
        max_time = m_time
        for i in range(15):
            data = loader.parse_file('data_200/200_9_0.1.xml')
            odd, even = Randomer().rand(data['set1'], data['start_seq'], data['len_DNA'])
            ann = Annealing(oli2seq(data['set1']),
                            data['set2'],
                            data['len_oligo'],
                            data['len_DNA'],
                            oli2seq(odd),
                            oli2seq(even),
                            max_time,
                            cooling_factor,
                            exponential_coef,
                            clustering_factor)
            cluster_odd, cluster_even, best_val = ann.run()
            sum_odd += cluster_odd
            sum_even += cluster_even
            sum_val += best_val
        sum_val = sum_val/15
        sum_odd = sum_odd/15
        sum_even = sum_even/15
        with open('results/quality_time.txt', 'a') as f:
            f.write('{} {} {} {}\n'.format(max_time, sum_odd, sum_even, sum_val))
    
def test_quality_over_pe(loader):
    max_time = 25
    cooling_factor = 0.9999
    exponential_coef = 1.0
    clustering_factor = 20
    data_files = os.listdir('data_200')
    for file_name in data_files:
        sum_odd, sum_even, sum_val = 0, 0, 0
        for i in range(15):
            data = loader.parse_file('data_200/'+file_name)
            odd, even = Randomer().rand(data['set1'], data['start_seq'], data['len_DNA'])
            ann = Annealing(oli2seq(data['set1']),
                            data['set2'],
                            data['len_oligo'],
                            data['len_DNA'],
                            oli2seq(odd),
                            oli2seq(even),
                            max_time,
                            cooling_factor,
                            exponential_coef,
                            clustering_factor)
            cluster_odd, cluster_even, best_val = ann.run()
            sum_odd += cluster_odd
            sum_even += cluster_even
            sum_val += best_val
        sum_val = sum_val/15
        sum_odd = sum_odd/15
        sum_even = sum_even/15
        with open('results/quality_pe.txt', 'a') as f:
            f.write('{} {} {} {}\n'.format(file_name, sum_odd, sum_even, sum_val))

def test_quality_over_cooling(loader):
    exponential_coef = 1.0
    clustering_factor = 20
    max_time = 25
    for c in [0.99, 0.995, 0.999, 0.9995, 0.9999, 0.99995, 0.99999]:
        cooling_factor = c
        sum_odd, sum_even, sum_val = 0, 0, 0
        for i in range(15):
            data = loader.parse_file('data_200/200_9_0.1.xml')
            odd, even = Randomer().rand(data['set1'], data['start_seq'], data['len_DNA'])
            ann = Annealing(oli2seq(data['set1']),
                            data['set2'],
                            data['len_oligo'],
                            data['len_DNA'],
                            oli2seq(odd),
                            oli2seq(even),
                            max_time,
                            cooling_factor,
                            exponential_coef,
                            clustering_factor)
            cluster_odd, cluster_even, best_val = ann.run()
            sum_odd += cluster_odd
            sum_even += cluster_even
            sum_val += best_val
        sum_val = sum_val/15
        sum_odd = sum_odd/15
        sum_even = sum_even/15
        with open('results/quality_cooling.txt', 'a') as f:
            f.write('{} {} {} {}\n'.format(cooling_factor, sum_odd, sum_even, sum_val))

def test_quality_over_clustering(loader):
    files = ['data/400_9_0.1.xml', 'data/200_9_0.1.xml']
    exponential_coef = 1.0
    max_time = 25
    cooling_factor = 0.9999
    for file_name in files:
        for c in [5, 10, 15, 20, 25, 30, 35, 40]:
            clustering_factor = c
            sum_odd, sum_even, sum_val = 0, 0, 0
            for i in range(15):
                data = loader.parse_file(file_name)
                odd, even = Randomer().rand(data['set1'], data['start_seq'], data['len_DNA'])
                ann = Annealing(oli2seq(data['set1']),
                                data['set2'],
                                data['len_oligo'],
                                data['len_DNA'],
                                oli2seq(odd),
                                oli2seq(even),
                                max_time,
                                cooling_factor,
                                exponential_coef,
                                clustering_factor)
                cluster_odd, cluster_even, best_val = ann.run()
                sum_odd += cluster_odd
                sum_even += cluster_even
                sum_val += best_val
            sum_val = sum_val/15
            sum_odd = sum_odd/15
            sum_even = sum_even/15
            with open('results/quality_clustering.txt', 'a') as f:
                f.write('{} {} {} {} {}\n'.format(file_name, clustering_factor, sum_odd, sum_even, sum_val))

def test_quality_over_size(loader):
    files = ['data_var_size/' + f for f in os.listdir('data_var_size')]
    exponential_coef = 1.0
    max_time = 25
    cooling_factor = 0.9999
    clustering_factor = 20
    for file_name in files:
        sum_odd, sum_even, sum_val = 0, 0, 0
        for i in range(15):
            data = loader.parse_file(file_name)
            odd, even = Randomer().rand(data['set1'], data['start_seq'], data['len_DNA'])
            ann = Annealing(oli2seq(data['set1']),
                            data['set2'],
                            data['len_oligo'],
                            data['len_DNA'],
                            oli2seq(odd),
                            oli2seq(even),
                            max_time,
                            cooling_factor,
                            exponential_coef,
                            clustering_factor)
            cluster_odd, cluster_even, best_val = ann.run()
            sum_odd += cluster_odd
            sum_even += cluster_even
            sum_val += best_val
            print(file_name, i, cluster_odd, cluster_even, best_val)
        sum_val = sum_val/15
        sum_odd = sum_odd/15
        sum_even = sum_even/15
        with open('results/quality_size.txt', 'a') as f:
            f.write('{} {} {} {}\n'.format(file_name, sum_odd, sum_even, sum_val))

def test_quality_over_exponential(loader):
    clustering_factor = 20
    cooling_factor = 0.9999
    max_time = 25
    for e in [0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]:
        exponential_coef = e
        sum_odd, sum_even, sum_val = 0, 0, 0
        for i in range(15):
            data = loader.parse_file('data_200/200_9_0.1.xml')
            odd, even = Randomer().rand(data['set1'], data['start_seq'], data['len_DNA'])
            ann = Annealing(oli2seq(data['set1']),
                            data['set2'],
                            data['len_oligo'],
                            data['len_DNA'],
                            oli2seq(odd),
                            oli2seq(even),
                            max_time,
                            cooling_factor,
                            exponential_coef,
                            clustering_factor)
            cluster_odd, cluster_even, best_val = ann.run()
            sum_odd += cluster_odd
            sum_even += cluster_even
            sum_val += best_val
        sum_val = sum_val/15
        sum_odd = sum_odd/15
        sum_even = sum_even/15
        with open('results/quality_exponential.txt', 'a') as f:
            f.write('{} {} {} {}\n'.format(exponential_coef, sum_odd, sum_even, sum_val))



if __name__ == '__main__':
    main()

