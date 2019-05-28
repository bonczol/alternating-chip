import sys
from Sequencer import Sequencer
from DataLoader import DataLoader


def main():
    URL = 'http://piotr.e.wawrzyniak.doctorate.put.poznan.pl/bio.php?n='+ sys.argv[1] + '&k=' + sys.argv[2] + '&mode=alternating&intensity=0&position=0&sqpe=' + sys.argv[3] + '&sqne=0&pose=0'
    print("URL:", URL)
    loader = DataLoader(URL)
    data = loader.getData()
    sequencer = Sequencer(data["set1"], data["set2"], data["len_DNA"], data["len_oligo"], data["start_seq"])
    print("Result:", sequencer.run()[0])


if __name__ == '__main__':
    main()
