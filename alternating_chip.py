import sys
from Sequencer import Sequencer
from DataLoader import DataLoader


def main():
	URL = 'http://piotr.e.wawrzyniak.doctorate.put.poznan.pl/bio.php?n=500&k=10&mode=alternating&intensity=0&position=0&sqpe=100&sqne=0&pose=0'
	loader = DataLoader(URL)
	data = loader.getData()
	print("set length:", len(data["set1"]))
	sequencer = Sequencer(data["set1"], data["set2"], data["len_DNA"], data["len_oligo"], data["start_seq"])

	#print(sequencer.run())


if __name__ == '__main__':
    main()