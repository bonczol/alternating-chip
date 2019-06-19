from Oligo import Oligo
import requests
from xml.etree import ElementTree


class DataLoader:
	def __init__(self, URL):
		self.URL = URL;


	def sendRequest(self):
		 return requests.get(self.URL)
	
	
	def getData(self):
		response = self.sendRequest()
		root = ElementTree.fromstring(response.content)

		data = dict()
		data["len_DNA"] = int(root.attrib["length"])
		data["start_seq"] = root.attrib["start"]
		data["len_oligo"] = int((len(data["start_seq"]) + 1)/ 2)
		data["set1"] = [ Oligo(child.text) for child in root[0]]
		data["set2"] = [ child.text for child in root[1]]

		return data

	def parse_file(self, file_name):
		with open(file_name) as file:
			root = ElementTree.fromstring(file.read())
			data = dict()
			data["len_DNA"] = int(root.attrib["length"])
			data["start_seq"] = root.attrib["start"]
			data["len_oligo"] = int((len(data["start_seq"]) + 1) / 2)
			data["set1"] = [Oligo(child.text) for child in root[0]]
			data["set2"] = [child.text for child in root[1]]
		return data



