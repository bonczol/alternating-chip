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
		data["start_seq"] = root.attrib["start"] + "A"
		data["len_oligo"] = int(len(data["start_seq"])/ 2)
		data["set1"] = [ Oligo(child.text) for child in root[0]]
		data["set2"] = [ Oligo(child.text) for child in root[1]]

		return data


