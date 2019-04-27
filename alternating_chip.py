import requests
from xml.etree import ElementTree as ET


def main():
    response = requests.get('http://piotr.e.wawrzyniak.doctorate.put.poznan.pl/bio.php?n=500&k=10&mode=alternating&intensity=0&position=0&sqpe=100&sqne=0&pose=0')
    tree = ET.parse(response.content)


if __name__ == '__main__':
    main()