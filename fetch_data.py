#!/usr/bin/python3
import requests
#ns = [100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600]
ns = [500, 550, 600]
ks = [8, 9, 10]
ps = [0.05, 0.1, 0.15, 0.2, 0.25]
base_url='http://piotr.e.wawrzyniak.doctorate.put.poznan.pl/bio.php'

for n in ns:
    for k in ks:
        for pe in ps:
            sqpe = int(n*pe)
            url = '{}?n={}&k={}&mode=alternating&intensity=0&position=0&sqpe={}&sqne=0&pose=0'.format(base_url, n, k, sqpe)
            t, lt, attempts = '', 0, 0
            while lt < 100 and attempts < 10:
                t = requests.get(url).text
                lt = len(t)
                attempts += 1
            if attempts < 10:
                with open('data/{}_{}_{}.xml'.format(n, k, pe), 'w+') as f:
                    f.write(t)
                print("Written {} {} {}".format(n, k, pe))



