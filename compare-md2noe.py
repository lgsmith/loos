import json
import numpy as np
import matplotlib.pyplot as plt

fn = 'sucrose-1.json'
with open(fn, 'r') as fp:
    d = json.load(fp)

glucose_H1 = '0GA2H1'
fructose_H11 = '2CU1H11'
fructose_H12 = '2CU1H12'

noes = d['noes']
for key in noes.keys():
    if glucose_H1 in key:
        if fructose_H12 in key or fructose_H11 in key:
            print(key)
            noe=np.array(list(map(float, noes[key])))
            plt.plot(d['mixing times'], noe)
            plt.show()