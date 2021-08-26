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
fig, ax = plt.subplots()
mixings = np.array(list(map(float, d['mixing times'])))
selected_buildups = np.zeros_like(mixings)
for key in noes.keys():
    if glucose_H1 in key:
        if fructose_H12 in key or fructose_H11 in key:
            print(key)
            noe = np.array(list(map(float, noes[key])))
            selected_buildups += noe
            ax.plot(mixings, noe, label=key)

ax.plot(mixings, selected_buildups/2, label='sum')
ax.legend()
ax.set_xlabel('mixing time, milliseconds')
ax.set_ylabel('NOE intensity')
fig.suptitle('Build-up curve, my tool, S.A. Chalmers et al. Fig. 5')
fig.savefig('buildup-fig5.pdf', transparent=True)
plt.show()