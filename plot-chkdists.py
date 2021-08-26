import numpy as np
import matplotlib.pyplot as plt
glucose_H1 = '0GA2H1'
fructose_H11 = '2CU1H11'
fructose_H12 = '2CU1H12'
data = np.loadtxt('chkdists.out')[::100]

fig, ax = plt.subplots()
ax.plot(data[:, 0], data[:, 1], label=glucose_H1+fructose_H11)
ax.plot(data[:, 0], data[:, 2], label=glucose_H1+fructose_H12)
ax.set_ylabel('angstroms')
ax.set_xlabel('frame no.')
fig.suptitle('Distance across glycosidic linkage, sample traj MD2NOE')
ax.legend()
fig.savefig('chkdists.pdf', transparent=True)