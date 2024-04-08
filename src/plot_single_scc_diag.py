import matplotlib
import matplotlib.pyplot as plt
import numpy as np


labels = [ '262x262', '300x300', '516x516', '524x524', '1032x1032']
l=[0,1,2,3,4]
real_ratio=[50, 45, 32, 40, 37]
cplx_ratio=[45, 48, 35, 32, 27]
x = np.arange(len(labels))  # the label locations

fig, ax = plt.subplots()

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('SYMM/DC time ratio (%) ')

ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.axhline(25)
#ax.plot(labels, real_ratio, 'bs', label="gamma point")
ax.plot(l, real_ratio, 'bs', label="gamma point")
ax.plot(l, cplx_ratio, 'g^', label=" high symm k point")

ax.legend()

fig.tight_layout()

plt.show()
