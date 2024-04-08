import matplotlib
import matplotlib.pyplot as plt
import numpy as np


labels = ['262x262', '300x300', '516x516', '524x524', '1032x1032']
dc_totals = [212, 251.3, 7828.8, 4796.7, 10425.8]
symm_totals = [95.7, 121.9, 2749.9, 1490.9, 2878.3]

x = np.arange(len(labels))  # the label locations
width = 0.35  # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(x - width/2, dc_totals, width, label='DivideAndConquer')
rects2 = ax.bar(x + width/2, symm_totals, width, label='Symmetry adapted')

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('SCC Times (ms)')

ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.legend()


def autolabel(rects):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')


autolabel(rects1)
autolabel(rects2)

fig.tight_layout()

plt.show()
