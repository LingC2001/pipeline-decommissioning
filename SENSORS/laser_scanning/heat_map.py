import matplotlib.pyplot as plt
import numpy as np
import random

heatmap = np.random.uniform(0.1, 0.2, (100, 100))

for i in range(heatmap.shape[0]):
    for j in range(heatmap.shape[1]):
        if np.sqrt((i-60)**2 + (j-50)**2) <= 35:
            heatmap[i, j] = random.uniform(0.9, 1)
        if np.sqrt((i-60)**2 + (j-50)**2) <= 20:
            heatmap[i, j] = random.uniform(0.4, 0.5)
        if np.round(np.sqrt((i-60)**2 + (j-50)**2)) == 32:
            heatmap[i, j] = 0.8
        if np.round(np.sqrt((i-60)**2 + (j-50)**2)) == 24:
            heatmap[i, j] = 0.8

plt.imshow(heatmap, cmap='bwr', interpolation='nearest')
plt.show()

