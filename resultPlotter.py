import ast
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

data = []

with open("outputOptimisation.txt", "r") as f:
    for line in f:
        row = ast.literal_eval(line.strip())  # converts string → list
        row = row[0:4] + row[5:]
        data.append(row)

data = np.array(data)  # shape: (n_samples, 7)

x = data[:, :6]
y1 = data[:, 6]
y2 = data[:, 7]

from sklearn.preprocessing import StandardScaler

X_scaled = StandardScaler().fit_transform(data)

pca = PCA(n_components=2)
X_2d = pca.fit_transform(X_scaled)
loadings = pca.components_

fig, ax = plt.subplots(figsize=(8, 6))

sc = ax.scatter(X_2d[:, 0], X_2d[:, 1], c=y2, cmap="viridis", alpha=0.6)
plt.colorbar(sc, ax=ax, label="Q_eNTU")

# Scale arrows to ~30% of the data range on each axis
x_scale = (X_2d[:, 0].max() - X_2d[:, 0].min()) * 0.3
y_scale = (X_2d[:, 1].max() - X_2d[:, 1].min()) * 0.3

textPositions = [(1.1, 1.1),
                 (2.5, 1.1),
                 (1.5, 1.1),
                 (1.3, 1.5),
                 (1.0, 1.1),
                 (0.9, 1.3)]

for i, name in enumerate(['tubes','baffles','lengths','pitches','tube passes','shell passes']):
    dx = loadings[0, i] * x_scale
    dy = loadings[1, i] * y_scale

    ax.annotate("", xy=(dx, dy), xytext=(0, 0),
                arrowprops=dict(arrowstyle="->", color="red", lw=2))
    ax.text(dx * textPositions[i][0], dy * textPositions[i][1], name, color="red", fontsize=10, fontweight="bold")

ax.axhline(0, color="grey", lw=0.5, linestyle="--")
ax.axvline(0, color="grey", lw=0.5, linestyle="--")
ax.set_xlabel("PC1")
ax.set_ylabel("PC2")
ax.set_title("Biplot: parameter directions in PCA space")
plt.tight_layout()

plt.savefig('PCA_plot.pdf')
plt.show()

