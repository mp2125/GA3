import ast
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler

# =========================
# GLOBAL STYLE CONTROL
# =========================
FONT_SIZE = 20

LABEL_FONT_SIZE = FONT_SIZE
TITLE_FONT_SIZE = FONT_SIZE + 2
TICK_FONT_SIZE = FONT_SIZE - 1
ARROW_LABEL_SIZE = FONT_SIZE

# NEW: arrow thickness scales with font size
ARROW_LW = max(1, FONT_SIZE / 6)   # keeps it reasonable for small/large fonts

data = []

with open("outputOptimisation.txt", "r") as f:
    for line in f:
        row = ast.literal_eval(line.strip())
        row = row[0:4] + row[5:]
        data.append(row)

data = np.array(data)

x = data[:, :6]
y1 = data[:, 6]
y2 = data[:, 7]

X_scaled = StandardScaler().fit_transform(data)

pca = PCA(n_components=2)
X_2d = pca.fit_transform(X_scaled)
loadings = pca.components_

fig, ax = plt.subplots(figsize=(8, 6))

sc = ax.scatter(
    X_2d[:, 0],
    X_2d[:, 1],
    c=y2,
    cmap="viridis",
    alpha=0.6
)

cbar = plt.colorbar(sc, ax=ax)
cbar.set_label("Q_eNTU", fontsize=LABEL_FONT_SIZE)
cbar.ax.tick_params(labelsize=TICK_FONT_SIZE)

# Scale arrows to ~30% of the data range on each axis
x_scale = (X_2d[:, 0].max() - X_2d[:, 0].min()) * 0.3
y_scale = (X_2d[:, 1].max() - X_2d[:, 1].min()) * 0.3

textPositions = [
    (1.1, 0.1), # tubes
    (6.5, 1.3), # baffles
    (2.5, 1.1), # lengths
    (1.6, 1.7), # pitches
    (1.0, 0.8), # tube passes
    (0.9, 1.3), # shell passes
]

for i, name in enumerate([
    'tubes','baffles','lengths','pitches','tube passes','shell passes'
]):
    dx = loadings[0, i] * x_scale
    dy = loadings[1, i] * y_scale

    ax.annotate(
        "",
        xy=(dx, dy),
        xytext=(0, 0),
        arrowprops=dict(
            arrowstyle="->",
            color="red",
            lw=ARROW_LW   # <-- scaled arrow thickness
        )
    )

    ax.text(
        dx * textPositions[i][0],
        dy * textPositions[i][1],
        name,
        color="red",
        fontsize=ARROW_LABEL_SIZE,
        fontweight="bold"
    )

ax.axhline(0, color="grey", lw=0.5, linestyle="--")
ax.axvline(0, color="grey", lw=0.5, linestyle="--")

ax.set_xlabel("PC1", fontsize=LABEL_FONT_SIZE)
ax.set_ylabel("PC2", fontsize=LABEL_FONT_SIZE)
# ax.set_title(
#     "Biplot: parameter directions in PCA space",
#     fontsize=TITLE_FONT_SIZE
# )

ax.tick_params(axis='both', labelsize=TICK_FONT_SIZE)

plt.tight_layout()
plt.savefig('PCA_plot.pdf')