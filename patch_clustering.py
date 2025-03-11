import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN


def get_positions_from_npz(paths):
    position_x = np.load(paths[0])
    position_y = np.load(paths[1])

    return position_x, position_y


paths_positions = ["generated/plants_last_x_positions.npz", "generated/plants_last_y_positions.npz"]
x, y = get_positions_from_npz(paths_positions)

points = np.column_stack((x, y))  # Fusionner x et y en une seule liste de points (N, 2)

# Appliquer DBSCAN
db = DBSCAN(eps=0.5, min_samples=2).fit(points)  # eps : distance max entre voisins
labels = db.labels_  # Les labels des clusters (-1 signifie bruit)

# Nombre de clusters
n_clusters = len(set(labels))

# Taille des clusters
unique_labels, counts = np.unique(labels, return_counts=True)
print(f"Nombre de clusters détectés : {n_clusters}")
print(f"Taille de chaque cluster : {counts}")

# Tracer le nuage de points avec les clusters
plt.figure()
plt.scatter(x, y, c=labels, cmap='viridis', edgecolor='k')
plt.colorbar(label="Cluster ID")
plt.title(f"Clustering DBSCAN - {n_clusters} clusters détectés")
plt.xlabel(r"$X$", size=14)
plt.ylabel(r"$Y$", size=14)
plt.savefig("figures/graphiques/Scénarios/2/patchs_clustering.png")

# Tracer l'histogramme du nombre de clusters en fonction de leur taille
# plt.figure()
# plt.hist(counts, edgecolor='k')
# plt.xlabel("Taille des patchs", size=14)
# plt.ylabel("Nombre de patchs", size=14)
# plt.savefig("figures/graphiques/Scénarios/histogramme_taille_clusters.png")
# plt.show()

# # Tracer la figure du nombre de cluster en fonction de leur taille
# plt.figure()
# plt.loglog(counts, unique_labels, 'o')
# plt.xlabel("Taille des patchs", size=14)
# plt.ylabel("Nombre de patchs", size=14)
# plt.savefig("figures/graphiques/Scénarios/power_law.png")
# plt.show()
