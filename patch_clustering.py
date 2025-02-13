import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN

# Exemple : Génération d'un nuage de points aléatoire
np.random.seed(42)
x = np.random.rand(100) * 10  # Positions x
y = np.random.rand(100) * 10  # Positions y
points = np.column_stack((x, y))  # Fusionner x et y en une seule liste de points (N, 2)

# Appliquer DBSCAN
db = DBSCAN(eps=0.8, min_samples=4).fit(points)  # eps : distance max entre voisins
labels = db.labels_  # Les labels des clusters (-1 signifie bruit)

# Nombre de clusters
n_clusters = len(set(labels))

# Taille des clusters
unique_labels, counts = np.unique(labels, return_counts=True)
print(unique_labels)
print(counts)
cluster_sizes = dict(zip(unique_labels, counts))

print(f"Nombre de clusters détectés : {n_clusters}")
print(f"Taille de chaque cluster : {cluster_sizes}")

# Tracer le nuage de points avec les clusters
plt.scatter(x, y, c=labels, cmap='viridis', edgecolor='k')
plt.colorbar(label="Cluster ID")
plt.title(f"Clustering DBSCAN - {n_clusters} clusters détectés")
plt.xlabel("X")
plt.ylabel("Y")
plt.show()
