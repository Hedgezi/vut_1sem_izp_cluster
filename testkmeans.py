import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import sys

x = []
y = []

with open('objekty', 'r') as f:
    f.readline()
    for i in f.readlines():
        x.append(int(i.split(' ')[1]))
        y.append(-(int(i.split(' ')[2])))

numtogen = int(sys.argv[1])

kmeans = KMeans(n_clusters=numtogen)
kmeans.fit(list(zip(x,y)))

plt.scatter(x, y, c=kmeans.labels_)
plt.show()