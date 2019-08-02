import numpy as np
from James_functions import add_param, add_cluster
from G_of_r import gofr, system_size, radbuild, KNN_search_pbc
import pandas as pd
import scipy.spatial as sp
from sklearn.cluster import AgglomerativeClustering, DBSCAN
import matplotlib.pyplot as plt
from sklearn import metrics
from sklearn.datasets.samples_generator import make_blobs
from sklearn.preprocessing import StandardScaler
import pylab as pl
from sklearn.cluster import DBSCAN
from scipy.spatial.distance import pdist, squareform
import matplotlib.cm as cm

def df_build(r, r2, Pe, N=10000, mac=True):
    if mac == False:
        x = np.array([])
        y = np.array([])
        frame = np.array([])

        for i in range(100999, 101999):
            with open("r2={}_Pe={}/{}.txt".format(r2, Pe, i), "r") as f:
                contents = f.read()

                contents = contents.split()
                contents = np.array(contents, dtype='float64')

            x = np.append(x, contents[0::2])
            y = np.append(y, contents[1::2])
            frame = np.append(frame, [i - 100999] * len(x))

        df = pd.DataFrame(x=x, y=y, frame=frame)

    else:

        with open("r2={}_Pe={}.txt".format(r2, Pe), "r") as f:
            contents = f.read()

            contents = contents.split()
            contents = np.array(contents, dtype='float64')

        x = contents[0::2]
        y = contents[1::2]
        points = np.vstack((x, y)).T
        frame = [0] * len(x)

    d = {'x': x, 'y': y, 'frame': frame}
    df = pd.DataFrame(d)

    return df, points


def Ward_cluster(points):

    clustering = AgglomerativeClustering(n_clusters=None, distance_threshold=2.5).fit(points)

    josh = clustering.labels_

    radii = radbuild(r, r2, N)
    plt.scatter(points[:, 0], points[:, 1], s=radii, c=josh)
    plt.show()


def DBSCAN_cluster(X, radii, thresh=2.1, min_samples=10, metric="euclidean", points=None, plot=True):

    db = DBSCAN(eps=thresh, min_samples=min_samples, metric=metric).fit(X)

    if metric != "euclidean":
        X = points
    else:
        pass

    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    labels = db.labels_

    unique_labels = set(labels)
    colors = [plt.cm.Spectral(each)
              for each in np.linspace(0, 1, len(unique_labels))]

    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise_ = list(labels).count(-1)

    if plot==True:
        for k, col in zip(unique_labels, colors):
            if k == -1:
                # Black used for noise.
                col = [0, 0, 0, 1]

            class_member_mask = (labels == k)

            xy = X[class_member_mask & core_samples_mask]
            rad1 = radii[class_member_mask & core_samples_mask]

            plt.scatter(xy[:, 0], xy[:, 1], s=rad1 * 5, c=tuple(col))

            xy = X[class_member_mask & ~core_samples_mask]
            rad2 = radii[class_member_mask & ~core_samples_mask]
            plt.scatter(xy[:, 0], xy[:, 1], s=rad2 * 5, c=tuple(col))

        plt.title('Estimated number of clusters: %d' % n_clusters_)
        plt.show()

    return labels


def square_distance(points, size):
    # 1) find the correct distance matrix
    for d in range(points.shape[1]):
        # find all 1-d distances
        pd = pdist(points[:, d].reshape(points.shape[0], 1))
        # apply boundary conditions
        pd[pd > size * 0.5] -= size

        try:

            # sum
            total += pd ** 2

        except Exception as e:

            # or define the sum if not previously defined
            total = pd ** 2
    # transform the condensed distance matrix...
    total = pl.sqrt(total)
    # ...into a square distance matrix
    square = squareform(total)

    return square

def threshfind(size):
    sf = 2.2 / (system_size(1, 0) ** 2)
    threshold = size * size * sf

    return threshold


#
# N = 10000
# r2 = 0.6
# r = 1
#
# size = system_size(r, r2, N)
# radii = radbuild(r, r2, N)
# Peclet = np.array([50, 20, 10, 5, 2, 1])
#
# Pe = Peclet[2]
#
# threshold = threshfind(size)
#
# data_frame, points = df_build(r, r2, Pe, N)
#
# data_frame2 = add_param(data_frame, r, r2, N)
#
# plt.figure(1)
# lab = DBSCAN_cluster(points, radii, thresh=threshold, min_samples=10)
#
# square = square_distance(points, size)
#
# plt.figure(2)
# lab2 = DBSCAN_cluster(square, radii, thresh=threshold, min_samples=10,
#                metric="precomputed", points=points)
#
# plt.figure(3)
# coloring = lab2[np.argwhere(lab2 >= 0)]
#
#
# plt.scatter(points[np.where(lab2 >=0),0],points[np.argwhere(lab2>=0),1],
#             s = radii, c = coloring)
#
#

