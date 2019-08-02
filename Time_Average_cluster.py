import numpy as np
from Cluster_find import system_size, square_distance, threshfind, radbuild


def time_cluster_labels(dataframe, r, r2, N=10000):
    size = system_size(r, r2, N)
    radii = radbuild(r, r2, N)

    dataframe['Cluster'] = 0

    for f in range(100999, 101999):

        x_data = dataframe.loc[dataframe['frame'] == f, 'x']
        y_data = dataframe.loc[dataframe['frame'] == f, 'y']

        points = np.vstack((x_data, y_data))
        points = np.transpose(points)

        square = square_distance(points, size)

        threshold = threshfind(size)

        cluster = DBSCAN_cluster(square, radii, thresh=threshold, min_samples=10,
                                 metric="precomputed", points=points,plot=False)

        dataframe.loc[dataframe['frame'] == f, 'Cluster'] = cluster

    return dataframe



N = 10000
r = 1
r2 = 0.2
Pe = 20

