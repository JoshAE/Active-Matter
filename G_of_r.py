#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 10:42:58 2019

@author: joshua
"""
import numpy as np
import matplotlib.pyplot as plt
import skimage.measure as skim
from skimage.measure import label
from scipy.spatial import cKDTree
import scipy
from scipy.signal import find_peaks


def radbuild(r, r2, N=10000):
    radii = np.zeros(N)
    a = 0
    while a < N:
        radii[a] = r
        radii[a + 1] = r
        radii[a + 2] = r2
        radii[a + 3] = r2
        a = a + 4

    return radii


def fileopen(r2, Pe):
    with open("r2={}_Pe={}.txt".format(r2, Pe), "r") as f:
        contents = f.read()

    contents = contents.split()
    contents = np.array(contents, dtype='float64')

    x = contents[0::2]
    y = contents[1::2]

    return x, y


def system_size(r, r2, N=10000):
    size = np.sqrt((((N / 2) * np.pi) + ((N / 2) * np.pi * r2 * r2)) / 0.5)
    return size


def binsize(r, r2, N=10000, sizescale=0.5):
    size = system_size(r, r2, N)
    bz = int(size * sizescale)

    return bz


def periodic_bc_clusters(intensity, r, r2):
    bz = binsize(r, r2)

    if intensity[0, 0] != intensity[-1, -1] and intensity[0, 0] > 0 and intensity[0, 0] > 0:
        val = intensity[-1, -1]
        intensity[intensity == val] = intensity[0, 0]

    if intensity[-1, 0] != intensity[0, -1] and intensity[-1, 0] > 0 and intensity[0, -1] > 0:
        val = intensity[0, -1]
        intensity[intensity == val] = intensity[-1, 0]

    for i in range(int(bz)):

        if intensity[i, -1] != intensity[i, 0] and intensity[i, -1] > 0 and intensity[i, 0] > 0:
            val = intensity[i, -1]
            intensity[intensity == val] = intensity[i, 0]

        if intensity[-1, i] != intensity[0, i] and intensity[-1, i] > 0 and intensity[0, i] > 0:
            val = intensity[-1, i]
            intensity[intensity == val] = intensity[0, i]

    count = 0
    while count in range(np.max(intensity)):
        tf = count in intensity

        if tf == False:
            intensity[intensity > count] -= 1
        else:
            count = count + 1

    return intensity


def cluster_find(x, y, r, r2, N=10000, sizescale=0.5):
    bz = binsize(r, r2, N, sizescale)

    bins, xb, yb = np.histogram2d(x, y, bins=bz)
    stst, (xe, ye), binnum = scipy.stats.binned_statistic_dd((x, y), np.ones(len(x)),
                                                             bins=int(bz), statistic='count', expand_binnumbers=True)

    binw = bins

    bins[bins > 0] = 1

    intensity = label(bins.astype(int))

    intensity = periodic_bc_clusters(intensity, r, r2)

    return intensity, binw, binnum


def isolate_large(intensity, binw, binnum, N=10000):
    props = skim.regionprops(intensity, intensity_image=binw)

    areas = np.zeros(np.max(intensity))
    ii = np.zeros(np.max(intensity))
    for b in range(np.max(intensity)):
        ii[b] = props[b].mean_intensity

        areas[b] = props[b].area

    cluster = areas
    cluster[np.argwhere(areas < max(areas))] = 0
    cluster[np.argwhere(areas == max(areas))] = 1

    cluster_num = np.float(np.argwhere(areas == max(areas))) + 1

    intensity[intensity != cluster_num] = 0
    intensity[intensity == cluster_num] = 1

    points_in_cluster = [True if intensity[binnum[0, i] - 1, binnum[1, i] - 1] == 1
                         else False for i in range(N)]

    return points_in_cluster


def KNN_search_pbc(xnew, ynew, r, r2, N=10000, srch = 10000):
    size = system_size(r, r2, N)
    test = np.vstack((xnew, ynew)).T

    testtl = test
    testt = test
    testtr = test
    testr = test
    testbr = test
    testb = test
    testbl = test
    testl = test

    sizex = np.zeros((len(xnew), 2))
    sizey = np.zeros((len(ynew), 2))

    sizex[:, 0] = np.ones(len(xnew)) * size
    sizey[:, 1] = np.ones(len(ynew)) * size

    testtl = testtl - sizex
    testtl = testtl + sizey

    testt = testt + sizey

    testtr = testtr + sizex
    testtr = testtr + sizey

    testr = testr + sizex

    testbr = testbr + sizex
    testbr = testbr - sizey

    testb = testb - sizey

    testbl = testbl - sizex
    testbl = testbl - sizey

    testl = test - sizex

    testtotal = np.vstack((test, testtl))
    testtotal = np.vstack((testtotal, testt))
    testtotal = np.vstack((testtotal, testtr))
    testtotal = np.vstack((testtotal, testr))
    testtotal = np.vstack((testtotal, testbr))
    testtotal = np.vstack((testtotal, testb))
    testtotal = np.vstack((testtotal, testbl))
    testtotal = np.vstack((testtotal, testl))

    tree = cKDTree(testtotal)
    dist, indx = tree.query(test, k=srch)

    indx = indx % N

    return indx, dist


def gofr(x, y, r, r2, size, N=10000):
    indx, dist = KNN_search_pbc(x, y, r, r2, N)

    n, bins = np.histogram(dist[:, 1:], bins=1000)

    rad = bins[:-1] + (bins[1] - bins[0]) / 2
    divisor = 2 * np.pi * rad * (rad[1] - rad[0])
    n = n / divisor
    peaks, _ = find_peaks(n, height=1000)

    plt.plot(rad, n)
    plt.xlim([1, 20])

    return peaks, rad, n


def cluster_show(r, r2, Pe, N=10000):
    x, y = fileopen(r2, Pe)

    size = system_size(r, r2, N)
    radii = radbuild(r, r2, N)

    clusters, bin_width, bin_number = cluster_find(x, y, r, r2)

    points_in_cluster = isolate_large(clusters, bin_width, bin_number, N)

    plt.subplot(131)
    plt.scatter(x[points_in_cluster], y[points_in_cluster], radii[points_in_cluster])
    plt.xlim([0, size])
    plt.ylim([0, size])

    plt.subplot(132)
    plt.scatter(x, y, radii)

    plt.subplot(133)
    plt.imshow(clusters)


def big_gofr_plot(r2, Pe, N=10000, r=1):
    x, y = fileopen(r2, Pe)

    size = system_size(r, r2, N)
    radii = radbuild(r, r2, N)

    clusters, bin_width, bin_number = cluster_find(x, y, r, r2)

    points_in_cluster = isolate_large(clusters, bin_width, bin_number, N)

    xbig = np.hstack((x[0::4], x[1::4]))
    ybig = np.hstack((y[0::4], y[1::4]))
    big_points_in_cluster = np.hstack((points_in_cluster[0::4], points_in_cluster[1::4]))

    xbig = xbig[big_points_in_cluster]
    ybig = ybig[big_points_in_cluster]

    gofr(xbig, ybig, r, r2, size)



def main():
    small_radii = np.array([0.1, 0.2, 0.25, 0.3])
    Peclet = np.array([50, 20, 10, 5, 2, 1])

    plt.figure(1)
    r = 1
    Pe = Peclet[2]

    for r2 in small_radii:
        big_gofr_plot(r2, Pe)

        length = (4 * r * r) / (r + r2)
        l2 = 2 * np.sqrt((4 * r * r) - (0.5 * length * 0.5 * length))

        plt.plot(np.array([length, length]), np.array([0, 10000]),'k-')
        plt.plot(np.array([l2,l2]),np.array([0,10000]),'k')


    plt.legend(['r2 = 0.1', 'r2 = 0.2', 'r2 = 0.25', 'r2 = 0.3'])
    plt.title('Peclet = {}'.format(Pe))
    plt.show()


if __name__ == '__main__':
    main()
