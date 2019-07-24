#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 12:20:10 2019

@author: joshua
"""
import numpy as np
import matplotlib.pyplot as plt


r2 = 0.8
r = 1
N= 1000 
size = np.sqrt((((N / 2) * np.pi) + ((N / 2) * np.pi * r2 * r2)) / 0.5)
binsize = int(size*0.9)
h = np.zeros((binsize,binsize))




for i in range(9):

    num = i 
    
    with open("{}.txt".format(num), "r") as f:
        contents = f.read()

    contents = contents.split()
    contents = np.array(contents, dtype='float64')
    
    x = contents[0::2]
    y = contents[1::2]

    
    histarr = np.histogram2d(x, y, bins = binsize)
   
    
    h += histarr[0]




plt.figure(3)
plt.hist(h)



