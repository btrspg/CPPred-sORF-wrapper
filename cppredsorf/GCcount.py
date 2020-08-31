#!/bin/python
# -*- coding:utf-8 -*- 
#########################################################################
# File Name: GCcount.py
# Author: Tong xiaoxue
# mail: xiaoxuetong@hust.edu.cn
# Created Time: Thu 14 Feb 2019 10:56:54 PM CST
#########################################################################
import sys

import numpy as np


def main():
    pattern = '>'

    with open(sys.argv[1]) as f:
        data = f.readlines()
        i = 0
    for line in data:
        if pattern in line:
            data.pop(i)
        i += 1

    wordsArr1 = []
    for i in range(0, len(data)):
        wordsArr1.append(data[i].strip('\n'))

    #	print wordsArr1
    gcCount = []

    for i in range(0, len(wordsArr1)):
        gcCount.append(wordsArr1[i].count('GC'))

    print("GC count: ", np.asarray(gcCount))


if __name__ == '__main__':
    main()
