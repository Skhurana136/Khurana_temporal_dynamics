# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 16:52:38 2020

@author: khurana
"""

import math
import numpy as np

# reading Tec files
# filename = 'model_domain_quad.tec'


def readTecfile(filename):
    numberOfHeaderLines = 3

    # pre-analysis to extract total number of lines in the data file which is
    # used further on to calculate number of time steps
    fid = open(filename, "r")

    # with open("t149_domain_quad.tec", "r") as fid
    tline_old = fid.readline()
    n = 1
    while len(tline_old) > 0:
        n = n + 1
        tline_old = fid.readline()

    # go to beginning of file;
    fid.seek(0)
    line_ex = (
        fid.readline()
    )  # check if readline includes or excludes new line character. chose function that removes new line character
    #    line_ex_a = line_ex.lstrip('VARIABLES  = "X","Y","Z",') #deleting this section of the variable string
    # Alternatively, slice string:
    line_ex_final = line_ex[12 : len(line_ex)].strip()
    Headers = line_ex_final.split(",")  # unpacking the variables
    #    Headers_new = [] REDUNDANT AFTER LINE 36
    #    for i in Headers:
    #        i = i.strip().replace('"','')
    #        Headers_new.append(i)

    line_ex = fid.readline()
    line_ex_final = line_ex.split(",")  # unpacking the variables
    fid.close()  # close the file or terminate the connection with the file system

    # Extracting number of data points
    line_ex_final = line_ex[line_ex.find("N=") :].partition(",")
    dataPoints = int(line_ex_final[0][2 : len(line_ex_final[0])])

    # Extracting useless data
    uslessDataLines = line_ex[line_ex.find("E=") :].partition(",")
    uslessDataLines = int(uslessDataLines[0][2 : len(uslessDataLines[0])])

    numberOfLinesInEachTimeStep = numberOfHeaderLines + dataPoints + uslessDataLines
    numberOfTimeSteps = int(math.floor(n / numberOfLinesInEachTimeStep))
    D = np.zeros(
        [dataPoints, len(Headers), numberOfTimeSteps]
    )  # Initializing np.darray

    # Reading select lines from the file
    count = 0
    for ii in range(numberOfTimeSteps):
        #        c = 0
        #        for i in range(750):
        #            line = linecache.getline(filename, (ii+1)*numberOfHeaderLines+ii*(uslessDataLines+dataPoints)+i+1)
        #            linex = line.split(' ')
        #            if (len(linex)!=len(Headers)+1):
        #                c = c + 1
        #        if (c > 0):
        #            print(ii, ":", c)
        #            count = count+1
        #            continue
        with open(filename) as fid:
            M = np.loadtxt(
                fid,
                dtype=str,
                delimiter=" ",
                skiprows=(ii + 1) * numberOfHeaderLines
                + ii * (uslessDataLines + dataPoints),
                max_rows=dataPoints,
            )
            if np.shape(M)[1] != len(Headers) + 1:
                count = count + 1
                continue
            elif np.shape(M)[0] != dataPoints:
                count = count + 1
                continue
        D[:, :, ii] = M[:, :-1]
    print("Total corrupted time steps: ", count)
    return dataPoints, numberOfTimeSteps, Headers, D

# converting loaded data to multi-dimensional array (timestepsx101rowsx61columns)
def Convertetoarray(D, datatype, vertnodes, horiznodes):
    if (datatype == "tec"):
        steps = np.shape(D)[2]
        Headers = np.shape(D)[1]
        size = np.shape(D)[0]
        df2 = np.ndarray([Headers - 3, steps, vertnodes, horiznodes])
        for i in range(Headers - 3):
            counter = 0
            for j in range(steps):
                #           print (Headers[i+3])
                while counter < (j + 1) * size:
                    for k in range(vertnodes):
                        for l in range(horiznodes):
                            df2[i, j, k, l] = D[counter - (j) * size, i + 3, j]
                            counter = counter + 1
        df = np.ndarray([Headers - 3, steps, vertnodes, horiznodes])
        for j in range(vertnodes):
            df[:, :, vertnodes - 1 - j, :] = df2[:, :, j, :]
    elif (datatype == "rates"):
        Headers = np.shape(D)[1]
        size = np.shape(D)[0]
        df2 = np.ndarray([Headers, vertnodes, horiznodes])
        for i in range(Headers):
            counter = 0
            while counter < size:
                for k in range(vertnodes):
                    for l in range(horiznodes):
                        df2[i, k, l] = D[counter - size, i]
                        counter = counter + 1
        df = np.ndarray([Headers, vertnodes, horiznodes])
        for j in range(vertnodes):
            df[:, vertnodes - 1 - j, :] = df2[:, j, :]        
    return df