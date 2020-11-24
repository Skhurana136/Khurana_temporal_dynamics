# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 16:54:58 2020

@author: khurana
"""
# Generating infiltration values based on precipitation, degree of saturation, porosity, depth of zone
def Qin(precip, ws, n, Zr):
    if n * Zr * (1 - ws) > precip:
        return precip
    elif n * Zr * (1 - ws) < precip:
        return n * Zr * (1 - ws)


# Generating recharge values based on precipitation, degree of saturation, porosity, depth of zone
def Qout(ws, wsmin, wsmax, Ks, c, ETmax, n, Zr, A):
    if ws <= wsmin:
        return 0
    elif ws >= wsmax:
        return Ks * (ws ** c) + ETmax
    elif ws < wsmax:
        return Ks * (ws ** c) + ws / (n * Zr * A)
