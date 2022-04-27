import numpy as np

def dblquad(x,y,f):
    X,Y = np.meshgrid(x,y)
    Z = f(X,Y)
    int1 = np.trapz(Z,x,axis=1)
    int2 = np.trapz(int1,y)
    return int2