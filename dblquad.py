import numpy as np
import matplotlib.pyplot as plt

def dblquad(x,y,f):
    X,Y = np.meshgrid(x,y)
    Y = np.tril(Y)
    Z = f(X,Y)
    int1 = np.trapz(Y,axis=1)
    int2 = np.trapz(int1,y)
    return int2

def triquad(x,y):
    m = -(y[-1]-y[0])/(x[-1]-x[0])
    b = y[-1]
    ytri = m*x+b
    # plt.plot(x,ytri)
    # plt.show()
    int1 = np.trapz(ytri,x)
    # int2 = np.trapz(np.ones_like(y),y)
    return int1

if __name__ == '__main__':
    f = 1
    x = np.linspace(0,4)
    y = np.linspace(0,3)
    print(triquad(x,y))