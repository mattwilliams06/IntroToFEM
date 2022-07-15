import numpy as np
import os
import matplotlib.pyplot as plt

def UniformMeshQuad4(sx,sy,ex,ey,nx,ny,return_nodes=False):
    dx = (ex-sx)/nx
    dy = (ey-sy)/ny
    x = np.arange(sx,ex+dx,dx)
    y = np.arange(sy,ey+dy,dy)
    X,Y = np.meshgrid(x,y)

    nodes = np.zeros(((nx+1)*(ny+1),3),dtype='object')
    elements = np.zeros((nx*ny,5),dtype='int')
    nids = np.zeros((ny+1,nx+1))

    k = 0
    for i in range(ny+1):
        for j in range(nx+1):
            nids[i,j] = k
            nodes[k,:] = int(k), X[i,j], Y[i,j]
            k += 1
    with open('nodes.txt','w') as f:
        f.write('\n'.join(' '.join(map(str,x)) for x in nodes))

    k = 0
    for i in range(ny):
        for j in range(nx):
            elements[k,:] = int(k), nids[i,j], nids[i,j+1], nids[i+1,j+1], nids[i+1,j]
            k += 1
    with open('elements.txt','w') as f:
        f.write('\n'.join(' '.join(map(str,x)) for x in elements))

    if return_nodes:
        return nodes, elements
    else:
        return None
    

if __name__ == '__main__':
    path = '/Users/mattjwilliams/Documents/PythonStuff/FEM/GangLi/Chapter6/'
    os.chdir(path)
    nodes, elements = UniformMeshQuad4(0,0,3,5,10,20,True)
    for i in range(len(elements)):
        nodez = elements[i][1:]
        n = len(nodez)
        for j in range(n):
            x1,y1 = nodes[nodez[j%n]][1:]
            x2,y2 = nodes[nodez[(j+1)%n]][1:]
            plt.plot([x1,x2],[y1,y2],'k')
    plt.show()