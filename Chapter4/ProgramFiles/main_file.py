import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt

class Lagrange:
    ''' A class for building quadratic Lagrange polynomials.
    
    Inputs
    ------
    x: array-like
    y: array-like

    Returns
    -------
    u: The Polynomial interpolation over the range covered
    by x and y
    '''
    def __init__(self):
        pass

    def _build_basis(self,x,idx):
        xarr = np.linspace(x[0],x[-1])
        temp = np.ones_like(xarr)
        for i in range(len(x)):
            if i != idx:
                basis = (xarr-x[i])/(x[idx]-x[i])
                temp *= basis
        return temp, xarr

    def __call__(self,x,y,return_x=False):
        self.Narr = []
        for i in range(len(x)):
            N, xarr = self._build_basis(x,i)
            self.Narr.append(N)
        u = np.zeros_like(xarr)
        for i, N in enumerate(self.Narr):
            u += N*y[i]
        if return_x:
            return u, xarr
        else:
            return u

# def create_1D_mesh(xstart,xend,n_elements):
#     dx = (xend-xstart)/n_elements
#     nodes = np.arange(xstart,xend+dx,dx)

def assemble_global_stiffness(K,F,k,f,eid,elements):
    gidx = [elements.iloc[eid]['node1'],elements.iloc[eid]['node2']]
    print(f'GIDX: {gidx}')
    for i,row in enumerate(gidx):
        F[row] += f[i]
        for j,col in enumerate(gidx):
            K[row,col] += k[i,j]
    return K, F

def lin_shape_iso(xi_vector):
    ''' Computes the 1D linear isoparametric shape functions and their derivatives.
    See page 118 - 120 in Li for more information.

    Inputs
    ------
    xi: array-like. A list of points [xi] in the master element between -1 and 1 
    where the shape functions should be evaluated.
    
    Returns
    -------
    N: matrix.  A matrix where the first column is N1(xi) and column 2 is N2(xi)
    Nx: matrix.  A matrix of the derivateves of N1 and N2, stored column-wise.
    '''
    n = len(xi_vector)
    N = np.zeros((2,n))
    Nx = np.zeros((2,n))
    N[0,:] = (1-xi_vector)/2
    N[1,:] = (1+xi_vector)/2
    Nx[0,:] = -0.5
    Nx[1,:] = 0.5
    return N, Nx

def jacobian(x_vector,Nx):
    ''' Computes the Jacobian dx/dxi, which is used in the derivative approximations
    using isoparametric mappings.  See page 121 of Li.

    Inputs
    ------
    x_vector: array.  Vector containing the x-coordinates of the nodes
    Nx: array.  Vectir containing the derivatives computed in lin_shape_iso
    '''
    J = x_vector@Nx
    return J

def get_1d_gauss(n_points):

    if n_points == 1:
        points = 0.0
        weights = 2.0
    elif n_points == 2:
        points = [-1/np.sqrt(3),1/np.sqrt(3)]
        weights = [1.0,1.0]
    elif n_points == 3:
        points = [-np.sqrt(3/5), 0.0, np.sqrt(3/5)]
        weights = [5/9,8/9,5/9]
    elif n_points == 4:
        points = [-0.861136, -0.339981, 0.339981, 0.861136]
        weights = [0.34785, 0.652145, 0.652145, 0.34785]
    else:
        raise ValueError('Error calling get_1D_gauss')
    return np.array(points), np.array(weights)


if __name__ == '__main__':
    # Load files into dataframes
    laptop = True
    if laptop:
        fpath = '/Users/mattwilliams/Documents/PythonProjects/FEM/GangLi/IntroToFEM/Chapter4/ProgramFiles'
    else:
        fpath = '/Users/mattjwilliams/Documents/PythonStuff/FEM/GangLi/Chapter4/ProgramFiles'
    fnames = ['elements.csv','ndprops.csv','nodes.csv','materials.txt','bfs.csv','bcs.csv']
    df_elements = pd.read_csv(os.path.join(fpath,fnames[0]),delimiter=' ')
    nelems = len(df_elements)
    df_ndprops = pd.read_csv(os.path.join(fpath,fnames[1]),delimiter=' ')
    df_nodes = pd.read_csv(os.path.join(fpath,fnames[2]),delimiter=' ')
    df_bfs = pd.read_csv(os.path.join(fpath,fnames[4]),delimiter=' ')
    df_bcs = pd.read_csv(os.path.join(fpath,fnames[5]),delimiter=' ')
    with open(os.path.join(fpath,fnames[3]),'r') as f:
        data = f.readlines()[1]
        E = float(data)

    # Create empty global matrices for later use
    K = np.zeros((nelems+1,nelems+1))
    F = np.zeros(nelems+1)
    # print(df_elements)
    # print(df_ndprops)
    # print(df_nodes)
    # print(E)

    # We're using linear nodes, so 2 nodes per element
    nodes_per_elem = 2
    n_nodes = len(df_nodes)
    gauss_xi, gauss_w = get_1d_gauss(nodes_per_elem)
    area = np.zeros(nodes_per_elem)
    sv = np.zeros(nodes_per_elem)
    for i in range(nelems):
        node1 = df_elements.iloc[i]['node1']
        node2 = df_elements.iloc[i]['node2']
        start_x = df_nodes.iloc[node1]['x-location']
        end_x = df_nodes.iloc[node2]['x-location']
        # print(f'Element {i}, x-values: ({start_x},{end_x})')
        N, Nx = lin_shape_iso(gauss_xi)
        J = jacobian(np.array([start_x,end_x]),Nx)
        for g in range(nodes_per_elem):
            sv[g] = df_bfs.iloc[i]['force']*N[0,g] + df_bfs.iloc[i+1]['force']*N[1,g]
            area[g] = df_ndprops.iloc[i]['area']*N[0,g] + df_ndprops.iloc[i+1]['area']*N[1,g]
        f = np.zeros(nodes_per_elem)
        for m in range(2):
            for g in range(nodes_per_elem):
                f[m] += sv[g]*N[m,g]*gauss_w[g]*J[g]
        k = np.zeros((nodes_per_elem,nodes_per_elem))
        for row in range(nodes_per_elem):
            for col in range(nodes_per_elem):
                for g in range(nodes_per_elem):
                    k[row,col] += area[g]*E*Nx[row,g]*Nx[col,g]*gauss_w[g]/J[g]
        K, F = assemble_global_stiffness(K,F,k,f,i,df_elements)

    for i in range(n_nodes):
        if df_bcs.iloc[i]['bc type'] == 1:
            # print(f'Type 1 BC on node {i}')
            penalty = abs(K[i,i]+1)*1e7
            K[i,i] = penalty
            print(K)
            F[i] += df_bcs.iloc[i]['bc val']*penalty
        if df_bcs.iloc[i]['bc type'] == 2:
            print(f'Type 2 BC on node {i}')
            F[i] += df_bcs.iloc[i]['bc val']*df_ndprops.iloc[i]['area']
    u = np.linalg.solve(K,F)
    print(u)
    plt.figure(1)
    plt.plot(u,marker='o',mfc='w',mec='k')
    plt.xlabel('Node number')
    plt.ylabel('Displacement (in)')
    plt.grid()
    plt.show()

    
