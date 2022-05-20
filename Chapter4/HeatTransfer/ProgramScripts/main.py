import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import csv

### NOTE: Need to make filepaths for the laptop now ###

def lagrange(X):
    '''
    Function for creating Lagrange basis functions.
    Inputs
    ------
    X: array-like. The elements nodal x-locations. 
    Returns
    -------
    List of basis functions and a list of their derivatives. This function supports linear
    and quadratic basis functions. 
    '''
    deg = len(X)-1
    if deg == 1:
        x0 = X[0]
        x1 = X[1]
        x = np.linspace(x0,x1)
        phi0 = (x-x1)/(x0-x1)
        phi1 = (x-x0)/(x1-x0)
        dphi0 = 1/(x0-x1)
        dphi1 = 1/(x1-x0)
        return [phi0, phi1], [dphi0*np.ones_like(x), dphi1*np.ones_like(x)], x
    elif deg ==2:
        x0, x1, x2 = X
        x = np.linspace(x0,x2)
        phi0 = (x-x1)*(x-x2)/((x0-x1)*(x0-x2))
        phi1 = (x-x0)*(x-x2)/((x1-x0)*(x1-x2))
        phi2 = (x-x0)*(x-x1)/((x2-x0)*(x2-x1))
        dphi0 = (2*x-x1-x2)/((x0-x1)*(x0-x2))
        dphi1 = (2*x-x0-x2)/((x1-x0)*(x1-x2))
        dphi2 = (2*x-x0-x1)/((x2-x0)*(x2-x1))
        return [phi0, phi1, phi2], [dphi0, dphi1, dphi2], x
    else:
        raise ValueError('Currently supporting only first or second order Lagrange basis functions.')

# def read_nodes(file):
#     path = os.getcwd()
#     if 'mattjwilliams' in path:
#         path = '/Users/mattjwilliams/Documents/PythonStuff/FEM/GangLi/Chapter4/HeatTransfer/ProgramFiles'
#     node_df = pd.read_csv(os.path.join(path,file),delimiter=' ')
#     return node_df

def read_file(file):
    path = os.getcwd()
    if 'mattjwilliams' in path:
        path = '/Users/mattjwilliams/Documents/PythonStuff/FEM/GangLi/Chapter4/HeatTransfer/ProgramFiles'
    else:
        path = '/Users/mattwilliams/Documents/PythonProjects/FEM/GangLi/IntroToFEM/Chapter4/HeatTransfer/ProgramFiles'
    df = pd.read_csv(os.path.join(path,file),delimiter=' ')
    return df

def mesher1D(xstart,xend,n_elements,deg=1):
    if deg == 1:
        n_nodes = n_elements + 1
        dx = (xend-xstart)/(n_elements)
        elements = np.arange(0,n_elements,dtype='int8')
        arr = np.zeros((n_elements,deg+2),dtype='object')
        node1 = np.arange(xstart,xend,dx)
        node2 = np.arange(xstart+dx,xend+dx,dx)
        arr[:,0] = elements
        arr[:,1] = node1
        arr[:,2] = node2
    else:
        n_nodes = 2*n_elements + 1
        dx = (xend-xstart)/(n_elements)
        elements = np.arange(0,n_elements,dtype='int8')
        arr = np.zeros((n_elements,deg+2),dtype='object')
        node1 = np.arange(xstart,xend,dx)
        node2 = np.arange(xstart+dx/2,xend+dx/2,dx)
        node3 = np.arange(xstart+dx,xend+dx,dx)
        arr[:,0] = elements
        arr[:,1] = node1
        arr[:,2] = node2
        arr[:,3] = node3
    path = os.getcwd()
    if 'mattjwilliams' in path:
        path = '/Users/mattjwilliams/Documents/PythonStuff/FEM/GangLi/Chapter4/HeatTransfer/ProgramFiles'
    else:
        path = '/Users/mattwilliams/Documents/PythonProjects/FEM/GangLi/IntroToFEM/Chapter4/HeatTransfer/ProgramFiles'
    fname = 'nodes.txt'
    if deg == 1:
        headers = ['element','node0','node1']
    else:
        headers = ['element','node0','node1','node2']
    with open(os.path.join(path,fname),'w',newline='') as f:
        mywriter = csv.writer(f,delimiter=' ')
        mywriter.writerow(headers)
        mywriter.writerows(arr)
    return None

def element_indexer(n_elements,deg):
    arr = np.zeros((n_elements,deg+2),dtype='int8')
    arr[:,0] = np.arange(n_elements)
    if deg == 1:
        arr[:,1] = np.arange(0,n_elements,deg)
        arr[:,2] = np.arange(1,n_elements+1,deg)
        headers = ['element','node0','node1']
    else:
        n_nodes = 2*n_elements + 1
        arr[:,1] = np.arange(0,n_nodes-1,deg)
        arr[:,2] = np.arange(1,n_nodes,deg)
        arr[:,3] = np.arange(2,n_nodes+1,deg)
        headers = ['element','node0','node1','node2']

    path = os.getcwd()
    if 'mattjwilliams' in path:
        path = '/Users/mattjwilliams/Documents/PythonStuff/FEM/GangLi/Chapter4/HeatTransfer/ProgramFiles'
    else:
        path = '/Users/mattwilliams/Documents/PythonProjects/FEM/GangLi/IntroToFEM/Chapter4/HeatTransfer/ProgramFiles'
    fname = 'elements.txt'
    with open(os.path.join(path,fname),'w',newline='') as f:
        mywriter = csv.writer(f,delimiter=' ')
        mywriter.writerow(headers)
        mywriter.writerows(arr)
    return None

def global_nodes(xstart,xend,n_elements,deg=1):
    if deg == 1:
        n_nodes = n_elements + 1
    else:
        n_nodes = 2*n_elements + 1
    arr = np.zeros((n_nodes,2),dtype='object')
    arr[:,0] = np.arange(n_nodes)
    dx = (xend-xstart)/(n_nodes-1)
    arr[:,1] = np.arange(xstart,xend+dx,dx)
    path = os.getcwd()
    if 'mattjwilliams' in path:
        path = '/Users/mattjwilliams/Documents/PythonStuff/FEM/GangLi/Chapter4/HeatTransfer/ProgramFiles'
    else:
        path = '/Users/mattwilliams/Documents/PythonProjects/FEM/GangLi/IntroToFEM/Chapter4/HeatTransfer/ProgramFiles'
    fname = 'global_nodes.txt'
    headers = ['global node','x location']
    with open(os.path.join(path,fname),'w',newline='') as f:
        mywriter = csv.writer(f,delimiter=' ')
        mywriter.writerow(headers)
        mywriter.writerows(arr)
    return None

def boundary_conditions(n_elements,bcs,deg=1):
    '''
    Boundary conditions passed as a list:
    Left end type, left end value, right end type, right end value

    Type 1 == Dirichlet, Type 2 == Neumann or flux
    '''
    type1, val1, type2, val2 = bcs
    if deg == 1:
        n_nodes = n_elements + 1
    else:
        n_nodes = 2*n_elements + 1
    arr = np.zeros((n_nodes,3),dtype='object')
    arr[:,0] = np.arange(n_nodes)
    arr[0,1] = type1
    arr[0,2] = val1
    arr[-1,1] = type2
    arr[-1,2] = val2
    path = os.getcwd()
    if 'mattjwilliams' in path:
        path = '/Users/mattjwilliams/Documents/PythonStuff/FEM/GangLi/Chapter4/HeatTransfer/ProgramFiles'
    else:
        path = '/Users/mattwilliams/Documents/PythonProjects/FEM/GangLi/IntroToFEM/Chapter4/HeatTransfer/ProgramFiles'
    fname = 'bcs.txt'
    headers = ['global node','type','value']
    with open(os.path.join(path,fname),'w',newline='') as f:
        mywriter = csv.writer(f,delimiter=' ')
        mywriter.writerow(headers)
        mywriter.writerows(arr)
    return None

def properties(n_elements,k1,k2):
    arr = np.zeros((n_elements,2),dtype='object')
    if n_elements%2 != 0:
        arr[:n_elements//2+1,1] = k1
        arr[n_elements//2+1:,1] = k2
        arr[:,0] = np.arange(n_elements)
    else:
        arr[:n_elements//2,1] = k1
        arr[n_elements//2:,1] = k2
        arr[:,0] = np.arange(n_elements)
    path = os.getcwd()
    if 'mattjwilliams' in path:
        path = '/Users/mattjwilliams/Documents/PythonStuff/FEM/GangLi/Chapter4/HeatTransfer/ProgramFiles'
    else:
        path = '/Users/mattwilliams/Documents/PythonProjects/FEM/GangLi/IntroToFEM/Chapter4/HeatTransfer/ProgramFiles'
    fname = 'elprops.txt'
    headers = ['element','kcond']
    with open(os.path.join(path,fname),'w',newline='') as f:
        mywriter = csv.writer(f,delimiter=' ')
        mywriter.writerow(headers)
        mywriter.writerows(arr)
    return None


if __name__ == '__main__':
    ### PREPROCESSING SECTION ###
    # rod is 100 cm long.
    xleft = 0.0
    xright = 100.0
    n_elements = 8
    # degree of basis functions
    deg = 2
    # set my local working directory for retrieving and saving program files
    path = os.getcwd()
    if 'mattjwilliams' in path:
        path = '/Users/mattjwilliams/Documents/PythonStuff/FEM/GangLi/Chapter4/HeatTransfer/ProgramFiles'
    else:
        path = '/Users/mattwilliams/Documents/PythonProjects/FEM/GangLi/IntroToFEM/Chapter4/HeatTransfer/ProgramFiles'
    fnames = ['nodes.txt','elements.txt','global_nodes.txt','bcs.txt','elprops.txt','globals.txt']
    # run functions that create the files containing the mesh information
    mesher1D(xleft,xright,n_elements,deg)
    element_indexer(n_elements,deg)
    global_nodes(xleft,xright,n_elements,deg)
    boundary_conditions(n_elements,[2,0.1,1,0.0],deg)
    properties(n_elements,0.92,0.12)
    # load information from program files into dataframes
    df_dict = {}
    df_names = ['df_nodes','df_elems','df_glob_nodes','df_bcs','df_elprops','df_globals']
    for fname, df_name in zip(fnames,df_names):
        df_dict[df_name] = read_file(fname)
    
    ### MAIN LOOP ###
    n_nodes = n_elements*deg + 1
    K = np.zeros((n_nodes,n_nodes))
    F = np.zeros(n_nodes)
    Tenv, r, kconv = df_dict['df_globals'].iloc[0]
    for e in range(n_elements):
        node_locs = df_dict['df_nodes'].iloc[e,1:]
        N, Nx, xelem = lagrange(node_locs)
        kcond = df_dict['df_elprops'].iloc[e]['kcond']
        n = len(Nx)
        k = np.zeros((n,n))
        f = np.zeros(n)
        for row in range(n):
            I = 2*kconv*Tenv/r*N[row]
            f[row] = np.trapz(I,x=xelem)
            for col in range(n):
                # This currently only has the first integral on the LHS
                I1 = Nx[row]*Nx[col]*kcond
                I2 = 2*kconv/r*N[row]*N[col]
                k[row,col] = np.trapz(I1,x=xelem) + np.trapz(I2,x=xelem)
        glob_nodes = df_dict['df_elems'].iloc[e,1:]
        for i,row in enumerate(glob_nodes):
            F[row] += f[i]
            for j,col in enumerate(glob_nodes):
                K[row,col] += k[i,j]
    for i in range(n_nodes):
        bc_type = df_dict['df_bcs'].iloc[i]['type']
        if bc_type == 2:
            F[i] += df_dict['df_bcs'].iloc[i]['value']
        if bc_type == 1:
            penalty = abs(K[i,i]+1)*1e7
            K[i,i] = penalty
            F[i] = df_dict['df_bcs'].iloc[i]['value']*penalty
    ### END MAIN LOOP ###

    ### SOLUTION AND POST PROCESSING ###
    temps = np.linalg.solve(K,F)
 
    sol = []
    xvals = []
    for e in range(n_elements):
        node_locs = df_dict['df_nodes'].iloc[e,1:]
        N, Nx, xelem = lagrange(node_locs)
        val_arr = np.zeros_like(xelem)
        xvals.append(xelem)
        for i in range(deg+1):
            val_arr += N[i]*temps[deg*e+i]
        sol.append(val_arr)


sol_temps = np.stack(sol).flatten()
sol_xvals = np.stack(xvals).flatten()
plt.plot(sol_xvals,sol_temps)
plt.grid()
plt.xlim(sol_xvals[0],sol_xvals[-1])
plt.ylim(0,30)
plt.xlabel('Distance along bar (cm)')
plt.ylabel('Temperature ($^oC$)')
plt.title(f'FEM Solution with {n_elements} Elements and Degree {deg} Interpolants')
plt.show()
 
