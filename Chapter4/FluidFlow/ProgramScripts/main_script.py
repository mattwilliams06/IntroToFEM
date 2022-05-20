import numpy as np
import matplotlib.pyplot as plt
import os
import csv
import pandas as pd

def lagrange(X):
    deg = len(X) - 1
    if deg == 1:
        x0, x1 = X
        x = np.linspace(x0,x1)
        N0 = (x-x1)/(x0-x1)
        N1 = (x-x0)/(x1-x0)
        dN0 = 1/(x0-x1)
        dN1 = 1/(x1-x0)
        return [N0,N1],[dN0,dN1], x
    else:
        x0,x1,x2 = X
        x = np.linspace(x0,x2)
        N0 = (x-x1)*(x-x2)/((x0-x1)*(x0-x2))
        N1 = (x-x0)*(x-x2)/((x1-x0)*(x1-x2))
        N2 = (x-x0)*(x-x1)/((x2-x0)*(x2-x1))
        dN0 = (2*x-x1-x2)/((x0-x1)*(x0-x2))
        dN1 = (2*x-x0-x2)/((x1-x0)*(x1-x2))
        dN2 = (2*x-x0-x1)/((x2-x0)*(x2-x1))
        return [N0,N1,N2],[dN0,dN1,dN2], x
    

def get_path():
    path = os.getcwd()
    if 'mattwilliams' in path:
        path = '/Users/mattwilliams/Documents/PythonProjects/FEM/GangLi/IntroToFEM/Chapter4/FluidFlow/ProgramFiles'
    return path

def read_file(path,fname):
    df = pd.read_csv(os.path.join(path,fname),delimiter=' ')
    return df

def mesher1D(xstart,xend,nelems,deg):
    if deg == 1:
        arr = np.zeros((nelems,2*(deg+1)+1),dtype='object')
        headers = ['element','node0','node1','x0','x1']
        elems = np.arange(nelems,dtype='int8')
        node0 = np.arange(nelems,dtype='int8')
        node1 = np.arange(1,nelems+1,dtype='int8')
        dx = (xend-xstart)/nelems
        x0 = np.arange(xstart,xend,dx)
        x1 = np.arange(xstart+dx,xend+dx,dx)
        arrs = [elems,node0,node1,x0,x1]
        fname = 'elements.csv'
        path = get_path()
        for i,a in enumerate(arrs):
            arr[:,i] = a
        with open(os.path.join(path,fname),'w') as f:
            mywriter = csv.writer(f,delimiter=' ')
            mywriter.writerow(headers)
            mywriter.writerows(arr)
    else:
        n_nodes = 2*nelems + 1
        arr = np.zeros((nelems,2*(deg+1)+1),dtype='object')
        headers = ['element','node0','node1','node2','x0','x1','x2']
        elems = np.arange(nelems,dtype='int8')
        node0 = np.arange(0,n_nodes-1,2,dtype='int8')
        node1 = np.arange(1,n_nodes,2,dtype='int8')
        node2 = np.arange(2,n_nodes+1,2,dtype='int8')
        dx = (xend-xstart)/nelems
        x0 = np.arange(xstart,xend,dx)
        x1 = np.arange(xstart+dx/2,xend+dx/2,dx)
        x2 = np.arange(xstart+dx,xend+dx,dx)
        arrs = [elems,node0,node1,node2,x0,x1,x2]
        fname = 'elements.csv'
        path = get_path()
        for i,a in enumerate(arrs):
            arr[:,i] = a
        with open(os.path.join(path,fname),'w') as f:
            mywriter = csv.writer(f,delimiter=' ')
            mywriter.writerow(headers)
            mywriter.writerows(arr)

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
        path = '/Users/mattwilliams/Documents/PythonProjects/FEM/GangLi/IntroToFEM/Chapter4/FluidFlow/ProgramFiles'
    fname = 'bcs.csv'
    headers = ['global node','type','value']
    with open(os.path.join(path,fname),'w') as f:
        mywriter = csv.writer(f,delimiter=' ')
        mywriter.writerow(headers)
        mywriter.writerows(arr)
    return None

def build_k(D,v,X):
    N, dN, xelem = lagrange(X)
    nrows = len(X)
    # k = np.zeros((nrows,nrows))
    Aarr = np.zeros((nrows,nrows))
    Darr = np.zeros((nrows,nrows))
    for row in range(nrows):
        for col in range(nrows):
            integrand1 = dN[row]*dN[col]
            if isinstance(integrand1,float):
                integrand1 *= np.ones_like(xelem)*D
            Aarr[row,col] += np.trapz(integrand1,x=xelem)
    print(Aarr)

    for row in range(nrows):
        for col in range(nrows):
            # if isinstance(dN[col],float):
            #     dN[col] *= np.ones_like(xelem)
            integrand2 = N[row]*dN[col]*v
            Darr[row,col] += np.trapz(integrand2,x=xelem)
    k = Aarr + Darr
    return k

def build_K(K,k,eid,df_elem):
    if k.size == 9:
        colnames = ['node0','node1','node2']
    else:
        colnames = ['node0','node1']
    idx = []
    for col in colnames:
        idx.append(int(df_elem.iloc[eid][col]))
    for i,row in enumerate(idx):
        for j,col in enumerate(idx):
            K[row,col] += k[i,j]
    return K

def build_f(s,X):
    N, Nx, xelem = lagrange(X)
    nrows = len(X)
    f = np.zeros(nrows)
    for row in range(nrows):
        integrand = N[row]*s
        f[row] += np.trapz(integrand,x=xelem)
    return f

def build_F(F,f,eid,df_elem):
    cols = df_elem.columns
    if 'node2' in cols:
        colnames = ['node0','node1','node2']
    else:
        colnames = ['node0','node1']
    idx = []
    for col in colnames:
        idx.append(int(df_elem.iloc[eid][col]))
    for i,row in enumerate(idx):
        F[row] += f[i]
    return F


if __name__ == '__main__':
    xl = 0.0
    xr = 1000
    nelems = 5
    deg = 1
    if deg == 1:
        node_names = ['node0','node1']
        x_names = ['x0','x1']
    else:
        node_names = ['node0','node1','node2']
        x_names = ['x0','x1','x2']
    mesher1D(xl,xr,nelems,deg)
    boundary_conditions(nelems,[1,0.01,1,0.0],deg)
    path = get_path()
    df_elem = read_file(path,'elements.csv')
    df_props = read_file(path,'props.csv')
    df_bcs = read_file(path,'bcs.csv')
    K = np.zeros((deg*nelems+1,deg*nelems+1))
    F = np.zeros(deg*nelems+1)
    D = df_props.loc[df_props['property']=='D']['val'].values[0]
    v = df_props.loc[df_props['property']=='v']['val'].values[0]
    s = df_props.loc[df_props['property']=='s']['val'].values[0]

    for e in range(nelems):
        nodes = []
        X = []
        for i in range(len(node_names)):
            nodes.append(int(df_elem.iloc[e][node_names[i]]))
            X.append(df_elem.iloc[e][x_names[i]])
        k = build_k(D,v,X)
        K = build_K(K,k,e,df_elem)
        f = build_f(s,X)
        F = build_F(F,f,e,df_elem)
    print(F)
    n_nodes = deg*nelems + 1
    for i in range(n_nodes):
        bc_type = df_bcs.iloc[i]['type']
        if bc_type == 2:
            F[i] += df_bcs.iloc[i]['value']
        if bc_type == 1:
            penalty = abs(K[i,i]+1)*1e7
            K[i,i] = penalty
            F[i] = df_bcs.iloc[i]['value']*penalty
    C = np.linalg.solve(K,F)
    
    sol = []
    xvals = []
    for e in range(nelems):
        node_locs = df_elem.iloc[e,deg+2:]
        N, Nx, xelem = lagrange(node_locs)
        val_arr = np.zeros_like(xelem)
        xvals.append(xelem)
        for i in range(deg+1):
            val_arr += N[i]*C[deg*e+i]
        sol.append(val_arr)

    sol_temps = np.stack(sol).flatten()
    sol_xvals = np.stack(xvals).flatten()
    plt.plot(sol_xvals,sol_temps)
    plt.grid()
    plt.xlim(sol_xvals[0],sol_xvals[-1])
    # plt.ylim(0,30)
    plt.xlabel('Distance (m)')
    plt.ylabel('Concentration ($mol/m^2$)')
    plt.title(f'FEM Solution with {nelems} Elements and Degree {deg} Interpolants')
    plt.show()

        