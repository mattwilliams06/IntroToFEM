import numpy as np
import csv
import os

from sympy import root

def CreateMesh1D(xstart,xend,n_elements):
    dx = (xend-xstart)/n_elements
    nodes = np.zeros((n_elements+1,2),dtype='object')
    node_locations = np.arange(xstart,xend+dx,dx)
    global_nodes = np.arange(len(node_locations),dtype='int8')
    nodes[:,0] = global_nodes
    nodes[:,1] = node_locations
    return nodes

def write_mesh(array,name='nodes.csv'):
    root_path = os.path.join(os.getcwd(),'GangLi','IntroToFEM','Chapter4')
    path = os.path.join(root_path,'ProgramFiles')
    # path = '/Users/mattjwilliams/Documents/PythonStuff/FEM/GangLi/Chapter4/ProgramFiles'
    fname = os.path.join(path,name)
    headers = ['node','x-location']
    with open(fname,'w',newline='') as f:
        mywriter = csv.writer(f,delimiter=' ')
        mywriter.writerow(headers)
        mywriter.writerows(array)
    return None

if __name__ == '__main__':
    print('Beginning mesher script.')
    print('A CSV file will be generated mapping node locations\n'
    'to global node numbers.')
    print('Enter the beginning x location, ending x location, and number of elements:')
    xstart = float(input('Beginning: '))
    xend = float(input('End: '))
    nelem = int(input('Number of elements: '))
    nodes = CreateMesh1D(xstart,xend,nelem)
    write_mesh(nodes)
    print('File creation complete.')
    