import numpy as np
import pandas as pd
import os
import csv

def GenerateElements(file):
    df = pd.read_csv(file)
    nelem = len(df)-1
    elements = np.zeros((nelem,3),dtype='int8')
    elem1 = np.arange(0,nelem)
    elements[:,0] = elem1
    elements[:,1] = elem1
    elem2 = np.arange(1,nelem+1)
    elements[:,2] = elem2
    return elements

def write_elements(array,name='elements.csv'):
    path = '/Users/mattjwilliams/Documents/PythonStuff/FEM/GangLi/Chapter4/ProgramFiles'
    fname = os.path.join(path,name)
    headers = ['element','node1','node2']
    with open(fname,'w',newline='') as f:
        mywriter = csv.writer(f,delimiter=' ')
        mywriter.writerow(headers)
        mywriter.writerows(array)
    return None

if __name__ == '__main__':
    fpath = '/Users/mattjwilliams/Documents/PythonStuff/FEM/GangLi/Chapter4/ProgramFiles/nodes.csv'
    if not os.path.exists(fpath):
        raise IOError('Nodes file does not yet exist. Run 1Dmesher first.')
    elems = GenerateElements(fpath)
    write_elements(elems)
    print('Element file creation complete.')