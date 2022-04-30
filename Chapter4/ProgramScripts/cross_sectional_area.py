import numpy as np
import os
import pandas as pd
import csv

def compute_area(path,file):
    df = pd.read_csv(os.path.join(path,file))
    area_start = 2.0
    area_end = 4.0
    length = 48.0
    area_slope = (area_end-area_start)/length
    area_func = lambda x: 2.0 + area_slope*x
    nelem = 5
    dx = length/nelem
    nodes = np.arange(nelem+1)
    ndprops = np.zeros((len(nodes),2),dtype='object')
    ndprops[:,0] = nodes
    ndprops[:,1] = area_func(ndprops[:,0]*dx)
    return ndprops

def write_elements(array,path,name='ndprops.csv'):
    fname = os.path.join(path,name)
    headers = ['node','area']
    with open(fname,'w',newline='') as f:
        mywriter = csv.writer(f,delimiter=' ')
        mywriter.writerow(headers)
        mywriter.writerows(array)
    return None

if __name__ == '__main__':
    root_path = os.path.join(os.getcwd(),'GangLi','IntroToFEM','Chapter4')
    path = os.path.join(root_path,'ProgramFiles')
    # path = '/Users/mattjwilliams/Documents/PythonStuff/FEM/GangLi/Chapter4/ProgramFiles'
    file = 'nodes.csv'
    arr = compute_area(path,file)
    write_elements(arr,path)
    print('File creation complete.')