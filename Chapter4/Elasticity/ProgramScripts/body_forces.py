import numpy as np
import csv
import os
import pandas as pd

def ComputeBodyForces(nodes_df):
    force_df = pd.DataFrame()
    forces = np.zeros(len(nodes_df))
    for i in range(len(nodes_df)):
        forces[i] = 500*(nodes_df.iloc[i]['x-location'])**(1/3)
    force_df['node'] = np.arange(len(nodes_df))
    force_df['force'] = forces
    return force_df

if __name__ == '__main__':
    root_path = os.path.join(os.getcwd(),'GangLi','IntroToFEM','Chapter4')
    path = os.path.join(root_path,'ProgramFiles')
    # path = '/Users/mattjwilliams/Documents/PythonStuff/FEM/GangLi/Chapter4/ProgramFiles'
    fname = 'nodes.csv'
    if not os.path.exists(os.path.join(path,fname)):
        raise Exception('Nodes file does not exist.')
    else:
        node_df = pd.read_csv(os.path.join(path,fname),delimiter=' ')
        df = ComputeBodyForces(node_df)
        df.to_csv(os.path.join(path,'bfs.csv'),sep=' ',index=False)
    