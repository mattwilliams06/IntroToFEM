import numpy as np
import csv
import os
import pandas as pd

program_files = '/Users/mattjwilliams/Documents/PythonStuff/FEM/GangLi/Chapter4/ProgramFiles'
fname = 'elements.csv'
df = pd.read_csv(os.path.join(program_files,fname),delimiter=' ')
nelem = len(df['element'])
nnodes = nelem + 1
nodes = np.arange(nnodes)
df_new = pd.DataFrame()
df_new['nodes'] = nodes
df_new['bc type'] = np.zeros_like(nodes)
df_new['bc val'] = np.zeros_like(nodes)

flag = True
print('This script stores the bounday conditions subject to the following:')
print('The first entry is the global node value, indexed to zero.')
print('The second entry is the boundary condition type: 0 for none, 1 for Dirichlet, 2 for Neumann.')
print('The third entry is the value.')
while flag:
    node = int(input('Please enter the node: '))
    bc_type = int(input('Please enter the type: '))
    val = float(input('Please enter the value: '))
    df_new.at[node,'bc type'] = bc_type
    df_new.at[node,'bc val'] = val
    answer = input('Any more boundary conditions to enter? (yes/no)')
    answer = answer.lower()
    if 'n' in answer:
        flag = False
df_new.to_csv(os.path.join(program_files,'bcs.csv'),sep=' ',index=False)

