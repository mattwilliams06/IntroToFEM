import matplotlib.pyplot as plt
import numpy as np
import os

def create_array(path):
    '''
    Text file at path must be in PSLG format
    '''
    os.chdir(path)
    data = []
    with open('shape.txt','r') as f:
        d = f.readlines()
        for i in d:
            k = i.rstrip().split(' ')
            data.append([int(j) for j in k])
    data = np.array(data)
    return data

def draw_shape(array):
    n_shapes = array[0,0]
    n_outer_pts = array[1,0]
    n_inner_pts = array[2,0]
    outer_x = array[3:3+n_outer_pts,0]
    outer_y = array[3:3+n_outer_pts,1]
    inner_x = array[3+n_outer_pts:,0]
    inner_y = array[3+n_outer_pts:,1]
    for i in range(n_outer_pts):
        x1,y1 = outer_x[i], outer_y[i]
        x2,y2 = outer_x[(i+1)%n_outer_pts],outer_y[(i+1)%n_outer_pts]
        plt.plot([x1,x2],[y1,y2],'k')
    for i in range(n_inner_pts):
        x1,y1 = inner_x[i], inner_y[i]
        x2,y2 = inner_x[(i+1)%n_inner_pts],inner_y[(i+1)%n_inner_pts]
        plt.plot([x1,x2],[y1,y2],'k')
    plt.show()

if __name__ == '__main__':
    path = '/Users/mattjwilliams/Documents/PythonStuff/FEM/GangLi/Chapter6/'
    arr = create_array(path)
    draw_shape(arr)

