import numpy as np
import matplotlib.pyplot as plt

class TriMesher:
    def __init__(self):
        self.edges = np.zeros((1,23))
        self.nodes = np.zeros((1,23))
        self.ray_angles = np.zeros((1,20))
        self.tris = np.zeros((1,6))
        self.edges = np.zeros((1,10))
        self.last_node_row = 0
        self.last_tri_row = 0

    def QueueCreate(self,n_col):
        # initial queue size == 10
        queue = np.zeros((10,n_col))
        head = 1
        tail = 0
        return queue, head, tail

    def QueueRemove(self,queue,head,tail):
        if head > tail:
            element = np.nan
        else:
            element = queue[head,:]
            head += 1
        return queue, head, tail, element

    def QueueAdd(self, queue, head, tail, element):
        q_size = queue.shape[0]
        length = tail - head
        if length == -1: # if queue is empty
            queue[0,:] = element
            head = 1
            tail = 0
        if tail == q_size:  # if tail is close to the end
            if tail - head < 0.5*q_size: # if the queue is less than half-filled
                queue[:tail-head+1,:] = queue[head:tail,:] # move queue forward
                head = 1
                tail = head + length
            else:
                t = np.zeros((2*q_size,queue.shape[1])) # if more than half filled, double the size
                t[:tail,:] = queue
                queue = t
        tail += 1
        queue[tail,:] = element # add the element
        return queue, head, tail

    def AddNode(self,x,y):
        row = self.last_node_row + 1 # indexes the next row
        if row > self.nodes.shape[0]: # if the nodes matrix is full, double the length
            tmp = nodes
            nodes = np.zeros((self.last_node_row*2,23))
            nodes[:self.last_node_row,:] = tmp
            tmp = self.ray_angles
            self.ray_angles = np.zeros((self.last_node_row*2,20))
            self.ray_angles[:self.last_node_row,:] = tmp

        nodes[row,:2] = [x,y]
        self.last_node_row += 1
        return None

    def AddTriangle(self,nids,eids):
        row = self.last_tri_row + 1
        if row > self.tris.shape[0]: # if the tris matrix is full, double the length
            tmp = self.tris
            self.tris = np.zeros((self.last_tri_row*2,6))
            self.tris[:self.last_tri_row,:] = tmp
        
        # Add triangle information to the last non-empty row
        self.tris[row,:3] = nids
        self.tris[row,3:6] = eids
        for i in range(3,6):
            eid = self.tris[row,i]
            if self.edges[eid,8] == 0:
                self.edges[eid,8] = row
            elif self.edges[eid,9] == 0:
                self.edges[eid,9] = row
            else:
                raise ValueError('Add tris error')

        self.last_tri_row += 1
        return None

    def IsExistingEdge(self,a,b):
        for i in range(self.nodes[a,2]):
            if self.edges[self.nodes[a,3+i],0] == b | self.edges[self.nodes[a,3+i],1] == b:
                re = self.nodes[a,3+i]
                return re
            else:
                return 0

    def NeighborEdges(self,a):
        oprev = self.edges[a,3]
        onext = self.edges[a,4]
        dprev = self.edges[a,5]
        dnext = self.edges[a,6]
        return oprev, onext, dprev, dnext

    def NeighborEdgesNodes(self,eid):
        oprev, onext, dprev, dnext = self.NeighborEdges(eid)
        o = self.edges[eid,0]
        d = self.edges[eid,1]
        nop = self.GetRayEndNode(o,oprev)
        non = self.GetRayEndNode(o,onext)
        ndp = self.GetRayEndNode(d,dprev)
        ndn = self.GetRayEndNode(d,dnext)

        if self.edges[eid,2] == 2: # flagged as a boundary edge
            nids = np.zeros(3)
            if nop == ndn: # neighbor edges
                eids = [oprev, dnext]
            else:
                eids = [dprev, onext]
        elif self.edges[eid,2] == 1: # flagged as an interior edge
            nids = np.zeros(4)
            eids = [oprev, dnext, dprev, onext]

        # get associated nodes
        nids[0] = o
        k = 1
        if nop == ndn:
            nids[k] = nop
            k += 1
        nids[k] = d
        k += 1
        if non == ndp:
            nids[k] = non
        return eids, nids

    def GetRayEndNode(self,nd,edge_id):
        # get the end node on the edge that is away from node nd
        if self.edges[edge_id,0] == nd:
            re = self.edges[edge_id,1]
        else:
            re = self.edges[edge_id,0]
        return re

    def GetNbrRayEndNode(self,center_nd,ray_nd_in,flag):
        n_rays = self.nodes[center_nd,2]
        for i in range(n_rays):
            edge_id = self.nodes[center_nd,i+3]
            if self.edges[edge_id,0] == ray_nd_in | self.edges[edge_id,1] == ray_nd_in:
                if flag == 1:
                    nbr_edge_id = self.nodes[center_nd,4+i%n_rays]
                else:
                    nbr_edge_id = self.nodes[center_nd,4+(i-2+n_rays)%n_rays]
                nid = self.GetRayEndNode[center_nd,nbr_edge_id]
                eid = nbr_edge_id
                return eid, nid
        nid = center_nd
        print('Failed ray operation\n')
        return None

    def InsertRay(self,nd,type,edge_row,edge_angle):
        n_rays = self.nodes[nd,2]
        # insert the ray in the node record
        if n_rays == 0:
            self.nodes[nd,2] = 1
            self.nodes[nd,3] = edge_row
            self.ray_angles[nd,0] = edge_angle
        else:
            tmp = np.array([self.nodes[nd,3:3+n_rays], edge_row]).reshape(-1,1)
            tmp2 = np.array([self.ray_angles[nd,:n_rays], edge_angle]).reshape(-1,1)
            t = np.concatenate((tmp,tmp2),axis=1)
            idx = np.argsort(t[:,1],axis=0)
            t = t[idx]
            self.nodes[nd,3:3+n_rays+1] = t[:,0]
            self.ray_angles[nd,:n_rays+1] = t[:,1]
            self.nodes[nd,2] = self.nodes[nd,2] + 1
        # set up the oprev, onext, dprev, dnext entries
        for i in range(n_rays):
            if edge_row == self.nodes[nd,3+i]:
                if self.edges[edge_row,2] == 1:
                    self.edges[edge_row,3] = self.nodes[nd,4+(i-2+n_rays)%n_rays]
                    self.edges[edge_row,4] = self.nodes[nd,4+i%n_rays]
                    orig = self.edges[edge_row,0]
                    opre_edge = self.edges[edge_row,3]
                    if self.edges[opre_edge,0] == orig:
                        self.edges[opre_edge,4] = edge_row
                    else:
                        self.edges[opre_edge,6] = edge_row
                    onext_edge = self.edges[edge_row,4]
                    if self.edges[onext_edge,0] == orig:
                        self.edges[onext_edge,3] = edge_row
                    else:
                        self.edges[onext_edge,5] = edge_row
                else:
                    self.edges[edge_row,5] = self.nodes[nd,4+(i-2+n_rays)%n_rays]
                    self.edges[edge_row,6] = self.nodes[nd,4+i%n_rays]
                    dest = self.edges[edge_row,1]
                    dprev_edge = self.edges[edge_row,5]
                    if self.edges[dprev_edge,0] == dest:
                        self.edges[dprev_edge,4] = edge_row
                    else:
                        self.edges[dprev_edge,6] = edge_row
                    dnext_edge = self.edges[edge_row,6]
                    if self.edges[dnext_edge,0] == dest:
                        self.edges[dnext_edge,3] = edge_row
                    else:
                        self.edges[dnext_edge,5] = edge_row 
