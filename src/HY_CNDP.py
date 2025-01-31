import time
from src import Params
from docplex.mp.model import Model
import polytope as pc
import numpy as np
import math


class HY_CNDP:
    
    def __init__(self, network):
        self.network = network
        
        
        self.g = {a:a.cost for a in self.network.links}
        
        for a in self.network.links2:
            a.y = 0
            a.add_cap = 0
            
        self.varlinks = []
        
        
        for a in self.network.links:
            if a.cost > 0:
                self.varlinks.append(a)
                
        
    
    def solve(self):
        time1 = time.time()
        self.solveSLP(1)
        time1 = time.time() - time1
        
        print("time", time1)
    
    
    def solveSLP(self, muk):
        
        t1_min = 0
        t1_max = 10
        
        t2 = time.time()
        polytope_r = {r: self.constructS0ExtremePoints(r) for r in self.network.origins}
        extreme_r = {r: pc.extreme(polytope_r[r]) for r in self.network.origins}
        
        t2 = time.time()-t2
        
        print("extreme points", t2)
        
        f = open("HY_time.txt", "w")
        f.write(str(t2))
        f.close()
        
        '''
        
        # step 1 choose tolerance factor delta k > 0
        delta_k = 1
        
        iteration = 0
        
        # let polytope S_o as defined in (6)
        # calculate penalty term H with all-or-nothing flow assignment
        while True:

            iteration += 1

            SLP_obj, tdotv, t1, best_x, best_y = self.evalExtremePoints(extreme_r, t1_min, t1_max)

                

            
            for a in self.network.links:
                print("\t", a, best_x[a])
            # if tdotv - t1 <= delta^k, terminate with z^k
            if tdotv - t1 <= delta_k:
                return SLP_obj
            
            break
                
            # compute multicutting plane
            # add linear constraints to 
    
        '''
            
    def evalExtremePoints(self, extreme_r, t1_min, t1_max):
        # construct combinations
        
        best_x = None
        best_y = None
        best_obj = 1e100
        best_tdotv = 0
        best_t1 = 0
            
        options = dict()
        combinations = 1
        
        for r in self.network.origins:
            options[r] = combinations
            combinations *= len(extreme_r[r])
            
            
        
        # y min, max
        combinations *= pow(2, len(self.varlinks))
        
        # t1 min, max
        combinations *= 2
        
        print("combinations", combinations, sum(len(extreme_r[r]) for r in self.network.origins))
        
        for i in range(0, combinations):
            # construct extreme point
            x = {a: 0 for a in self.network.links}
            y = {a: 0 for a in self.varlinks}
            
            
            for r in self.network.origins:
                # divide out y choices and t1 choices
                red = math.floor(i / 2 / pow(2, len(self.varlinks)))

                red = math.floor(red / options[r])
                
                
                col = red % len(extreme_r[r])
                
                xc = extreme_r[r][col]
                
                for ij in self.network.links:
                    x[ij] += xc[ij.id]
            
            
            for j in range(0, len(self.varlinks)):
                # at minimum, divide by 2
                # e.g. 0,1,2,3 -> 0 4,5,6,7 -> 1
                # e.g. 0,1 -> 1  2,3-> 1
                red = math.floor(i / pow(2, len(self.varlinks) - j))
                
                if red % 2 == 0:
                    y[self.varlinks[j]] = 0
                else:
                    y[self.varlinks[j]] = self.varlinks[j].max_add_cap
            
            if i % 2 == 0:
                t1 = t1_min
            else:
                t1 = t1_max
                
            print(col, y, t1)
            
            for a in self.network.links:
                print("\t", a, x[a])
            
            SLP_obj, tdotv, t1 = self.evaluate(x, y, t1)  
                
            if SLP_obj < best_obj:
                best_x = x
                best_y = y
                best_obj = SLP_obj
                best_tdotv = tdotv
                best_t1 = t1
                
            
            
        return best_obj, best_tdotv, best_t1, best_x, best_y
            
            
    def evaluate(self, x, y, t1):
        for a in self.network.links:
            a.x = x[a]
            
        for a in self.varlinks:
            a.add_cap = y[a]
            
 
        
        tdotv = sum(a.getTravelTimeC(a.x, a.add_cap, "UE") * a.x for a in self.network.links)
        
        output = tdotv
        return output, tdotv, t1
        
    
    def constructS0ExtremePoints(self, r):
        # this is link flow extreme points for a specific origin
        # splitting by origin to reduce the problem size
        
        xvarslen = len(self.network.links) # x per origin
  
        
        num_cons = len(self.network.nodes) 
        cols = xvarslen
        
        
        rows = num_cons + xvarslen * 2
        
        # v first
        
        A = np.zeros((rows, cols))
        b = np.zeros(rows)
        
        rowidx = 0
        
        for j in self.network.nodes:



            row1 = A[rowidx]
            #row2 = A[rowidx+1]



            dem = 0
            if j.id == r.id:
                dem = - sum(r.getDemand(s) for s in self.network.zones)                
            elif isinstance(j, type(r)) == True:
                dem = r.getDemand(j)




            #print(rowidx, "node", j, "origin", r.id)

            if dem <= 0:
                for ij in j.incoming:
                    #print("\tin", ij.id, ij)
                    varidx = ij.id
                    row1[varidx] = 1
                    #row2[varidx] = -1
                for jk in j.outgoing:
                    #print("\tout", jk.id, jk)
                    varidx = jk.id 
                    row1[varidx] = -1
                    #row2[varidx] = 1

                b[rowidx] = dem
                rowidx += 1
            elif dem > 0:    
                for ij in j.incoming:
                    print("\tin", ij.id, ij)
                    varidx = ij.id
                    row1[varidx] = -1
                    #row2[varidx] = -1
                for jk in j.outgoing:
                    print("\tout", jk.id, jk)
                    varidx = jk.id 
                    row1[varidx] = 1
                    #row2[varidx] = 1
                b[rowidx] = -dem
                rowidx += 1

            #b[rowidx+1] = -dem
                
                
        for j in self.network.links:
                
            col = j.id
            A[rowidx][col] = 1
            b[rowidx] = r.totaldemand
            
            rowidx+= 1
            A[rowidx][col] = -1
            
            rowidx += 1
            
        
        p = pc.Polytope(A, b)
        
        return p