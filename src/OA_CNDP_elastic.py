import time
from src import Params
from docplex.mp.model import Model

class OA_CNDP_elastic:
    
    def __init__(self, network):
        self.network = network
        
        
        for a in self.network.links2:
            a.y = 0
            a.add_cap = 0
            
        self.varlinks = self.network.links
        
        self.x_target = dict()
        self.q_target = dict()
        
        for line in open("data/"+self.network.name+"/linkflows_1.txt", "r"):
            data = line.split()
            i = int(data[0])
            j = int(data[1])
            x = float(data[2])
            self.x_target[self.network.findLink(i, j)] = x
            
        for line in open("data/"+self.network.name+"/demand_1.txt", "r"):
            data = line.split()
            r = int(data[0])
            s = int(data[1])
            q = float(data[2])
            self.q_target[(self.network.findNode(r), self.network.findNode(s))] = q    
        
 
        
        self.network.params.equilibrate_demand = True
    
    def solve(self):
       
        timelimit = 3600
        iteration = 0
        starttime = time.time()
        ub = 1e15
        lb = 0
        gap = 1
        cutoff = 0.01
        
        last_x_f = {a:0 for a in self.network.links}
        last_y_l = {a:0 for a in self.varlinks}
        last_x_l = {a:0 for a in self.network.links}
        
        self.initRMP()
        
        while gap > cutoff:
            iteration += 1
            
            # solve RMP -> y, LB
            x_l, y_l, q_l, obj_l = self.solveRMP()
            lb = obj_l
            ll_l = self.calcBeckmann(x_l, y_l)
            
            # solve TAP -> x, UB
            x_f, ll_f = self.TAP(y_l)

            # add VF cut
            self.addVFCut(x_l, q_l)
            # add TSTT cut
            self.addObjCut(x_l, q_l)
            

            
            elapsed = time.time() - starttime
            if lb > 0:
                gap = (ub - lb)/lb
            else:
                gap = 1
            
            print(iteration, lb, ub, gap, elapsed, ll_f, ll_l - ll_f)
            
            #for a in self.varlinks:
            #    print("\t", a, yhat[a], a.C/2)
            
 
            if elapsed > timelimit:
                break
                
            last_x_l = x_l
            last_y_l = y_l
            last_x_l = x_l
    
    def calcOFV(self):
        output = self.network.getTSTT("UE")
        
        for a in self.varlinks:
            output += self.g[a] * a.add_cap
        
        return output

    def calcBeckmann(self, xhat, yhat):
        total = 0
        
        for a in self.network.links:
            if a in self.varlinks:
                total += a.getPrimitiveTravelTimeC(xhat[a], yhat[a])
            else:
                total += a.getPrimitiveTravelTimeC(xhat[a], 0)
                
        return total
            
    def addVFCut(self, x_l, xhat, yhat):
        yhat_ext = dict()
        
        B1 = self.calcBeckmann(x_l, yhat)
        B2 = self.calcBeckmann(xhat, yhat)
        
        for a in self.network.links:
            if a in self.varlinks:
                yhat_ext[a] = yhat[a]
            else:
                yhat_ext[a] = 0
                
        self.rmp.add_constraint(B1-B2 + sum( (self.rmp.x[a] - x_l[a]) * a.getTravelTimeC(x_l[a], yhat_ext[a], "UE") for a in self.network.links) + sum( (self.rmp.y[a] - yhat[a]) * (a.intdtdy(x_l[a], yhat[a]) - a.intdtdy(xhat[a], yhat[a])) for a in self.varlinks) <= 0)
    

      
    def addObjCut(self, xhat, qhat):
        for a in self.network.links:
            self.rmp.add_constraint(mu_a[a] >= (xhat[a]-x_target[a]) * (xhat[a]-x_target[a]) + (self.rmp.x[a] - xhat[a])* 2*(xhat[a]-x_target[a]))
        
    
        for r in self.network.origins:
            for s in r.getDests():
                self.rmp.add_constraint(mu_w[(r,s)] >= (qhat[(r,s)]-q_target[(r,s)]) * (qhat[(r,s)]-q_target[(r,s)]) + (self.rmp.q[(r,s)] - qhat[(r,s)])* 2*(qhat[(r,s)]-q_target[(r,s)]))
      
    def initRMP(self):   
        self.rmp = Model()
        self.rmp.mu_a = {a:self.rmp.continuous_var(lb=0,ub=1e10) for a in self.network.links}
        self.rmp.mu_w = {(r,s):self.rmp.continuous_var(lb=0,ub=1e10) for r in self.network.origins for s in r.getDests()}
        #self.rmp.eta = {a:self.rmp.continuous_var(lb=-1e10,ub=1e10) for a in self.network.links}
        self.rmp.y = {a:self.rmp.continuous_var(lb=0, ub=a.max_add_cap) for a in self.varlinks}
        self.rmp.x = {a:self.rmp.continuous_var(lb=0, ub=self.network.TD) for a in self.network.links}
        self.rmp.xc = {(a,r):self.rmp.continuous_var(lb=0, ub=r.totaldemand) for a in self.network.links for r in self.network.origins}
        self.rmp.q = {(r,s): self.rmp.continuous_var(lb=0, ub=2*r.getDemand(s)) for r in self.network.origins for s in r.getDests()}
        self.rmp.y = {(r,s): self.rmp.continuous_var(lb=0, ub=2*r.getDemand(s)) for r in self.network.origins for s in r.getDests()}
        self.rmp.beta = {a:self.rmp.continuous_var(lb=0,ub=1e10) for a in self.network.links}
        self.rmp.rho = {(r,s):self.rmp.continuous_var() for r in self.network.origins for s in r.getDests()}
        
        
        
        for a in self.network.links:
            self.rmp.add_constraint(sum(self.rmp.xc[(a,r)] for r in self.network.origins) == self.rmp.x[a])
            
        for i in self.network.nodes:                    
            for r in self.network.origins:            

                if i.id == r.id:
                    dem = - sum(self.rmp.q[(r,s)] for s in self.network.zones)                
                elif isinstance(i, type(r)) == True:
                    dem = self.rmp.q[(r,i)]
                else:
                    dem = 0

                self.rmp.add_constraint(sum(self.rmp.xc[(a,r)] for a in i.incoming) - sum(self.rmp.xc[(a,r)] for a in i.outgoing) == dem)
        
        
        self.rmp.minimize(sum(self.rmp.mu_a[a] for a in self.network.links) + sum(self.rmp.mu_w[(r,s)] for r in self.network.origins for s in r.getDests()))
        
    def solveRMP(self):
    
        t_solve = time.time()
        self.rmp.solve(log_output=False)
        t_solve = time.time() - t_solve
        
        yhat = {a:self.rmp.y[a].solution_value for a in self.varlinks}
        x_l = {a:self.rmp.x[a].solution_value for a in self.network.links}
        obj_l = self.rmp.objective_value
        
        return x_l, yhat, obj_l
        
    def TAP(self, y):
    
        for a in self.varlinks:
            a.add_cap = y[a]
            
        self.network.tapas("UE", None)
        xhat = {a:a.x for a in self.network.links}
        obj_f = self.calcOFV()
        
        return xhat, obj_f
 