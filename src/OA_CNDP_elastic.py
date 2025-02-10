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
        
        

        self.vfcut_x = list()
        self.vfcut_q = list()
    
    def solve(self):
       
        timelimit = 3600
        max_iter = 30
        iteration = 0
        starttime = time.time()
        ub = 1e15
        lb = 0
        gap = 1
        cutoff = 0.01
        
        last_x_f = {a:0 for a in self.network.links}
        last_y_l = {(r,s):0 for r in self.network.origins for s in r.getDests()}
        last_x_l = {a:0 for a in self.network.links}
        
        self.initRMP()
        
        while gap > cutoff and iteration < max_iter:
            iteration += 1
            
            # solve RMP -> y, LB
            print("\tsolving RMP")
            x_l, q_l, y_l, obj_l = self.solveRMP()

            lb = obj_l

            ll_l = self.calcLLobj(x_l, q_l, y_l)
            
            #print(y_l)
            
            # solve TAP -> x, UB
            print("\tTAP")
            x_f, q_f, obj_f = self.TAP(y_l)
            ll_f = self.calcLLobj(x_f, q_f, y_l)

            # add VF cut
            self.addVFCut(x_l, x_f, q_l, q_f, y_l)
            # add TSTT cut
            #self.addObjCut(x_l, q_l)
            
            ub = min(ub, obj_f)
            
            elapsed = time.time() - starttime
            if lb > 0:
                gap = (ub - lb)/lb
            else:
                gap = 1
            
            #print(obj_f)
            print(iteration, lb, ub, gap, elapsed, ll_l, ll_f)
            
            #for a in self.varlinks:
            #    print("\t", a, yhat[a], a.C/2)
            
 
            if elapsed > timelimit:
                break
                
                
            delta_y = 0
            total_y = 0
            for r in self.network.origins:
                for s in r.getDests():
                    delta_y += abs(y_l[(r,s)] - last_y_l[(r,s)])
                    total_y += y_l[(r,s)]
            
            delta_y_pct = delta_y / total_y
            print("\tdelta_y", delta_y_pct)   
            last_x_l = x_l
            last_y_l = y_l
            last_x_l = x_l
            
            
    
    def calcOFV(self):
        output = sum((self.rmp.x[a].solution_value - self.x_target[a]) ** 2 for a in self.network.links) 
        output +=  sum( (self.rmp.q[(r,s)].solution_value - self.q_target[(r,s)]) ** 2 for r in self.network.origins for s in r.getDests())
        

        return output

    def calcLLobj(self, xhat, qhat, yhat):
        total = 0
        
    
                    
        total += sum(a.getPrimitiveTravelTimeC(xhat[a], 0) for a in self.network.links)
        total += sum(r.intDinv(s, qhat[(r,s)], yhat[(r,s)]) for r in self.network.origins for s in r.getDests())  
        
        
        return total
            
    def addVFCut(self, x_l, x_f, q_l, q_f, y_l):
    
        '''
        lhs_val = 0
        rhs_val = 0
        
        lhs_val += sum(a.getPrimitiveTravelTime(x_l[a]) for a in self.network.links)
        lhs_val += sum(-r.intDinv(s, q_l[(r,s)], y_l[(r,s)]) for r in self.network.origins for s in r.getDests())
        
        rhs_val += sum(a.getPrimitiveTravelTime(x_f[a]) for a in self.network.links)
        rhs_val += sum(-r.intDinv(s, q_f[(r,s)], y_l[(r,s)]) for r in self.network.origins for s in r.getDests())
        
        print("compare ", lhs_val, rhs_val)
        '''
        
        approx = 0
        actual = 0
        
        '''
        for a in self.network.links:
            actual += a.getPrimitiveTravelTime(x_f[a])
            oa_point = x_f[a] * 0.8
            approx += a.getPrimitiveTravelTime(oa_point) + (x_f[a] - oa_point) * a.getTravelTime(oa_point, self.network.type)
        '''
        
        '''
        for r in self.network.origins:
            for s in r.getDests():
                actual += (-r.intDinv(s, q_f[(r,s)], y_l[(r,s)]))
                
                approx_y = y_l[(r,s)]*1.2
                approx_q = q_f[(r,s)]
                
                rho_oa = -r.intDinv(s, approx_q, approx_y)
                rho_oa -= (q_f[(r,s)] - approx_q) * r.Dinv(s, approx_q, approx_y)
                rho_oa -= (y_l[(r,s)] - approx_y) * r.intDerivDinv(s, approx_q, approx_y)
                approx += rho_oa
        '''
        
        '''
        for r in self.network.origins:
            for s in r.getDests():
                actual += (-r.intDinv(s, q_f[(r,s)], y_l[(r,s)]))
                y_approx = y_l[(r,s)]
                
                theta = (y_approx - self.rmp.y_lb[(r,s)]) / (self.rmp.y_ub[(r,s)] - self.rmp.y_lb[(r,s)])
                approx += -(theta * r.intDinv(s, q_f[(r,s)], self.rmp.y_ub[(r,s)]) + (1-theta) * r.intDinv(s, q_f[(r,s)], self.rmp.y_lb[(r,s)]))
      
            
        print("check approx=", approx, "actual=", actual)
        '''
        
        # OA on LHS from x_l, q_l
        for a in self.network.links:
            self.rmp.add_constraint(self.rmp.beta[a] >= a.getPrimitiveTravelTime(x_l[a]) + (self.rmp.x[a] - x_l[a]) * a.getTravelTime(x_l[a], self.network.type))
        
        for r in self.network.origins:
            for s in r.getDests():
                rho_oa = -r.intDinv(s, q_l[(r,s)], y_l[(r,s)])
                rho_oa -= (self.rmp.q[(r,s)] - q_l[(r,s)]) * r.Dinv(s, q_l[(r,s)], y_l[(r,s)])
                rho_oa -= (self.rmp.y[(r,s)] - y_l[(r,s)]) * r.intDerivDinv(s, q_l[(r,s)], y_l[(r,s)])
                self.rmp.add_constraint(self.rmp.rho[(r,s)] >= rho_oa)
            
        self.vfcut_x.append(x_f)
        self.vfcut_q.append(q_f)

        rhs = self.calcVFRHS(len(self.vfcut_x)-1)
        self.rmp.vfcut.append(self.rmp.add_constraint(sum(self.rmp.beta[a] for a in self.network.links) + sum(self.rmp.rho[(r,s)] for r in self.network.origins for s in r.getDests()) <= rhs))
      
        '''
        for a in self.network.links:
            self.rmp.add_constraint(self.rmp.x[a] == x_f[a])
            
        for r in self.network.origins:
            for s in r.getDests():
                self.rmp.add_constraint(self.rmp.q[(r,s)] >= q_f[(r,s)] - 0.1)
                self.rmp.add_constraint(self.rmp.q[(r,s)] <= q_f[(r,s)] + 0.1)
        '''



        '''
        lhsc = 0
        rhsc = 0
        
        for a in self.network.links:
            lhsc += a.getPrimitiveTravelTime(x_l[a]) + (x_f[a] - x_l[a]) * a.getTravelTime(x_l[a], self.network.type)
        
        for r in self.network.origins:
            for s in r.getDests():
                rho_oa = -r.intDinv(s, q_l[(r,s)], y_l[(r,s)])
                rho_oa -= (q_f[(r,s)] - q_l[(r,s)]) * r.Dinv(s, q_l[(r,s)], y_l[(r,s)])
                rho_oa -= (y_l[(r,s)] - y_l[(r,s)]) * r.intDerivDinv(s, q_l[(r,s)], y_l[(r,s)])
                lhsc += rho_oa
        
        for r in self.network.origins:
            for s in r.getDests():        
                theta = (y_l[(r,s)] - self.rmp.y_lb[(r,s)]) / (self.rmp.y_ub[(r,s)] - self.rmp.y_lb[(r,s)])
                rhsc += -theta * r.intDinv(s, q_f[(r,s)], self.rmp.y_ub[(r,s)])
                rhsc += -(1-theta) * r.intDinv(s, q_f[(r,s)], self.rmp.y_lb[(r,s)])
        
        rhsc += sum(a.getPrimitiveTravelTime(x_f[a]) for a in self.network.links)         
        
        print("vfcut", lhsc, rhsc)
        '''
      
    def addObjCut(self, xhat, qhat):
        for a in self.network.links:
            self.rmp.add_constraint(self.rmp.mu_a[a] >= (xhat[a]-self.x_target[a]) * (xhat[a]-self.x_target[a]) + (self.rmp.x[a] - xhat[a])* 2*(xhat[a]-self.x_target[a]))
        
    
        for r in self.network.origins:
            for s in r.getDests():
                self.rmp.add_constraint(self.rmp.mu_w[(r,s)] >= (qhat[(r,s)]-self.q_target[(r,s)]) * (qhat[(r,s)]-self.q_target[(r,s)]) + (self.rmp.q[(r,s)] - qhat[(r,s)])* 2*(qhat[(r,s)]-self.q_target[(r,s)]))
    
    def updateYbounds(self):
        for r in self.network.origins:
            for s in r.getDests():
                self.rmp.y_theta.rhs = self.rmp.theta[(r,s)] * (self.rmp.y_ub[(r,s)] - self.rmp.y_lb[(r,s)]) + self.rmp.y_lb[(r,s)]
        
        for idx in range(0, len(self.rmp.vfcut)):
            self.rmp.vfcut[idx].rhs = self.calcVFRHS()
    
    def calcVFRHS(self, idx):
        newrhs = 0
        newrhs += sum(-self.rmp.theta[(r,s)] * r.intDinv(s, self.vfcut_q[idx][(r,s)], self.rmp.y_ub[(r,s)]) for r in self.network.origins for s in r.getDests())
        newrhs += sum(-(1-self.rmp.theta[(r,s)]) * r.intDinv(s, self.vfcut_q[idx][(r,s)], self.rmp.y_lb[(r,s)]) for r in self.network.origins for s in r.getDests())
        
        newrhs += sum(a.getPrimitiveTravelTime(self.vfcut_x[idx][a]) for a in self.network.links) 
        return newrhs
        #return 1e15

    def initRMP(self):   
        self.rmp = Model()
    
        self.rmp.y_lb = dict()
        self.rmp.y_ub = dict()

        for r in self.network.origins:
            for s in r.getDests():
                self.rmp.y_lb[(r,s)] = 1
                self.rmp.y_ub[(r,s)] = 10 # change this!

        
        self.rmp.mu_a = {a:self.rmp.continuous_var(lb=0) for a in self.network.links}
        self.rmp.mu_w = {(r,s):self.rmp.continuous_var(lb=0) for r in self.network.origins for s in r.getDests()}
        #self.rmp.eta = {a:self.rmp.continuous_var(lb=-1e10,ub=1e10) for a in self.network.links}

        self.rmp.x = {a:self.rmp.continuous_var(lb=0, ub=self.network.TD) for a in self.network.links}
        self.rmp.xc = {(a,r):self.rmp.continuous_var(lb=0, ub=r.totaldemand) for a in self.network.links for r in self.network.origins}
        self.rmp.q = {(r,s): self.rmp.continuous_var(lb=0, ub=10*r.getDemand(s)) for r in self.network.origins for s in r.getDests()}
        
        
        self.rmp.y = {(r,s): self.rmp.continuous_var(lb=self.rmp.y_lb[(r,s)], ub=self.rmp.y_ub[(r,s)]) for r in self.network.origins for s in r.getDests()}
        self.rmp.beta = {a:self.rmp.continuous_var(lb=0) for a in self.network.links}
        self.rmp.rho = {(r,s):self.rmp.continuous_var(lb=-1e10) for r in self.network.origins for s in r.getDests()}
        self.rmp.theta = {(r,s):self.rmp.continuous_var(lb=0,ub=1) for r in self.network.origins for s in r.getDests()}
        
        # self.rmp.y.ub = new_upper_bound
        # self.rmp.y.lb = new_lower_bound
        
        
        self.rmp.y_theta = dict()
        self.rmp.vfcut = list()
        
        
        
        
        qzero = {(r,s): 0 for r in self.network.origins for s in r.getDests()}
        xzero = {a:0 for a in self.network.links}
        
        self.addObjCut(xzero, qzero)
        
        for r in self.network.origins:
            for s in r.getDests():
                self.rmp.y_theta = self.rmp.add_constraint(self.rmp.y[(r,s)] == self.rmp.theta[(r,s)] * (self.rmp.y_ub[(r,s)] - self.rmp.y_lb[(r,s)]) + self.rmp.y_lb[(r,s)])
        
        
        for a in self.network.links:
            self.rmp.add_constraint(sum(self.rmp.xc[(a,r)] for r in self.network.origins) == self.rmp.x[a])
            self.rmp.add_constraint(self.rmp.beta[a] >= a.t_ff * self.rmp.x[a])
            
        for i in self.network.nodes:                    
            for r in self.network.origins:            

                if i.id == r.id:
                    dem = - sum(self.rmp.q[(r,s)] for s in r.getDests())                
                elif isinstance(i, type(r)) == True and i in r.getDests():
                    dem = self.rmp.q[(r,i)]
                else:
                    dem = 0

                self.rmp.add_constraint(sum(self.rmp.xc[(a,r)] for a in i.incoming) - sum(self.rmp.xc[(a,r)] for a in i.outgoing) == dem)
        
        
        
        self.rmp.minimize(sum(self.rmp.mu_a[a] for a in self.network.links) + sum(self.rmp.mu_w[(r,s)] for r in self.network.origins for s in r.getDests()))
        
    def solveRMP(self):
    
        t_solve = time.time()
        self.rmp.solve(log_output=False)
        t_solve = time.time() - t_solve
        
        x_l = {a:self.rmp.x[a].solution_value for a in self.network.links}
        q_l = {(r,s): self.rmp.q[(r,s)].solution_value for r in self.network.origins for s in r.getDests()}
        y_l = {(r,s): self.rmp.y[(r,s)].solution_value for r in self.network.origins for s in r.getDests()}
        obj_l = self.rmp.objective_value
        
    
        
        return x_l, q_l, y_l, obj_l
        
    def TAP(self, y):
    
        for r in self.network.origins:
            for s in r.getDests():
                r.y[s] = y[(r,s)]
            
        self.network.tapas("UE", None)
        xhat = {a:a.x for a in self.network.links}
        qhat = {(r,s): r.bush.demand[s] for r in self.network.origins for s in r.getDests()}
        obj_f = self.calcOFV()
        
        return xhat, qhat, obj_f
 