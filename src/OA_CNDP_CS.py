import time
from src import Params
from docplex.mp.model import Model

class OA_CNDP_CS:
    
    def __init__(self, network):
        self.network = network
        self.CG_tol = 1e-4
        self.inf = 1e+9
        self.cs_error = 0.01
        
        self.solveSO_only = True
        
        self.last_xhat = None
        self.last_obj_f = 1000000000
        
        self.sameycuts = []
        self.sameyvars = []
        
        self.g = {a:100*a.cost for a in self.network.links}
        
        for a in self.network.links2:
            a.y = 0
            a.add_cap = 0
            
        self.varlinks = []
        
        for a in self.network.links:
            if a.cost > 0:
                self.varlinks.append(a)
                
        print("yvars", len(self.varlinks))
        
        self.paths = {r:{s:[] for s in self.network.zones} for r in self.network.origins}
        self.link_cons = dict()
        self.dem_cons = dict()
        self.path_use_cons = dict()
        
        self.rt_pricing = 0.0
        self.params = Params.Params()
        self.pathcount = 0
        
    def solve(self):
        # initialize paths
        self.initRMP()
        
        
    
        iteration = 0
        best_ub = self.inf
        best_y = None
        
        start_time = time.time()
        
        while(True):
            iteration += 1
            
            print("solving RMP", iteration)
            
            RMP_status, ofv, x, y  = self.solveRMP()
            
            print(RMP_status)
            
            self.addTSTTCut(x, y)
            self.addLinkTTCut(x, y)
            unusedpaths, newpaths = self.addPaths(x, y)
            x_f, obj_f = self.TAP(y)
            
            if obj_f < best_ub:
                best_ub = obj_f
                best_y = y
                
            if iteration >= 10:
                break
            
            elapsed_time = time.time() - start_time    
            print(iteration, ofv, obj_f, best_ub, self.pathcount, unusedpaths, newpaths, elapsed_time)
        
    
        
    def addPaths(self, x, y):
        # search for new paths to add
        for a in self.network.links:
            a.x = x[a]
           
        for a in self.varlinks:
            a.add_cap = y[a]
            
        newpaths = False
        unusedpaths = False
        
        for r in self.network.origins:
            self.network.dijkstras(r, "UE")
            for s in self.network.zones:
                if r.getDemand(s) > 0:
                    min_tt_rs = s.cost
                    for p in self.paths[r][s]:
                        # we should be using this path, but we aren't
                        if self.rmp.phi[p].solution_value == 0 and p.getTravelTime("UE") <= min_tt_rs + self.cs_error:
                            unusedpaths = True
                        
                    if self.addPath(self.network.trace(r, s)):
                        newpaths = True
        
        return unusedpaths, newpaths
    
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
            
    
        
    
    def addTSTTCut(self, xhat, yhat):
        for a in self.network.links:
        

            if a in self.varlinks:
                firstterm = xhat[a] * a.getTravelTimeC(xhat[a], yhat[a], "UE")
                secondterm = a.getTravelTimeC(xhat[a], yhat[a], "UE") + xhat[a] * a.getDerivativeTravelTimeCx(xhat[a], yhat[a])
                thirdterm = xhat[a] * a.getDerivativeTravelTimeCy(xhat[a], yhat[a])
            
                self.rmp.add_constraint(self.rmp.zeta[a] >= firstterm + secondterm * (self.rmp.x[a] - xhat[a]) + thirdterm * (self.rmp.y[a] - yhat[a]))
            else:
                firstterm = xhat[a] * a.getTravelTimeC(xhat[a], 0, "UE")
                secondterm = a.getTravelTimeC(xhat[a], 0, "UE") + xhat[a] * a.getDerivativeTravelTimeCx(xhat[a], 0)

                self.rmp.add_constraint(self.rmp.zeta[a] >= firstterm + secondterm * (self.rmp.x[a] - xhat[a]))
        
    def addLinkTTCut(self, xhat, yhat):
        for a in self.network.links:
            yterm = 0
            
            y_ext = 0
            if a in self.varlinks:
                y_ext = yhat[a]
                yterm = (self.rmp.y[a] - y_ext) * a.getDerivativeTravelTimeCy(xhat[a], y_ext)
                
            txy = a.getTravelTimeC(xhat[a], y_ext, "UE")
            xterm = (self.rmp.x[a] - xhat[a]) * a.getDerivativeTravelTimeCx(xhat[a], y_ext)
            
            self.rmp.add_constraint(self.rmp.tau[a] >= txy + xterm + yterm)
      
    def initRMP(self):   
    
        # init paths
        
        
                    
                
        self.rmp = Model()
        
        
        self.rmp.parameters.read.scale = -1
        
        
        self.rmp.zeta = {a:self.rmp.continuous_var(lb=0,ub=1e10) for a in self.network.links}
        self.rmp.beta = {a:self.rmp.continuous_var(lb=0) for a in self.network.links}
        
        self.rmp.y = {a:self.rmp.continuous_var(lb=0, ub=a.max_add_cap) for a in self.varlinks}
        self.rmp.x = {a:self.rmp.continuous_var(lb=0, ub=self.network.TD) for a in self.network.links}
        self.rmp.h = dict()
        self.rmp.mu = {(r,s): self.rmp.continuous_var(lb=0) for r in self.network.origins for s in self.network.zones}
        self.rmp.phi = dict()
        self.rmp.theta = dict()
        self.rmp.tau = {a:self.rmp.continuous_var(lb=a.t_ff) for a in self.network.links}
        
        for r in self.network.origins:
            if r.totaldemand > 0:
                self.network.dijkstras(r, "UE")

                self.paths[r] = dict()
                
                for s in self.network.zones:
                    if r.getDemand(s) > 0:
                        p = self.network.trace(r, s)
                        self.paths[r][s] = list()
                        self.paths[r][s].append(p)
                        
                        self.createPathVars(p)
                        
                        self.path_use_cons[(r,s)] = self.rmp.add_constraint(self.rmp.phi[p] + 0 >= 1)
                        self.pathcount += 1
        
        
        
        # KKT >= 0
        for a in self.network.links:
            self.rmp.add_constraint(self.rmp.zeta[a] >= self.rmp.x[a] * a.t_ff)
        
    
        for a in self.network.links:
            self.link_cons[a] = self.rmp.add_constraint(self.rmp.x[a] - sum(self.rmp.h[p] for p in self.getPaths() if a in p.links) == 0)
            
        for r in self.network.origins:
            for s in self.network.zones:
                if r.getDemand(s) > 0:
                    self.dem_cons[(r,s)] = self.rmp.add_constraint(sum(self.rmp.h[p] for p in self.paths[r][s]) == r.getDemand(s))
        
        self.rmp.minimize(sum(self.rmp.zeta[a] for a in self.network.links) + sum(self.g[a] * self.rmp.y[a] for a in self.varlinks))
        
        
    def equals(self, p1, p2):
        if p1.r != p2.r or p1.s != p2.s or len(p1.links) != len(p2.links):
            return False
        
        for a in p1.links:
            if a not in p2.links:
                return False
                
        return True
        
    def addPath(self, p):
    
        # first make sure the path doesn't exist already
        for r in self.network.origins:
            for s in self.paths[r].keys():
                for p_old in self.paths[r][s]:
                    if self.equals(p, p_old):
                        return False
        
        self.pathcount += 1
        self.paths[p.r][p.s].append(p)
        
        self.createPathVars(p)
        
        
        self.dem_cons[(p.r, p.s)].lhs.add_term(self.rmp.h[p], 1)
        
        for a in p.links:
            self.link_cons[a].lhs.add_term(self.rmp.h[p], -1)
        
        self.path_use_cons[(r,s)].lhs.add_term(self.rmp.phi[p], 1)

        return True
    
    def createPathVars(self, p):
        self.rmp.h[p] = self.rmp.continuous_var(lb=0)
        self.rmp.phi[p] = self.rmp.binary_var()
        self.rmp.theta[p] = self.rmp.continuous_var(lb=0)

        self.rmp.add_constraint(self.rmp.theta[p] == sum(self.rmp.tau[a] for a in p.links))

        # KKT >=
        self.rmp.add_constraint(self.rmp.mu[(p.r,p.s)] <= self.cs_error + self.rmp.theta[p])
        
        # complementary slackness constraint
        self.rmp.add_constraint(self.rmp.h[p] <= p.r.getDemand(p.s) * self.rmp.phi[p])
        M = sum(a.getTravelTimeC(self.network.getTotalTrips(), 0, "UE") for a in p.links)
        self.rmp.add_constraint(self.rmp.theta[p] - self.rmp.mu[(p.r,p.s)] <= M * (1-self.rmp.phi[p]))
        
        
        
    def solveRMP(self):
    
        t_solve = time.time()
        self.rmp.solve(log_output=False)
        t_solve = time.time() - t_solve
        

        if self.rmp.solve_details.status == 'infeasible' or self.rmp.solve_details.status == 'integer infeasible':
            return 'infeasible',self.inf, dict(), dict()
        RMP_status = self.rmp.solve_details.status
        
        
        x = {a:self.rmp.x[a].solution_value for a in self.network.links}
        y = {a:self.rmp.y[a].solution_value for a in self.varlinks}
        #print(RMP_status)

        OFV = self.rmp.objective_value
        
        return RMP_status, OFV, x, y
        
    def isYDifferent(self, y, lasty):
        
        for a in self.varlinks:
            if abs(y[a] - lasty[a]) > 1e-6:
                return True
        return False
        
    def TAP(self, y):
    
        for a in self.varlinks:
            a.add_cap = y[a]
            
        self.network.tapas("UE", None)
        xhat = {a:a.x for a in self.network.links}
        obj_f = self.calcOFV()

        return xhat, obj_f
        
        
    def getPaths(self):
        all_paths = []
        for r in self.network.origins:
            for s in self.paths[r].keys():
                for p in self.paths[r][s]:
                    all_paths.append(p)
        return all_paths     
 