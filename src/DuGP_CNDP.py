import time
from src import Params
from docplex.mp.model import Model
import math



class DuGP_CNDP:
    
    def __init__(self, network):
        self.network = network
        
        self.inf = 1e+9
        self.solveSO_only = True
        
        self.g = {a:a.cost for a in self.network.links}
        
        for a in self.network.links2:
            a.y = 0
            a.add_cap = 0
            
        self.varlinks = []
        
        
        #for a in self.network.links:
        #    if a.cost > 0:
        #        self.varlinks.append(a)
                
        self.params = Params.Params()
        
        
        self.delta_a1 = {a:1 for a in self.network.links}
        self.delta_a2 = {a:1 for a in self.network.links}
        self.delta_a3 = {a:1 for a in self.varlinks}
                
    def calcBeckmann(self, xhat, yhat):
        total = 0
        
        for a in self.network.links:
            if a in self.varlinks:
                total += a.getPrimitiveTravelTimeC(xhat[a], yhat[a])
            else:
                total += a.getPrimitiveTravelTimeC(xhat[a], 0)
                
        return total
            
    def solve(self):
        timelimit = 3600
        iteration = 0
        starttime = time.time()
        ub = 1e100
        lb = 0
        
        last_lb = 0
        
        gap = 1
        cutoff = 0.01
        tap_time = 0
        
        last_xhat = {a:-1 for a in self.network.links}
        last_yhat = {a:-1 for a in self.varlinks}
        last_x_l = {a:-1 for a in self.network.links}
        best_y = {a:-1 for a in self.network.links}
        
        yhat = None
        xhat = None
        x_l = None
        
        B_f = 10000000
        obj_f = 100000000
        
        self.initRMP()
        
        while gap > cutoff:
            iteration += 1
            
            # solve RMP -> y, LB

            CG_status, obj_l, x_l, yhat = self.solveRMP()
            
            if CG_status == "infeasible":
                print('infeasible')
                break
            
            
            lb = obj_l
            B_l = self.calcBeckmann(x_l, yhat)
            # solve TAP -> x, UB

            
            
                
                
            #self.checkVFCut(x_l, yhat, x_l, xhat, yhat)
            # add TSTT cut
            if not self.solveSO_only: 
                # add VF cut
                if self.isYDifferent(yhat, last_yhat):
                    if self.params.PRINT_BB_INFO:
                        print("\tSolving TAP")
                    t1 = time.time()
                    
   
                    
                    xhat, obj_f = self.TAP(yhat, last_yhat)
                    

                    t1 = time.time()-t1
                    tap_time += t1
                    B_f = self.calcBeckmann(xhat, yhat)


                    if ub > obj_f:
                        ub = obj_f
                        best_y = yhat

                    #self.addVFCut(x_l, xhat, yhat)
                    self.addVFCut2(x_l, xhat, yhat)
                else:
                    if self.params.PRINT_BB_INFO:
                        print("\tSkipping TAP")
                    #self.addVFCutSameY(x_l, xhat, yhat)
                

                self.addTSTTCut(xhat, yhat)
            else:
                self.addTSTTCut(x_l, yhat, last_x_l)
            
            
            
            elapsed = time.time() - starttime
            if self.solveSO_only:
                gap = (lb - last_lb)/lb
                gap = 1
            else:
                if lb > 0:
                    gap = (ub - lb)/lb
                else:
                    gap = 1
            
            print(iteration, lb, ub, gap, elapsed, tap_time)
            
            if self.params.PRINT_BB_INFO:
                print("\t", B_l, B_f, B_l-B_f)
            
            #for a in self.varlinks:
            #    print("\t", a, yhat[a], last_yhat[a], best_y[a], a.C/2)
            
            if iteration > 2:
                break
                
            if elapsed > timelimit:
                break
                
            last_xhat = xhat
            last_yhat = yhat
            last_x_l = x_l
            last_lb = lb
        
        
        
        
    
    def initRMP(self):
        self.rmp = Model()

        self.rmp.logy = {a:self.rmp.continuous_var(lb=-1e100, ub=math.log(a.max_add_cap)) for a in self.varlinks}
        self.rmp.logx = {a:self.rmp.continuous_var(lb=-1e100, ub=math.log(self.network.TD)) for a in self.network.links}
        
        self.rmp.logZ = self.rmp.continuous_var()

        
        self.rmp.term1_lb = {a:1 for a in self.network.links}
        self.rmp.term2_lb = {a:1 for a in self.network.links}
        self.rmp.term3_lb = {a:1 for a in self.varlinks}
        
        eps_flow = 0.01
        eps_y = 0.0001
        
        for a in self.network.links:
            self.rmp.term1_lb[a] = a.t_ff * eps_flow
            self.rmp.term2_lb[a] = a.t_ff * a.alpha * pow(eps_flow / (a.C + eps_y), a.beta) * eps_flow
            
        for a in self.varlinks:
            self.rmp.term3_lb[a] = self.g[a] * eps_y
        
        
        self.rmp.logGamma1 = {a:self.rmp.continuous_var(lb=0) for a in self.network.links}
        self.rmp.logGamma2 = {a:self.rmp.continuous_var(lb=0) for a in self.network.links}
        self.rmp.logGamma3 = {a:self.rmp.continuous_var(lb=0) for a in self.varlinks}
        
        
        '''
        for a in self.network.links:
            self.rmp.add_constraint(sum(self.rmp.xc[(a,r)] for r in self.network.origins) == self.rmp.x[a])
            
        for i in self.network.nodes:                    
            for r in self.network.origins:            

                if i.id == r.id:
                    dem = - sum(r.getDemand(s) for s in self.network.zones)                
                elif isinstance(i, type(r)) == True:
                    dem = r.getDemand(i)
                else:
                    dem = 0

                self.rmp.add_constraint(sum(self.rmp.xc[(a,r)] for a in i.incoming) - sum(self.rmp.xc[(a,r)] for a in i.outgoing) == dem)
                
        self.rmp.minimize(self.rmp.logZ)  
        '''
        
        
    def solveRMP(self):
    
        t_solve = time.time()
        self.rmp.solve(log_output=False)
        t_solve = time.time() - t_solve
        
        '''
        if self.rmp.solve_details.status == "optimal with unscaled infeasibilities":
            print(self.rmp.solve_details.status, "solve again")
            self.rmp.parameters.read.scale = -1
            t_solve = time.time()
            self.rmp.solve(log_output=False)
            t_solve = time.time() - t_solve
        '''
        
        
        if self.rmp.solve_details.status == 'infeasible' or self.rmp.solve_details.status == 'integer infeasible':
            return 'infeasible',self.inf, dict(), dict()
        RMP_status = self.rmp.solve_details.status
        
        x = {a: math.exp(self.rmp.logx[a].solution_value) for a in self.network.links}
        y = {a: math.exp(self.rmp.logy[a].solution_value) for a in self.varlinks}
        
        OFV = self.calcOFV(x, y)
        
        print("\t",self.rmp.logZ.solution_value, math.exp(self.rmp.logZ.solution_value), OFV)
        
        return RMP_status, OFV, x, y 
    
    
    def calcOFV(self, x, y):
        output = 0
        
        for a in self.network.links:
            y_ext = 0
            if a in self.varlinks:
                y_ext = y[a]
                
            output += x[a] * a.getTravelTimeC(x[a], y_ext, "UE")
        
        for a in self.varlinks:
            output += self.g[a] * y[a]
        
        return output
        
        
    def addTSTTCut(self, xp, yp, lastx):
        # initial constraints
        
        
        
        lhs1chk = sum(self.delta_a1[a] * math.log(self.rmp.term1_lb[a]) + self.delta_a2[a] * math.log(self.rmp.term2_lb[a]) for a in self.network.links)
        lhs2chk = sum(self.delta_a1[a]*self.rmp.logGamma1[a].solution_value + self.delta_a2[a]*self.rmp.logGamma2[a].solution_value for a in self.network.links)
            
        print("logZ", self.rmp.logZ.solution_value, lhs1chk, lhs2chk)
        
        #for a in self.network.links:
        #    print("\tgamma1", a, self.rmp.logGamma1[a].solution_value, self.delta_a1[a], self.delta_a1[a] * math.log(self.rmp.gamma1_lb[a]))
        #    print("\t\tgamma2", self.rmp.logGamma2[a].solution_value, self.delta_a2[a], self.delta_a2[a] * math.log(self.rmp.gamma2_lb[a]))

        
        
        
        for a in self.network.links:
            y_ext = 0
            
            if a in self.varlinks:
                y_ext = yp[a]
                
            obj_p = self.calcOFV(xp, yp) + sum(self.g[a] * a.C for a in self.network.links)
            self.delta_a1[a] = (a.t_ff * xp[a]) / obj_p
            
            print("\tgamma1", xp[a], lastx[a], self.rmp.logGamma1[a].solution_value, math.log((a.t_ff * xp[a]) / self.delta_a1[a]))
            print("\t\tlb1", math.log(self.rmp.term1_lb[a]), math.log(self.delta_a1[a]))
            
            
            self.delta_a2[a] = (a.t_ff * a.alpha * pow( xp[a] / (a.C + y_ext), a.beta) * xp[a]) / obj_p
            
            delta_a3_ext = 0
            delta3_term = 0
            delta3_lb_term = 0
            
            if a in self.varlinks:
                self.delta_a3[a] = self.g[a] * (a.C + yp[a]) / obj_p
                delta_a3_ext = self.delta_a3[a]
                delta3_lb_term = delta_a3_ext * math.log(self.rmp.gamma3_lb[a])
                delta3_term = delta_a3_ext*self.rmp.logGamma3[a]
            
            lhs1 = sum(self.delta_a1[a] * math.log(self.rmp.term1_lb[a]) + self.delta_a2[a] * math.log(self.rmp.term2_lb[a]) + delta3_lb_term for a in self.network.links)
            lhs2 = sum(self.delta_a1[a]*self.rmp.logGamma1[a] + self.delta_a2[a]*self.rmp.logGamma2[a] + delta3_term for a in self.network.links)
            self.rmp.add_constraint(lhs1+lhs2-self.rmp.logZ <= 0)
            
            # new lower bounds on gamma
            
            #print(a, math.log(self.rmp.gamma1_lb[a])-math.log(self.delta_a1[a]), math.log(self.rmp.gamma2_lb[a])-math.log(self.delta_a2[a]))
            self.rmp.add_constraint(self.rmp.logGamma1[a] + math.log(self.delta_a1[a]) >= math.log(self.rmp.term1_lb[a]))
            self.rmp.add_constraint(self.rmp.logGamma2[a] + math.log(self.delta_a2[a]) >= math.log(self.rmp.term2_lb[a]))
            
            if a in self.varlinks:
                #print(a, self.delta_a3[a], self.rmp.gamma3_lb[a])
                self.rmp.add_constraint(self.rmp.logGamma3[a] + math.log(self.delta_a3[a]) >= math.log(self.rmp.term3_lb[a]))
                self.rmp.add_constraint(self.rmp.logGamma3[a] * self.delta_a3[a] == math.log(self.g[a]) + self.logy[a])
            
            self.rmp.add_constraint(self.rmp.logGamma1[a] * self.delta_a1[a] == self.rmp.logx[a] + math.log(a.t_ff))
        
        # loop and add constraint on violated condensation
        
    
    
    