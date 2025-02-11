import time
from src import Params
from src import BB_node
from src import Network
from src import Zone
from src import Link
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
        self.vfcuts = list()
        
        self.ub = 1e15
        self.best_q = dict()
        self.best_x = dict()
        self.best_y = dict()
    
    def solve(self):
        timelimit = 3600
        starttime = time.time()
        
        
        self.initRMP()
        
        for r in self.network.origins:
            for s in r.getDests():
                self.rmp.add_constraint(self.rmp.y[(r,s)] == 10.193556856721418)
        
        self.ub = 1e15
        lb = 0
        
        bb_nodes = []
        
        root = BB_node.BB_node(self.rmp.y_lb.copy(), self.rmp.y_ub.copy(), lb)
        
        bb_nodes.append(root)
        
        max_node_iter = 10
        min_gap = -1e-2

        global_lb = 0

        max_iter = 1
        iter = 0
        
        while len(bb_nodes) > 0 and iter < max_iter:
        
            print("\tavail nodes")
            for n in bb_nodes:
                print("\t\t", n.lb)
                
            bb_node = bb_nodes.pop(0)
            iter += 1
            
            if bb_node.lb > self.ub:
                print(iter, "FATHOM")
                continue
            
            local_lb = self.solveNode(bb_node, max_node_iter, timelimit, starttime)
            
            
            if local_lb is not None:
                global_lb = local_lb

            if len(bb_nodes) > 0:
                for n in bb_nodes:
                    global_lb = min(global_lb, n.lb)

            gap = min(1, self.ub)
            if global_lb > 0:
                gap = (self.ub - global_lb) / global_lb
            

            elapsed_time = time.time() - starttime

            print(iter, global_lb, self.ub, local_lb, gap, elapsed_time)
            
            print("\tsolved node", bb_node.lb, local_lb)
            print("\t\t", bb_node.y_lb)
            print("\t\t", bb_node.y_ub)
            
            if gap < min_gap:
                break

            if elapsed_time > timelimit:
                break

            if local_lb is not None and local_lb < self.ub:
                # iterate through y, find largest gap...
                worst = None
                worst_gap = 0
                
                total_gap = 0

                for r in self.network.origins:
                    for s in r.getDests():
                        gap = self.calcVFgapPct(r, s)

                        total_gap += gap

                        if gap > worst_gap:
                            worst_gap = gap
                            worst = (r,s)

                # I missed the case where worst=0 but there is a gap between obj_l and obj_f due to weak vf cuts. 
                #In that case, probably want to add cuts!
                
                print("\tvf gap", total_gap, "worst", worst_gap, worst)    
        
                # y cannot be at ub or lb, or the gap would be 0...
                # branch on y directly will prevent further gap there
                
                # do not branch unless gap exists
                if worst is not None:
                    y_lb_1 = bb_node.y_lb.copy()
                    y_lb_2 = bb_node.y_lb.copy()
                    y_ub_1 = bb_node.y_ub.copy()
                    y_ub_2 = bb_node.y_ub.copy()

                    y_lb_2[worst] = self.rmp.y[worst].solution_value
                    y_ub_1[worst] = self.rmp.y[worst].solution_value

                    left = BB_node.BB_node(y_lb_1, y_ub_1, local_lb)
                    right = BB_node.BB_node(y_lb_2, y_ub_2, local_lb)
                    bb_nodes.append(left)
                    bb_nodes.append(right)
                    
                    print("\tleft node", local_lb)
                    print("\t\t", left.y_lb)
                    print("\t\t", left.y_ub)
                    
                    print("\tright node", local_lb)
                    print("\t\t", right.y_lb)
                    print("\t\t", right.y_ub)
                
            #max_node_iter = 1
        
        
        self.printSolution(self.best_x, self.best_q)
        
            
        
        
    def solveNode(self, bbnode, max_iter, timelimit, starttime):
       
        lb = bbnode.lb
        
        
       
        min_gap = 1e-2
            
        self.rmp.y_ub = bbnode.y_ub
        self.rmp.y_lb = bbnode.y_lb
        self.updateYbounds()
        
        iteration = 0
        
        
        gap = 1
        cutoff = 0.01
        
        last_x_f = {a:0 for a in self.network.links}
        last_y_l = {(r,s):0 for r in self.network.origins for s in r.getDests()}
        last_x_l = {a:0 for a in self.network.links}
        
        
        
        
        while gap > cutoff and iteration < max_iter:
            iteration += 1
            
            # solve RMP -> y, LB
            #print("\t\tsolving RMP")
            status, x_l, q_l, y_l, obj_l = self.solveRMP()
            
            if status == 'infeasible':
                if self.network.params.PRINT_BB_BASIC:
                    print(status)
                return None

            lb = obj_l

            ll_l = self.calcLLobj(x_l, q_l, y_l)
            
            #print(y_l)
            
            # solve TAP -> x, UB
            #print("\t\tTAP")
            x_f, q_f, obj_f = self.TAP(y_l)
            ll_f = self.calcLLobj(x_f, q_f, y_l)
            
            
            #self.printSolution(x_l, q_l)

            # add VF cut
            self.addVFCut(x_l, x_f, q_l, q_f, y_l)
            # add TSTT cut
            self.addObjCut(x_l, q_l)
            
         
            if self.ub > obj_f:
                self.ub = obj_f
                self.best_x = x_f
                self.best_q = q_f
                self.best_y = y_l
            
            elapsed = time.time() - starttime
            if lb > 0:
                gap = (self.ub - lb)/lb
            elif self.ub == 0:
                gap = 0
            else:
                gap = 1
            
            #print(obj_f)
            if self.network.params.PRINT_BB_BASIC:
                print("\tBB node", iteration, lb, self.ub, gap, elapsed, ll_l, ll_f)
            
            #for a in self.varlinks:
            #    print("\t", a, yhat[a], a.C/2)
            
 
            if gap < min_gap:
                break
                
            if elapsed > timelimit:
                break
                
            '''    
            delta_y = 0
            total_y = 0
            for r in self.network.origins:
                for s in r.getDests():
                    delta_y += abs(y_l[(r,s)] - last_y_l[(r,s)])
                    total_y += y_l[(r,s)]
                    
            
            
            delta_y_pct = delta_y / total_y
            print("\t\tdelta_y", delta_y, delta_y_pct)  
            
            total_gap = sum(self.calcVFgap(r, s) for r in self.network.origins for s in r.getDests())
            
            print("\t\tVF gap", total_gap) 
            '''

            last_x_l = x_l
            last_y_l = y_l
            last_x_l = x_l
            
        
        return lb    
    
    def calcOFV(self, x, q):
        output = sum((x[a] - self.x_target[a]) ** 2 for a in self.network.links) 
        output += sum( (q[(r,s)] - self.q_target[(r,s)]) ** 2 for r in self.network.origins for s in r.getDests())
        
        #if output == 0:
        #    self.printSolution(x, q)


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
        self.vfcuts.append(dict())
        idx = len(self.vfcuts)-1
        self.vfcuts[idx][(r,s)] = self.rmp.add_constraint(sum(self.rmp.beta[a] for a in self.network.links) + sum(self.rmp.rho[(r,s)] for r in self.network.origins for s in r.getDests()) <= rhs)
      
        
        
        



        
        lhsc = 0
        rhsc = 0
        
        for a in self.network.links:
            lhsc += a.getPrimitiveTravelTime(x_l[a]) + (self.x_target[a] - x_l[a]) * a.getTravelTime(x_l[a], self.network.type)
        
        for r in self.network.origins:
            for s in r.getDests():
                rho_oa = -r.intDinv(s, q_l[(r,s)], y_l[(r,s)])
                rho_oa -= (self.q_target[(r,s)] - q_l[(r,s)]) * r.Dinv(s, q_l[(r,s)], y_l[(r,s)])
                rho_oa -= (y_l[(r,s)] - y_l[(r,s)]) * r.intDerivDinv(s, q_l[(r,s)], y_l[(r,s)])
                lhsc += rho_oa
        
        for r in self.network.origins:
            for s in r.getDests():        
                theta = (y_l[(r,s)] - self.rmp.y_lb[(r,s)]) / (self.rmp.y_ub[(r,s)] - self.rmp.y_lb[(r,s)])
                rhsc += -theta * r.intDinv(s, self.q_target[(r,s)], self.rmp.y_ub[(r,s)])
                rhsc += -(1-theta) * r.intDinv(s, self.q_target[(r,s)], self.rmp.y_lb[(r,s)])
        
        rhsc += sum(a.getPrimitiveTravelTime(x_f[a]) for a in self.network.links)         
        
        print("vfcut-check", lhsc, rhsc)
        
        
    
      
    def addObjCut(self, xhat, qhat):
        for a in self.network.links:
            self.rmp.add_constraint(self.rmp.mu_a[a] >= (xhat[a]-self.x_target[a]) ** 2 + (self.rmp.x[a] - xhat[a])* 2*(xhat[a]-self.x_target[a]))
        
    
        for r in self.network.origins:
            for s in r.getDests():
                self.rmp.add_constraint(self.rmp.mu_w[(r,s)] >= (qhat[(r,s)]-self.q_target[(r,s)]) * (qhat[(r,s)]-self.q_target[(r,s)]) + (self.rmp.q[(r,s)] - qhat[(r,s)])* 2*(qhat[(r,s)]-self.q_target[(r,s)]))
    
    def updateYbounds(self):
        for r in self.network.origins:
            for s in r.getDests():
                self.rmp.y[(r,s)].lb = self.rmp.y_lb[(r,s)]
                self.rmp.y[(r,s)].ub = self.rmp.y_ub[(r,s)]
                self.rmp.remove_constraint(self.rmp.y_theta[(r,s)])
                self.rmp.y_theta[(r,s)] = self.rmp.add_constraint(self.rmp.y[(r,s)] == self.rmp.theta[(r,s)] * self.rmp.y_ub[(r,s)] + (1-self.rmp.theta[(r,s)]) * self.rmp.y_lb[(r,s)])
                
        for idx in range(0, len(self.rmp.vfcut)):
            self.rmp.remove_constraint(self.vfcuts[idx][(r,s)])
            rhs = self.calcVFRHS(idx)
            self.vfcuts[idx][(r,s)] = self.rmp.add_constraint(sum(self.rmp.beta[a] for a in self.network.links) + sum(self.rmp.rho[(r,s)] for r in self.network.origins for s in r.getDests()) <= rhs)
    
    def calcVFRHS(self, idx):
        newrhs = 0
        newrhs += sum(-self.rmp.theta[(r,s)] * r.intDinv(s, self.vfcut_q[idx][(r,s)], self.rmp.y_ub[(r,s)]) for r in self.network.origins for s in r.getDests())
        newrhs += sum(-(1-self.rmp.theta[(r,s)]) * r.intDinv(s, self.vfcut_q[idx][(r,s)], self.rmp.y_lb[(r,s)]) for r in self.network.origins for s in r.getDests())
        
        newrhs += sum(a.getPrimitiveTravelTime(self.vfcut_x[idx][a]) for a in self.network.links) 
        return newrhs
        #return 1e15
        
    def calcVFgap(self, r, s):
        theta = self.rmp.theta[(r,s)].solution_value
        
        output = 1e15
        
        for idx in range(0, len(self.vfcut_q)):
            actual = -r.intDinv(s, self.vfcut_q[idx][(r,s)], self.rmp.y[(r,s)].solution_value)
            approx = -(theta * r.intDinv(s, self.vfcut_q[idx][(r,s)], self.rmp.y_ub[(r,s)]) + (1-theta)*r.intDinv(s, self.vfcut_q[idx][(r,s)], self.rmp.y_lb[(r,s)]))
            diff = approx - actual
            
            output = min(output, diff)
            
        return output
        
    def calcVFgapPct(self, r, s):
        theta = self.rmp.theta[(r,s)].solution_value
        
        output = 1e15
        
        for idx in range(0, len(self.vfcut_q)):
            actual = -r.intDinv(s, self.vfcut_q[idx][(r,s)], self.rmp.y[(r,s)].solution_value)
            approx = -(theta * r.intDinv(s, self.vfcut_q[idx][(r,s)], self.rmp.y_ub[(r,s)]) + (1-theta)*r.intDinv(s, self.vfcut_q[idx][(r,s)], self.rmp.y_lb[(r,s)]))
            diff = approx - actual
            
            #print("\t\tgap", r, s, actual, approx, theta, self.rmp.y[(r,s)].solution_value, self.rmp.y[(r,s)].lb, self.rmp.y[(r,s)].ub)
            output = min(output, diff/abs(actual))
            
        return output

    def initRMP(self):   
        self.rmp = Model()
    
        self.rmp.y_lb = dict()
        self.rmp.y_ub = dict()
        
        self.rmp.parameters.read.scale = -1

        for r in self.network.origins:
            for s in r.getDests():
                self.rmp.y_lb[(r,s)] = 1
                self.rmp.y_ub[(r,s)] = 20 # change this!

        
        self.rmp.mu_a = {a:self.rmp.continuous_var(lb=0) for a in self.network.links}
        self.rmp.mu_w = {(r,s):self.rmp.continuous_var(lb=0) for r in self.network.origins for s in r.getDests()}
        #self.rmp.eta = {a:self.rmp.continuous_var(lb=-1e10,ub=1e10) for a in self.network.links}

        self.rmp.x = {a:self.rmp.continuous_var(lb=0) for a in self.network.links}
        self.rmp.xc = {(a,r):self.rmp.continuous_var(lb=0) for a in self.network.links for r in self.network.origins}
        self.rmp.q = {(r,s): self.rmp.continuous_var(lb=0) for r in self.network.origins for s in r.getDests()}
        
        
        self.rmp.y = {(r,s): self.rmp.continuous_var(lb=self.rmp.y_lb[(r,s)], ub=self.rmp.y_ub[(r,s)]) for r in self.network.origins for s in r.getDests()}
        self.rmp.beta = {a:self.rmp.continuous_var(lb=0) for a in self.network.links}
        self.rmp.rho = {(r,s):self.rmp.continuous_var(lb=-1e10) for r in self.network.origins for s in r.getDests()}
        self.rmp.theta = {(r,s):self.rmp.continuous_var(lb=0,ub=1) for r in self.network.origins for s in r.getDests()}
        
        # self.rmp.y.ub = new_upper_bound
        # self.rmp.y.lb = new_lower_bound
        
        
        self.rmp.y_theta = dict()
        self.rmp.vfcut = list()
        
        '''
        for a in self.network.links:
            self.rmp.add_constraint(self.rmp.x[a] >= self.x_target[a]-0.1)
            self.rmp.add_constraint(self.rmp.x[a] <= self.x_target[a]+0.1)
            print(a, self.x_target[a])
            
        for r in self.network.origins:
            for s in r.getDests():
                self.rmp.add_constraint(self.rmp.q[(r,s)] == self.q_target[(r,s)])
                #self.rmp.add_constraint(self.rmp.q[(r,s)] <= self.q_target[(r,s)] + 0.1)
                print(r, s, self.q_target[(r,s)])
        '''
        
        
        qzero = {(r,s): 0 for r in self.network.origins for s in r.getDests()}
        xzero = {a:0 for a in self.network.links}
        
        self.addObjCut(xzero, qzero)
        
        self.rmp.y_theta = dict()
        
        for r in self.network.origins:
            for s in r.getDests():
                self.rmp.y_theta[(r,s)] = self.rmp.add_constraint(self.rmp.y[(r,s)] == self.rmp.theta[(r,s)] * self.rmp.y_ub[(r,s)] + (1-self.rmp.theta[(r,s)]) * self.rmp.y_lb[(r,s)])
        
        
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
                    
                #print(i, r, dem)
                    
                self.rmp.add_constraint(sum(self.rmp.xc[(a,r)] for a in i.incoming) - sum(self.rmp.xc[(a,r)] for a in i.outgoing) == dem)
        
        
        
        self.rmp.minimize(sum(self.rmp.mu_a[a] for a in self.network.links) + sum(self.rmp.mu_w[(r,s)] for r in self.network.origins for s in r.getDests()))
        
    def solveRMP(self):
    
        t_solve = time.time()
        self.rmp.solve(log_output=False)
        t_solve = time.time() - t_solve
        
        if self.rmp.solve_details.status == 'infeasible' or self.rmp.solve_details.status == 'integer infeasible':
            return 'infeasible', dict(), dict(), dict(), 1e15
        
        if self.rmp.solve_details.status == 'optimal with unscaled infeasibilities':
            print(self.rmp.solve.details.status)
            return 'infeasible', dict(), dict(), dict(), 1e15

        
        x_l = {a:self.rmp.x[a].solution_value for a in self.network.links}
        q_l = {(r,s): self.rmp.q[(r,s)].solution_value for r in self.network.origins for s in r.getDests()}
        y_l = {(r,s): self.rmp.y[(r,s)].solution_value for r in self.network.origins for s in r.getDests()}
        obj_l = self.rmp.objective_value
        
    
        
        return "solved", x_l, q_l, y_l, obj_l
        
    def TAP(self, y):
    
        for r in self.network.origins:
            for s in r.getDests():
                r.y[s] = y[(r,s)]
            
        self.network.tapas("UE", None)
        xhat = {a:a.x for a in self.network.links}
        qhat = {(r,s): r.bush.demand[s] for r in self.network.origins for s in r.getDests()}
        obj_f = self.calcOFV(xhat, qhat)
        
        return xhat, qhat, obj_f
        
    def printSolution(self, x, q):
        print("link flows")
        for a in self.network.links:
            print("\t", a, x[a], self.x_target[a])
           
        print("demand")
        for r in self.network.origins:
            for s in r.getDests():
                print("\t", r, s, q[(r,s)], self.q_target[(r,s)])
 