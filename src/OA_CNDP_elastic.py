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
        self.y_target = dict()
        
        scenario = "0"

        for line in open("data/"+self.network.name+"/linkflows_"+scenario+".txt", "r"):
            data = line.split()
            i = int(data[0])
            j = int(data[1])
            x = float(data[2])
            self.x_target[self.network.findLink(i, j)] = x
            
        for line in open("data/"+self.network.name+"/demand_"+scenario+".txt", "r"):
            data = line.split()
            r = int(data[0])
            s = int(data[1])
            q = float(data[2])
            self.q_target[(self.network.findNode(r), self.network.findNode(s))] = q   

        for line in open("data/"+self.network.name+"/demand_y_"+scenario+".txt", "r"):
            data = line.split()
            r = int(data[0])
            s = int(data[1])
            q = float(data[2])
            self.y_target[(self.network.findNode(r), self.network.findNode(s))] = q    
        
        
 
        
        self.network.params.equilibrate_demand = True
        
        

        self.vfcut_x = list()
        self.vfcut_q = list()
        self.vfcuts = list()
        
        self.oacut_x = dict()
        self.oacut_q = dict()
        self.oacut_y = dict()
        
        for a in self.network.links:
            self.oacut_x[a] = list()
            
        for r in self.network.origins:
            for s in r.getDests():
                self.oacut_q[(r,s)] = list()
                self.oacut_y[(r,s)] = list()
        
        
        
        self.ub = 1e15
        self.best_q = dict()
        self.best_x = dict()
        self.best_y = dict()
    
    def solve(self):
        timelimit = 3600
        starttime = time.time()
        
        
        self.initRMP()
        
        '''
        for r in self.network.origins:
            for s in r.getDests():
                self.rmp.add_constraint(self.rmp.y[(r,s)] == 10.193556856721418)
        '''
        
        
        self.ub = 1e15
        lb = 0
        
        bb_nodes = []
        
        root = BB_node.BB_node(self.rmp.y_lb.copy(), self.rmp.y_ub.copy(), lb, 1e15)
        
        bb_nodes.append(root)
        
        max_node_iter = 10
        min_gap = 1e-2

        global_lb = 0

        max_iter = 100
        iter = 0
        
        print("iter", "global_lb", "best ub", "local_lb", "gap", "elapsed_time")
        
        while len(bb_nodes) > 0 and iter < max_iter:
        
            if self.network.params.PRINT_BB_INFO:
                print("\tavail nodes")
                for n in bb_nodes:
                    print("\t\t", n.lb, n.ub)

                    print("\t\t\t", n.y_lb)
                    print("\t\t\t", n.y_ub)
            
            
            
            # node selection logic
            idx = 0
            best_lb = 1e15
            best_ub = 1e15
            #best_ub = 0
            
            for n in range(0, len(bb_nodes)):
                if bb_nodes[n].lb < best_lb or (bb_nodes[n].lb == best_lb and bb_nodes[n].ub < best_ub):
                    best_lb = bb_nodes[n].lb
                    best_ub = bb_nodes[n].ub
                    idx = n
                    
            
            bb_node = bb_nodes.pop(idx)
            iter += 1
            
            if bb_node.lb > self.ub:
                if self.network.params.PRINT_BB_BASIC:
                    print(iter, "FATHOM")
                continue
            
            if self.network.params.PRINT_BB_INFO:
                print("--------------------")
            
            local_lb, local_ub = self.solveNode(bb_node, max_node_iter, timelimit, starttime)
            
            
            if local_lb is not None:
                global_lb = local_lb

            if len(bb_nodes) > 0:
                for n in bb_nodes:
                    global_lb = min(global_lb, n.lb)

            global_lb = max(0, global_lb) # numerical errors
            gap = min(1, self.ub)
            #gap = 1
            if global_lb > 0:
                gap = (self.ub - global_lb) / global_lb
            

            elapsed_time = time.time() - starttime

            print(iter, f"{global_lb:.3f}", f"{self.ub:.3f}", f"{local_lb:.3f}", f"{gap:.3f}", f"{elapsed_time:.2f}")
            
            if self.network.params.PRINT_BB_INFO:
                print("\tsolved node", bb_node.lb, local_lb, local_ub)
                print("\t\t", bb_node.y_lb)
                print("\t\t", bb_node.y_ub)
            
            '''
            if not self.validateFeasible():
                print("TARGET INFEAS")
                exit()
            '''
            
            if gap < min_gap:
                break

            if elapsed_time > timelimit:
                break

            if local_lb is not None and local_lb <= self.ub:
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
                
                #print("\tvf gap", total_gap, "worst", worst_gap, worst)    
        
                # y cannot be at ub or lb, or the gap would be 0...
                # branch on y directly will prevent further gap there
                
                # do not branch unless gap exists
                branch = True
                
                if branch:
                    y_lb_1 = bb_node.y_lb.copy()
                    y_lb_2 = bb_node.y_lb.copy()
                    y_ub_1 = bb_node.y_ub.copy()
                    y_ub_2 = bb_node.y_ub.copy()

                    if worst is None:
                        worst_gap = 0

                        for r in self.network.origins:
                            for s in r.getDests():
                                gap = bb_node.y_ub[(r,s)] - bb_node.y_lb[(r,s)] 

                                if gap > worst_gap:
                                    worst_gap = gap
                                    worst = (r,s)
                                    
                        mid = (bb_node.y_ub[worst] + bb_node.y_lb[worst])/2
                        y_lb_2[worst] = mid
                        y_ub_1[worst] = mid
                    else:
                        #mid = self.rmp.y[worst].solution_value
                        mid = (bb_node.y_ub[worst] + bb_node.y_lb[worst])/2
                        y_lb_2[worst] = mid
                        y_ub_1[worst] = mid
  
                        

                    left = BB_node.BB_node(y_lb_1, y_ub_1, local_lb, local_ub)
                    right = BB_node.BB_node(y_lb_2, y_ub_2, local_lb, local_ub)
                    
                    
                    # REMOVE THIS
                    '''
                    if self.checkTargetValidity(left.y_lb, left.y_ub):
                        bb_nodes.append(left)
                        
                        if self.network.params.PRINT_BB_INFO:
                            print("\tleft node", local_lb)
                            print("\t\t", left.y_lb)
                            print("\t\t", left.y_ub)
                        
                    if self.checkTargetValidity(right.y_lb, right.y_ub):
                        bb_nodes.append(right)
                    
                        if self.network.params.PRINT_BB_INFO:
                            print("\tright node", local_lb)
                            print("\t\t", right.y_lb)
                            print("\t\t", right.y_ub)
                    '''
                    
                    
                    
                    bb_nodes.append(left)
                        
                    if self.network.params.PRINT_BB_INFO:
                        print("\tleft node", local_lb)
                        print("\t\t", left.y_lb)
                        print("\t\t", left.y_ub)
                    
                    bb_nodes.append(right)
                    
                    if self.network.params.PRINT_BB_INFO:
                        print("\tright node", local_lb)
                        print("\t\t", right.y_lb)
                        print("\t\t", right.y_ub)
                    
            #max_node_iter = 1
        
        
        print("best obj", self.calcOFV(self.best_x, self.best_q))
        self.printSolution(self.best_x, self.best_q, self.best_y)
        
            
        
        
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
        last_x_l = {a:0 for a in self.network.links}
        last_y_l = {(r,s):0 for r in self.network.origins for s in r.getDests()}
        last_q_l = {(r,s):0 for r in self.network.origins for s in r.getDests()}
        
        last_gap = 1e14
        
        node_ub = 1e15
        
        while gap > cutoff and iteration < max_iter:
            iteration += 1
            
            # solve RMP -> y, LB
            #print("\t\tsolving RMP")
            status, x_l, q_l, y_l, obj_l = self.solveRMP()
            
            
            
            if status == 'infeasible':
                if self.network.params.PRINT_BB_BASIC:
                    print(status)
                    
                
                # REMOVE THIS
                valid = self.checkTargetValidity(self.rmp.y_lb, self.rmp.y_ub)
                            
                if valid:
                    print("INFEAS EXIT")
                    exit()
                
                return None, None

            lb = obj_l

            ll_l = self.calcLLobj(x_l, q_l, y_l)
            
            #print(y_l)
            
            # solve TAP -> x, UB
            #print("\t\tTAP")
            x_f, q_f, obj_f = self.TAP(y_l)
            ll_f = self.calcLLobj(x_f, q_f, y_l)
            
            
            #self.printSolution(x_l, q_l, y_l)


            self.addCuts(x_l, x_f, q_l, q_f, y_l)

            
         
            node_ub = min(node_ub, obj_f)
            
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
                print("\tBB node", iteration, lb, self.ub, f"{gap:.2f}", f"{elapsed:.2f}", f"{ll_l:.2f}", f"{ll_f:.2f}")
                
                '''
                if not self.validateFeasible():
                    print("TARGET INFEAS")
                    exit()
                '''
            
            #for a in self.varlinks:
            #    print("\t", a, yhat[a], a.C/2)
            
            if ll_l <= ll_f:
                break
 
            if gap < min_gap:
                break
                
            # lb is worse than best ub
            if lb > self.ub:
                break
                
            # no improvement due to weak vf cut
            #if gap == last_gap:
            #    break
            
            # if everything is the same, then adding another cut is pointless
            if self.isSameSolution(last_x_l, x_l, last_q_l, q_l, last_y_l, y_l):
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
            last_q_l = q_l
            last_gap = gap
            
        
        return lb, node_ub  
        
    def isSameSolution(self, last_x, x, last_q, q, last_y, y):
        for a in self.network.links:
            if abs(last_x[a] - x[a]) > self.network.params.SOL_TOL_X:
                return False
        
        for r in self.network.origins:
            for s in r.getDests():
                if abs(last_q[(r,s)] - q[(r,s)]) > self.network.params.SOL_TOL_Q or abs(last_y[(r,s)] - y[(r,s)]) > self.network.params.SOL_TOL_Y:
                    return False
                    
        return True
                    
    
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
        
    def addOApoint_x(self, a, x_l):
        return self.addOApoint_helper(self.oacut_x[a], x_l, self.network.params.OA_TOL_X)
    
    def addOApoint_q(self, r, s, q_l):
        return self.addOApoint_helper(self.oacut_q[(r,s)], q_l, self.network.params.OA_TOL_Q)
        
    def addOApoint_y(self, r, s, y_l):
        return self.addOApoint_helper(self.oacut_y[(r,s)], y_l, self.network.params.OA_TOL_Y)    
        
    def addOApoint_helper(self, list, pt, tol):
        # decide whether to add this point and return true or false
        similar = False
        
        # this is inefficient, I should change it!
        for p in list:
            if abs(pt - p) < tol:
                similar = True
                break
                
        if not similar:
            list.append(pt)
            return len(list)-1
        else:
            return -1
            
        
    def addOACut_x(self, a, x_l, idx):
        self.rmp.add_constraint(self.rmp.beta[a] >= a.getPrimitiveTravelTime(x_l) + (self.rmp.x[a] - x_l) * a.getTravelTime(x_l, self.network.type), ctname="oa_beta_"+str(idx)+"-"+str(a))
            
            
            
    def addOACut_q(self, r, s, q_l, y_l, idx):
        calc = -r.intDinv(s, q_l, y_l)
        calc += -(self.rmp.q[(r,s)] - q_l) * r.Dinv(s, q_l, y_l)
        calc += -(self.rmp.y[(r,s)] - y_l) * r.intDerivDinv(s, q_l, y_l)

        #print("check calc")
        #print("\t", q_l[(r,s)], y_l[(r,s)])
        #print("\t", r.intDinv(s, q_l[(r,s)], y_l[(r,s)]), r.Dinv(s, q_l[(r,s)], y_l[(r,s)]), r.intDerivDinv(s, q_l[(r,s)], y_l[(r,s)]))
        self.rmp.add_constraint(self.rmp.rho[(r,s)] >= calc, ctname="oa-rho_"+str(idx)+"-"+str(r)+"_"+str(s))
    
    def addObjCut_x(self, a, xhat, idx):
        self.rmp.add_constraint(self.rmp.mu_a[a] >= (xhat-self.x_target[a]) ** 2 + (self.rmp.x[a] - xhat)* 2*(xhat-self.x_target[a]), ctname="oa_mu_x_"+str(idx)+"-"+str(a))
        
    
    def addObjCut_q(self, r, s, qhat, idx):
        self.rmp.add_constraint(self.rmp.mu_w[(r,s)] >= (qhat-self.q_target[(r,s)]) * (qhat-self.q_target[(r,s)]) + (self.rmp.q[(r,s)] - qhat)* 2*(qhat-self.q_target[(r,s)]), ctname="oa_mu_q_"+str(idx)+"-"+str(r)+"_"+str(s))
     
        
    def addCuts(self, x_l, x_f, q_l, q_f, y_l):
        for a in self.network.links:
            idx = self.addOApoint_x(a, x_l[a])
            if idx >= 0:
                self.addObjCut_x(a, x_l[a], idx)
                self.addOACut_x(a, x_l[a], idx)
        
        for r in self.network.origins:
            for s in r.getDests():
                idx = self.addOApoint_q(r, s, q_l[(r,s)])
                if idx >= 0:
                    self.addObjCut_q(r, s, q_l[(r,s)], idx)
                    self.addOACut_q(r, s, q_l[(r,s)], y_l[(r,s)], idx)
                    
        self.vfcut_x.append(x_f)
        self.vfcut_q.append(q_f)


        rhs = self.calcVFRHS(x_f, q_f, self.rmp.theta)
        
        self.vfcuts.append(self.rmp.add_constraint(sum(self.rmp.beta[a] for a in self.network.links) + sum(self.rmp.rho[(r,s)] for r in self.network.origins for s in r.getDests()) <= rhs, ctname="vfcut_"+str(idx)))
      
        

    def checkVFcut(self):

        
        lhsc = 0
        rhsc = 1e15
        lhs_actual = 0
        rhs_actual = 1e15

        valid = self.checkTargetValidity(self.rmp.y_lb, self.rmp.y_ub)
        
        

        for a in self.network.links:
            beta_a = -1e15
            for x_oa in self.oacut_x:
                beta_a = max(beta_a, a.getPrimitiveTravelTime(x_oa[a]) + (self.x_target[a] - x_oa[a]) * a.getTravelTime(x_oa[a], self.network.type))
            lhsc += beta_a
        
        for r in self.network.origins:
            for s in r.getDests():
                rho_oa = -1e15
                for idx in range(0, len(self.oacut_q)):
                    q_oa = self.oacut_q[idx]
                    y_oa = self.oacut_y[idx]
                    
                    calc = -r.intDinv(s, q_oa[(r,s)], y_oa[(r,s)])
                    calc += -(self.q_target[(r,s)] - q_oa[(r,s)]) * r.Dinv(s, q_oa[(r,s)], y_oa[(r,s)])
                    calc += -(self.y_target[(r,s)] - y_oa[(r,s)]) * r.intDerivDinv(s, q_oa[(r,s)], y_oa[(r,s)])

                    rho_oa = max(calc, rho_oa)
                lhsc += rho_oa
                
        for a in self.network.links:
            lhs_actual += a.getPrimitiveTravelTime(self.x_target[a])

        for r in self.network.origins:
            for s in r.getDests():
                lhs_actual += -r.intDinv(s, self.q_target[(r,s)], self.y_target[(r,s)])
          
        
        for idx in range(0, len(self.vfcut_q)):  
            rhs_term = 0
            for a in self.network.links:
                rhs_term += a.getPrimitiveTravelTime(self.vfcut_x[idx][a])
            
            for r in self.network.origins:
                for s in r.getDests():
                    rhs_term += -r.intDinv(s, self.vfcut_q[idx][(r,s)], self.y_target[(r,s)])
            
            #print(rhs_term)    
            rhs_actual = min(rhs_actual, rhs_term)
           
           
        theta = dict()
        
        for r in self.network.origins:
            for s in r.getDests():   
                 theta[(r,s)] = (self.y_target[(r,s)] - self.rmp.y_lb[(r,s)]) / (self.rmp.y_ub[(r,s)] - self.rmp.y_lb[(r,s)])
                 
                 #print("CALC-theta", r, s, theta)
            
        for idx in range(0, len(self.vfcut_q)):
            rhs_term = 0
            
    
                  
                
            q_vf = self.vfcut_q[idx]
            x_vf = self.vfcut_x[idx]

            rhs_term = self.calcVFRHS(x_vf, q_vf, theta)
                    
            #print("CALC-RHS", idx, rhs_term)         
                
            rhsc = min(rhsc, rhs_term)
        
      
        if valid:
            #print("num vf cuts:", len(self.vfcut_x))
            #print("vfcut-check", lhsc, lhs_actual, rhs_actual, rhsc)

            if lhsc > lhs_actual or lhs_actual > rhs_actual or rhs_actual > rhsc:
                exit()
        
    def calcVFRHS(self, x_vf, q_vf, theta):
    
        rhs_term = 0
            
        for r in self.network.origins:
            for s in r.getDests():   


                #print(r, s, theta, self.rmp.y_lb[(r,s)], self.rmp.y_ub[(r,s)], self.y_target[(r,s)])

                side_l = -r.intDinv(s, q_vf[(r,s)], self.rmp.y_ub[(r,s)])
                side_r = -r.intDinv(s, q_vf[(r,s)], self.rmp.y_lb[(r,s)])

                rhs_term += theta[(r,s)]*side_l + (1-theta[(r,s)])*side_r

                #print(theta, theta*side_l + (1-theta)*side_r, side_l, side_r)

        beck = sum(a.getPrimitiveTravelTime(x_vf[a]) for a in self.network.links) 
        rhs_term += beck

        return rhs_term
        #return 1e15
    
    '''  
    def addObjCut(self, xhat, qhat):
        for a in self.network.links:
            self.rmp.add_constraint(self.rmp.mu_a[a] >= (xhat[a]-self.x_target[a]) ** 2 + (self.rmp.x[a] - xhat[a])* 2*(xhat[a]-self.x_target[a]))
        
    
        for r in self.network.origins:
            for s in r.getDests():
                self.rmp.add_constraint(self.rmp.mu_w[(r,s)] >= (qhat[(r,s)]-self.q_target[(r,s)]) * (qhat[(r,s)]-self.q_target[(r,s)]) + (self.rmp.q[(r,s)] - qhat[(r,s)])* 2*(qhat[(r,s)]-self.q_target[(r,s)]))
    '''
    def updateYbounds(self):
        for r in self.network.origins:
            for s in r.getDests():
                self.rmp.y[(r,s)].lb = self.rmp.y_lb[(r,s)]
                self.rmp.y[(r,s)].ub = self.rmp.y_ub[(r,s)]
                self.rmp.remove_constraint(self.rmp.y_theta[(r,s)])
                self.rmp.y_theta[(r,s)] = self.rmp.add_constraint(self.rmp.y[(r,s)] == self.rmp.theta[(r,s)] * self.rmp.y_ub[(r,s)] + (1-self.rmp.theta[(r,s)]) * self.rmp.y_lb[(r,s)], ctname="theta_"+str(r)+"_"+str(s))
                
        for idx in range(0, len(self.vfcuts)):
            self.rmp.remove_constraint(self.vfcuts[idx])
            rhs = self.calcVFRHS(self.vfcut_x[idx], self.vfcut_q[idx], self.rmp.theta)
            self.vfcuts[idx] = self.rmp.add_constraint(sum(self.rmp.beta[a] for a in self.network.links) + sum(self.rmp.rho[(r,s)] for r in self.network.origins for s in r.getDests()) <= rhs, ctname="vfcut-"+str(idx))
    
        #print("UPDATE y BOUNDS")
        
    
        
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
                
        '''  
        for r in self.network.origins:
            for s in r.getDests():
                if r.id == 1 and s.id == 2:
                    self.rmp.y_lb[(r,s)] = 9.172430041411648
                    self.rmp.y_ub[(r,s)] = 13.46110257294304
                elif r.id == 2 and s.id == 1:
                    self.rmp.y_lb[(r,s)] = 7.472245410284032
                    self.rmp.y_ub[(r,s)] = 10.497948507811676
        '''
        
        self.rmp.mu_a = {a:self.rmp.continuous_var(lb=0) for a in self.network.links}
        self.rmp.mu_w = {(r,s):self.rmp.continuous_var(lb=0) for r in self.network.origins for s in r.getDests()}
        #self.rmp.eta = {a:self.rmp.continuous_var(lb=-1e10,ub=1e10) for a in self.network.links}

        self.rmp.x = {a:self.rmp.continuous_var(lb=0) for a in self.network.links}
        self.rmp.xc = {(a,r):self.rmp.continuous_var(lb=0) for a in self.network.links for r in self.network.origins}
        self.rmp.q = {(r,s): self.rmp.continuous_var(lb=0, ub=10) for r in self.network.origins for s in r.getDests()}
        
        
        self.rmp.y = {(r,s): self.rmp.continuous_var(lb=self.rmp.y_lb[(r,s)], ub=self.rmp.y_ub[(r,s)]) for r in self.network.origins for s in r.getDests()}
        self.rmp.beta = {a:self.rmp.continuous_var(lb=0) for a in self.network.links}
        self.rmp.rho = {(r,s):self.rmp.continuous_var(lb=-1e15) for r in self.network.origins for s in r.getDests()}
        self.rmp.theta = {(r,s):self.rmp.continuous_var(lb=0,ub=1) for r in self.network.origins for s in r.getDests()}
        
        # self.rmp.y.ub = new_upper_bound
        # self.rmp.y.lb = new_lower_bound
        
        
        self.rmp.y_theta = dict()
        
        
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
        
        
        # avoid the x=0 solution
        for a in self.network.links:
            if self.addOApoint_x(a, 0) >= 0:
                self.addObjCut_x(a, 0, 0)
        
        for r in self.network.origins:
            for s in r.getDests():
                if self.addOApoint_q(r, s, 0) >= 0:
                    self.addObjCut_q(r, s, 0, 0)
        
        
        self.rmp.y_theta = dict()
        
        for r in self.network.origins:
            for s in r.getDests():
                self.rmp.y_theta[(r,s)] = self.rmp.add_constraint(self.rmp.y[(r,s)] == self.rmp.theta[(r,s)] * self.rmp.y_ub[(r,s)] + (1-self.rmp.theta[(r,s)]) * self.rmp.y_lb[(r,s)], ctname="theta_"+str(r)+"_"+str(s))
        
        
        for a in self.network.links:
            self.rmp.add_constraint(sum(self.rmp.xc[(a,r)] for r in self.network.origins) == self.rmp.x[a], ctname="xc-eq_"+str(a)+"_"+str(r))
            self.rmp.add_constraint(self.rmp.beta[a] >= a.t_ff * self.rmp.x[a], ctname="beta_"+str(a)+"0")
            
        for i in self.network.nodes:                    
            for r in self.network.origins:            

                if i.id == r.id:
                    dem = - sum(self.rmp.q[(r,s)] for s in r.getDests())                
                elif isinstance(i, type(r)) == True and i in r.getDests():
                    dem = self.rmp.q[(r,i)]
                else:
                    dem = 0
                    
                #print(i, r, dem)
                    
                self.rmp.add_constraint(sum(self.rmp.xc[(a,r)] for a in i.incoming) - sum(self.rmp.xc[(a,r)] for a in i.outgoing) == dem, ctname="cons_"+str(i)+"_"+str(r))
        
        
        self.rmp.objfunc = sum(self.rmp.mu_a[a] for a in self.network.links) + sum(self.rmp.mu_w[(r,s)] for r in self.network.origins for s in r.getDests())
        self.rmp.minimize(self.rmp.objfunc)
        
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
        
    def printSolution(self, x, q, y):
        print("link flows")
        for a in self.network.links:
            print("\t", a, x[a], self.x_target[a])
           
        print("demand")
        for r in self.network.origins:
            for s in r.getDests():
                print("\t", r, s, q[(r,s)], self.q_target[(r,s)], y[(r,s)])
 
    def checkTargetValidity(self, y_lb, y_ub):

        for r in self.network.origins:
            for s in r.getDests():   

                if self.y_target[(r,s)] < y_lb[(r,s)] or self.y_target[(r,s)] > y_ub[(r,s)]:
                    return False
        return True
        
    def validateFeasible(self):
        valid = self.checkTargetValidity(self.rmp.y_lb, self.rmp.y_ub)
        
        if not valid:
            return True
        
        solution = self.rmp.solution
      
        target_solution = self.rmp.new_solution()
        
        for a in self.network.links:
            target_solution.add_var_value(self.rmp.x[a], self.x_target[a]) 
            
            for r in self.network.origins:
                if self.x_target[a] == 5.118124 and r.id == 1:
                    target_solution.add_var_value(self.rmp.xc[(a,r)], self.x_target[a]) 
                elif self.x_target[a] == 5.223466 and r.id == 2:
                    target_solution.add_var_value(self.rmp.xc[(a,r)], self.x_target[a]) 
                else:
                    target_solution.add_var_value(self.rmp.xc[(a,r)], 0) 
                    
                    
            target_solution.add_var_value(self.rmp.beta[a], a.getPrimitiveTravelTime(self.x_target[a])+0.001)
            target_solution.add_var_value(self.rmp.mu_a[a], 0)
        
        for r in self.network.origins:
            for s in r.getDests():
                target_solution.add_var_value(self.rmp.q[(r,s)], self.q_target[(r,s)])
                
                target_solution.add_var_value(self.rmp.y[(r,s)], self.y_target[(r,s)])
                target_solution.add_var_value(self.rmp.rho[(r,s)], -r.intDinv(s, self.q_target[(r,s)], self.y_target[(r,s)]))
                target_solution.add_var_value(self.rmp.mu_w[(r,s)], 0)
                target_solution.add_var_value(self.rmp.theta[(r,s)], (self.y_target[(r,s)] - self.rmp.y[(r,s)].lb) / (self.rmp.y[(r,s)].ub - self.rmp.y[(r,s)].lb))
        
                #print("TARGET-THETA", r, s, target_solution.get_value(self.rmp.theta[(r,s)]))
        output = target_solution.is_feasible_solution()

        
        print("TEST", output, target_solution.get_value(self.rmp.objfunc))
        
        unsatisfied = target_solution.find_unsatisfied_constraints(self.rmp)
        
        for i in unsatisfied:
            print(i)
            print("\t", target_solution.get_value(i.lhs), target_solution.get_value(i.rhs))
            print("\t", solution.get_value(i.lhs), solution.get_value(i.rhs))
            
            
        if output == False:    
            self.checkVFcut()
            exit()
                
        return output
            
            
        
                