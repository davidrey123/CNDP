import time
from src import Params
from src import BB_node
from src import Network
from src import Zone
from src import Link
from docplex.mp.model import Model
import math


class OA_CNDP_elastic_CG:
    
    def __init__(self, network, useCG):
        self.network = network
        
        self.useCG = useCG
        
        for a in self.network.links2:
            a.y = 0
            a.add_cap = 0
            
        self.varlinks = self.network.links
        
        self.x_target = dict()
        self.q_target = dict()
        self.y_target = dict()
        
        self.params = network.params
        
        scenario = "0"

        for line in open("data/"+self.network.name+"/linkflows_"+scenario+".txt", "r"):
            data = line.split()
            if len(data) > 0:
                i = int(data[0])
                j = int(data[1])
                x = float(data[2])
                self.x_target[self.network.findLink(i, j)] = x
            
        for line in open("data/"+self.network.name+"/demand_"+scenario+".txt", "r"):
            data = line.split()
            if len(data) > 0:
                r = int(data[0])
                s = int(data[1])
                q = float(data[2])
                self.q_target[(self.network.findNode(r), self.network.findNode(s))] = q   

        for line in open("data/"+self.network.name+"/demand_y_"+scenario+".txt", "r"):
            data = line.split()
            if len(data) > 0:
                r = int(data[0])
                s = int(data[1])
                y = float(data[2])
                self.y_target[(self.network.findNode(r), self.network.findNode(s))] = y    

        self.x_base = dict()
        self.q_base = dict()
        self.y_base = dict()
        self.xc_base = dict()
        
        if scenario != "0":
            scenario = scenario + "sol"
        
        for line in open("data/"+self.network.name+"/linkflows_"+scenario+".txt", "r"):
            data = line.split()
            if len(data) > 0:
                i = int(data[0])
                j = int(data[1])
                x = float(data[2])
                self.x_base[self.network.findLink(i, j)] = x
                
        for line in open("data/"+self.network.name+"/linkflowsC_"+scenario+".txt", "r"):
            data = line.split()
            if len(data) > 0:
                i = int(data[0])
                j = int(data[1])
                r = int(data[2])
                x = float(data[3])
                self.xc_base[(self.network.findLink(i, j), self.network.findNode(r))] = x
            
        for line in open("data/"+self.network.name+"/demand_"+scenario+".txt", "r"):
            data = line.split()
            if len(data) > 0:
                r = int(data[0])
                s = int(data[1])
                q = float(data[2])
                self.q_base[(self.network.findNode(r), self.network.findNode(s))] = q  

        for line in open("data/"+self.network.name+"/demand_y_"+scenario+".txt", "r"):
            data = line.split()
            if len(data) > 0:
                r = int(data[0])
                s = int(data[1])
                y = float(data[2])
                self.y_base[(self.network.findNode(r), self.network.findNode(s))] = y
 
        
        self.network.params.equilibrate_demand = True
        
        

        self.vfcut_x = list()
        self.vfcut_q = list()
        self.vfcuts = list()
        
        self.oacut_x = dict()
        self.oacut_q = dict()
        self.oacut_y = dict()
        
        self.q_lb = dict()
        self.q_ub = dict()
        self.mu_lb = dict()
        self.mu_ub = dict()
        
        for a in self.network.links:
            self.oacut_x[a] = list()
            
        for r in self.network.origins:
            for s in r.getDests():
                self.oacut_q[(r,s)] = list()
                self.oacut_y[(r,s)] = list()
        
        
        self.paths = {r:{s:[] for s in self.network.zones} for r in self.network.origins}
        
        
        
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
        
        root = BB_node.BB_node(self.rmp.y_lb.copy(), self.rmp.y_ub.copy(), lb, 1e15, self.mu_lb, self.mu_ub)
        
        
        bb_nodes.append(root)
        
        max_node_iter = 30
        min_gap = 1e-2

        global_lb = 0

        max_iter = 1
        iter = 0
        
        print("iter", "global_lb", "best ub", "local_lb", "gap", "elapsed_time")
        
        while len(bb_nodes) > 0 and iter < max_iter:
            
            
            
            
            
            # node selection logic
            idx = 0
            best_lb = 1e15
            largest_diff = 1e15
            #best_ub = 0
            
            # fathoming nodes with worse lb to avoid saving them
            new_bb_nodes = list()
            for n in range(0, len(bb_nodes)):
                if bb_nodes[n].lb < self.ub:
                    new_bb_nodes.append(bb_nodes[n])
                                    
                    diff = 0
                    
                    for r in self.network.origins:
                        for s in r.getDests():
                            diff = max(diff, bb_nodes[n].y_ub[(r,s)] - bb_nodes[n].y_lb[(r,s)])

                    
                    if bb_nodes[n].lb < best_lb or diff > largest_diff:
                        best_lb = bb_nodes[n].lb
                        largest_diff = diff
                        idx = len(new_bb_nodes)-1
                
                
                    
            bb_nodes = new_bb_nodes
            
            
            if self.network.params.PRINT_BB_INFO:
                print("\tavail nodes")
                for n in bb_nodes:
                    print("\t\t", n.lb, n.ub)

                    #print("\t\t\t", n.y_lb)
                    #print("\t\t\t", n.y_ub)
                    
                    
            bb_node = bb_nodes.pop(idx)
            
            
            iter += 1
            
            if bb_node.lb > self.ub:
                if self.network.params.PRINT_BB_BASIC:
                    print(iter, "FATHOM")
                continue
            
            if self.network.params.PRINT_BB_INFO:
                print("--------------------")
            
            print("solving node")
            status, local_lb, local_ub, local_y = self.solveNode(bb_node, max_node_iter, timelimit, starttime)


            if status == "infeasible":
                print(iter, "solved node --- fathom infeasible")
                continue
                 
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
                #print("\t\tlb", bb_node.y_lb)
                #print("\t\tub", bb_node.y_ub)
                #print("\t\tsol", local_y)

                
                if self.params.VALIDATE_BASE and not self.validateFeasible():
                    print("TARGET INFEAS")
                    exit()
                

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

                #I'm not sure this ever occurs
                '''
                if self.equalsY(local_y, bb_node.y_lb) or self.equalsY(local_y, bb_node.y_ub):
                    print("BRANCH check")
                '''

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

                        
                    
                    mid = 0
                    
                    if self.params.branching_strategy == self.params.BRANCH_MIDPOINT:
                        mid = (bb_node.y_ub[worst] + bb_node.y_lb[worst])/2
                    elif self.params.branching_strategy == self.params.BRANCH_SOL:
                        mid = self.rmp.y[worst].solution_value
                    elif self.params.branching_strategy == self.params.BRANCH_MAX_GAP:
                        largest_diff = 0
                        best_y = 0
                        
                        for q in self.oacut_q[worst]:
                            if q == 0:
                                continue
                            y_ub = bb_node.y_ub[worst]
                            y_lb = bb_node.y_lb[worst]
                            
                            
                            a = r.a[s]
                            
                            
                            #theta = (a * (y_ub-y_lb) / (a * math.log(y_ub/q) - a * math.log(y_lb/q)) - y_lb) / (y_ub-y_lb)
                            
                            theta = (-r.intDinv(s, q, y_ub) + r.intDinv(s, q, y_lb) - 2 * a * (y_ub-y_lb)/q ) / (2 * a * (y_ub - y_lb) * (y_ub - y_lb)/q)
                            
                            if theta < 0 or theta > 1:
                                theta = 0.5
                                
                            y = theta * y_ub + (1-theta) * y_lb
                            
                            diff = -(theta * r.intDinv(s, q, y_ub) + (1-theta) * r.intDinv(s, q, y_lb)) + r.intDinv(s, q, y)
                            
                            #print("best y check", theta, diff, y)
                            
                            
                            if diff > largest_diff:
                                largest_diff = diff
                                best_y = y
                        
                        mid = best_y
                    else:
                        print("UNKNOWN branch strategy")
                        exit()
                        
                    print("branching", worst, mid)
                    y_lb_2[worst] = mid
                    y_ub_1[worst] = mid

                    mu_ub1 = self.calcTTs(y_ub_1)
                    mu_lb_2 = self.calcTTs(y_lb_2)

                    left = BB_node.BB_node(y_lb_1, y_ub_1, local_lb, local_ub, bb_node.mu_lb, mu_ub1)
                    right = BB_node.BB_node(y_lb_2, y_ub_2, local_lb, local_ub, mu_lb_2, bb_node.mu_ub)


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
                        #print("\t\t", left.y_lb)
                        #print("\t\t", left.y_ub)

                    bb_nodes.append(right)

                    if self.network.params.PRINT_BB_INFO:
                        print("\tright node", local_lb)
                        #print("\t\t", right.y_lb)
                        #print("\t\t", right.y_ub)

                #max_node_iter = 1

        if self.network.params.PRINT_BB_INFO:
            print("best obj", self.calcOFV(self.best_x, self.best_q))
            #self.printSolution(self.best_x, self.best_q, self.best_y)
        
            
        
        
    def solveNode(self, bbnode, max_iter, timelimit, starttime):
       
        lb = bbnode.lb
        
        
        self.mu_lb = bbnode.mu_lb
        self.mu_ub = bbnode.mu_ub
       
        min_gap = 1e-2
            
        self.rmp.y_ub = bbnode.y_ub
        self.rmp.y_lb = bbnode.y_lb
        self.updateYbounds()
        
        best_ub = 1e15
        
        y_sol = None
        iteration = 0
        
        
        gap = 1
        cutoff = 0.01
        
        last_x_f = {a:0 for a in self.network.links}
        last_x_l = {a:0 for a in self.network.links}
        last_y_l = {(r,s):0 for r in self.network.origins for s in r.getDests()}
        last_q_l = {(r,s):0 for r in self.network.origins for s in r.getDests()}
        status = "infeasible"
        x_l = None
        q_l = None
        y_l = None
        obj_l = 1e15
        
        last_gap = 1e14
        
        node_ub = 1e15
        
        while gap > cutoff and iteration < max_iter:
            iteration += 1
            
            
            if self.params.VALIDATE_BASE and not self.validateFeasible():
                print("TARGET INFEAS")
                exit()
            
            
            # solve RMP -> y, LB
            #print("\t\tsolving RMP")
            
            if self.useCG:
                status, x_l, q_l, y_l, obj_l = self.CG()
            else:
                if self.params.PRINT_BB_INFO:
                    print("\tsolving RMP")
                status, x_l, q_l, y_l, obj_l = self.solveRMP()
            
            
            
            if status == 'infeasible':
                if self.network.params.PRINT_BB_BASIC:
                    print(status)
                    
                
                # REMOVE THIS
                '''
                valid = self.checkTargetValidity(self.rmp.y_lb, self.rmp.y_ub)
                            
                if valid:
                    print("INFEAS EXIT")
                    exit()
                '''
                
                return status, None, None, None
            
            y_test = self.network.calcY(x_l, q_l)
            
            for r in self.network.origins:
                for s in r.destSet:
                    if y_test[(r,s)] > self.rmp.y[(r,s)].ub:
                       y_test[(r,s)] =  self.rmp.y[(r,s)].ub
                    elif y_test[(r,s)] < self.rmp.y[(r,s)].lb:
                       y_test[(r,s)] =  self.rmp.y[(r,s)].lb
            
            
            ''' 
            for a in self.network.links:
                print(a, x_l[a], self.x_target[a])
            '''
            

            
               
            
            
            

            lb = obj_l

            ll_l = self.calcLLobj(x_l, q_l, y_l)
            
            #print(y_l)
            
            # solve TAP -> x, UB
            #print("\t\tTAP")
            x_f, q_f, obj_f = self.TAP(y_test)
            ll_f = self.calcLLobj(x_f, q_f, y_l)
            
            
            for r in self.network.origins:
                for s in r.destSet:
                    print( (r,s), q_l[(r,s)], q_f[(r,s)], self.q_target[(r,s)], y_test[(r,s)], self.y_target[(r,s)])
            
            
            #for a in self.network.links:
            #    print(a, x_f[a], self.x_base[a])
            
            '''
            y_test2 = self.network.calcY(self.x_target, self.q_target)
            
            for r in self.network.origins:
                for s in r.destSet:
                    print("test y", (r,s), y_test[(r,s)], y_test2[(r,s)], self.y_target[(r,s)], q_l[(r,s)], self.q_target[(r,s)])
            '''
            
            if self.params.PRINT_BB_INFO:
                print("rmp obj ", obj_l, self.calcOFV(x_l, q_l))
            
            #self.printSolution(x_l, q_l, y_l)


            self.addCuts(x_l, x_f, q_l, q_f, y_l)

            if self.params.PRINT_BB_INFO:
                print("obj f", obj_f, "ll-obj l", ll_l, "ll-obj f", ll_f)
            
         
            node_ub = min(node_ub, obj_f)
            
            if self.ub > obj_f:
                self.ub = obj_f
                self.best_x = x_f
                self.best_q = q_f
                self.best_y = y_test
                
                '''
                print("\n\n*****", obj_f)
                for r in self.network.origins:
                    for s in r.getDests():
                        print(r, s, r.y[s], r.a[s])
                self.printSolution(self.best_x, self.best_q, self.best_y)
                print("\n\n*****")
                '''
                
            if best_ub > obj_f:
                best_ub = obj_f
                y_sol = y_l
            
            elapsed = time.time() - starttime
            if lb > 0:
                gap = (self.ub - lb)/lb
            elif self.ub == 0:
                gap = 0
            else:
                gap = 1
            
            #print(obj_f)
            if self.network.params.PRINT_BB_BASIC:
                print("\tBB node", iteration, lb, obj_f, self.ub, f"{gap:.2f}", f"{elapsed:.2f}", f"{ll_l:.2f}", f"{ll_f:.2f}")
                
                if self.params.PRINT_BB_INFO:
                    #print("\t", q_l)
                    print("\t\tl", self.calcLLobj(x_l, q_l, self.rmp.y_lb), self.calcLLobj(x_l, q_l, self.rmp.y_ub))
                    #print("\t\tf", q_f)
                    print("\t\tdem gap", self.calcGap(x_l, q_l, y_l))
                
                '''
                for r in self.network.origins:
                    self.network.dijkstras(r, "FF")
                    for s in r.getDests():
                        print("\t\t", r, s, self.rmp.rho[(r,s)].solution_value, -s.cost * self.rmp.q[(r,s)].solution_value)
                '''
                
            
            #for a in self.varlinks:
            #    print("\t", a, yhat[a], a.C/2)
            
            
 
            if gap < min_gap:
                if self.params.PRINT_BB_INFO:
                    print("end by low gap", gap, self.ub, lb)
                break
                
            if ll_l <= ll_f:
                if self.params.PRINT_BB_INFO:
                    print("end by ll")
                break
                
            # lb is worse than best ub
            if lb > self.ub:
                if self.params.PRINT_BB_INFO:
                    print("end b/c lb > ub")
                break
                
            # no improvement due to weak vf cut
            #if gap == last_gap:
            #    break
            
            # if everything is the same, then adding another cut is pointless
            if self.isSameSolution(last_x_l, x_l, last_q_l, q_l, last_y_l, y_l):
                if self.params.PRINT_BB_INFO:
                    print("end due to same sol")
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
            
        
        return "solved", lb, node_ub, y_sol
        
        

        
    def equalsY(self, y_lhs, y_rhs):
        tol = 1e-4
        for w in y_lhs.keys():
            if abs(y_lhs[w] - y_rhs[w]) > tol:
                return False
        return True
        
    def calcGap(self, x, q, y):
    
        beck_diff = 0
        dem_diff = 0
        
        for a in self.network.links:
            a.x = x[a]
            beck_diff += abs(a.getPrimitiveTravelTime(x[a]) - self.rmp.beta[a].solution_value)
        
        for r in self.network.origins:
            for s in r.getDests():
                r.bush.demand[s] = q[(r,s)]
                r.y[s] = y[(r,s)]
                dem_diff += abs(self.rmp.rho[(r,s)].solution_value - (-r.intDinv(s, q[(r,s)], y[(r,s)])))
        
        sptt, tmf, totaldemand = self.network.getSPTT(self.network.type)
        
        rhs_val = sum(self.rmp.beta[a].solution_value for a in self.network.links) + sum(self.rmp.rho[(r,s)].solution_value for r in self.network.origins for s in r.getDests())
        
        return tmf/totaldemand, beck_diff, dem_diff, rhs_val
        
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
        total -= sum(r.intDinv(s, qhat[(r,s)], yhat[(r,s)]) for r in self.network.origins for s in r.getDests())  
        
        
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
            if abs(pt - p) < 1e-5: # was tol
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
            
            add = True
            for x_save in self.x_saved[a]:
                if abs(x_l[a] - x_save) < a.C * 0.001:
                    add = False
                    break
            
            if add:
                self.x_saved[a].append(x_l[a])
                idx = self.addOApoint_x(a, x_l[a])
                if idx >= 0:
                    self.addObjCut_x(a, x_l[a], idx)
                    self.addOACut_x(a, x_l[a], idx)
        
        for r in self.network.origins:
            for s in r.getDests():
                
                add = True
                for q_save in self.q_saved[(r,s)]:
                    if abs(q_save - q_l[(r,s)]) < 0.001:
                        add = False
                        break
                
                if add:
                    self.q_saved[(r,s)].append(q_l[(r,s)])
                    idx = self.addOApoint_q(r, s, q_l[(r,s)])
                    if idx >= 0:
                        self.addObjCut_q(r, s, q_l[(r,s)], idx)
                        self.addOACut_q(r, s, q_l[(r,s)], y_l[(r,s)], idx)
                    
        self.vfcut_x.append(x_f)
        self.vfcut_q.append(q_f)


        rhs = self.calcVFRHS(x_f, q_f, self.rmp.theta)
        
        self.vfcuts.append(self.rmp.add_constraint(sum(self.rmp.beta[a] for a in self.network.links) + sum(self.rmp.rho[(r,s)] for r in self.network.origins for s in r.getDests()) <= rhs, ctname="vfcut_"+str(self.cut_idx)))
        self.cut_idx += 1
        

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
        
        
        for r in self.network.origins:
            for s in r.getDests():
                if (r,s) in self.q_lb:
                    self.rmp.remove_constraint(self.q_lb[(r,s)])
                    self.rmp.remove_constraint(self.q_ub[(r,s)])

                    self.q_ub[(r,s)] = self.rmp.add_constraint(self.rmp.q[(r,s)] <= self.rmp.y[(r,s)] * math.sqrt(r.a[s] / self.mu_lb[(r,s)]), ctname="q_ub"+str((r,s)))
                    self.q_lb[(r,s)] = self.rmp.add_constraint(self.rmp.q[(r,s)] >= self.rmp.y[(r,s)] * math.sqrt(r.a[s] / self.mu_ub[(r,s)]), ctname="q_lb"+str((r,s)))

                    #print("q-lb", r, s, self.mu_lb[(r,s)], math.exp(-self.mu_lb[(r,s)] / r.a[s]), self.rmp.y_lb[(r,s)], math.exp(-self.mu_lb[(r,s)] / r.a[s]) *self.rmp.y_lb[(r,s)])
                    #print("q-ub", r, s, self.mu_ub[(r,s)], math.exp(-self.mu_ub[(r,s)] / r.a[s]), self.rmp.y_ub[(r,s)], math.exp(-self.mu_ub[(r,s)] / r.a[s]) *self.rmp.y_ub[(r,s)])
                    
                    #print("q-ub", r, s, r.a[s], self.mu_lb[(r,s)], math.exp(-self.mu_lb[(r,s)] / r.a[s]), self.rmp.y_ub[(r,s)], math.exp(-self.mu_lb[(r,s)] / r.a[s]) *self.y_base[(r,s)])
                    #print("q-lb", r, s, self.mu_ub[(r,s)], math.exp(-self.mu_ub[(r,s)] / r.a[s]), self.rmp.y_lb[(r,s)], math.exp(-self.mu_ub[(r,s)] / r.a[s]) *self.y_base[(r,s)])
                    
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
        
        
    def calcTTs(self, y):
        output = dict()
        
        for r in self.network.origins:
            for s in r.getDests():
                r.y[s] = y[(r,s)]
            
        self.network.tapas("UE", None)
        
        for r in self.network.origins:
            self.network.dijkstras(r, "UE")
            for s in r.getDests():
                output[(r,s)] = s.cost
                
        return output
        
    def initRMP(self):   
        self.rmp = Model()
        
        self.cut_idx = 0
        self.x_saved = dict()
        self.q_saved = dict()
        
        for a in self.network.links:
            self.x_saved[a] = list()

        for r in self.network.origins:
            for s in r.destSet:
                self.q_saved[(r,s)] = list()
    
        self.rmp.y_lb = dict()
        self.rmp.y_ub = dict()
        
        self.rmp.parameters.read.scale = -1

        for r in self.network.origins:
            for s in r.getDests():
                self.rmp.y_lb[(r,s)] = 1
                self.rmp.y_ub[(r,s)] = 1000 # change this!
        
        
        self.mu_lb = self.calcTTs(self.rmp.y_lb)
        self.mu_ub = self.calcTTs(self.rmp.y_ub)
                
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

        self.rmp.x = {a:self.rmp.continuous_var(lb=0, name="x_"+str(a)) for a in self.network.links}
        
        self.rmp.q = {(r,s): self.rmp.continuous_var(lb=0, ub=40) for r in self.network.origins for s in r.getDests()}
        
        print("CG", self.useCG)
        
        if self.useCG:
            
            # initialize path set
            for r in self.network.origins:
                self.network.dijkstras(r, "UE")
                for s in r.getDests():
                    p = self.network.trace(r, s)
                    self.paths[r][s] = list()
                    self.paths[r][s].append(p)

            self.rmp.h = {p:self.rmp.continuous_var(lb=0) for p in self.getPaths()}
            
            self.link_consLE = dict()
            self.link_consGE = dict()
            self.dem_consLE = dict()
            self.dem_consGE = dict()
            
        else:
            self.rmp.xc = {(a,r):self.rmp.continuous_var(lb=0, name="xc_"+str(a)+"_"+str(r)) for a in self.network.links for r in self.network.origins}
        
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
        
        
        for r in self.network.origins:
            for s in r.getDests():
                self.q_lb[(r,s)] = self.rmp.add_constraint(self.rmp.q[(r,s)] >= self.rmp.y[(r,s)] * math.exp(-self.mu_lb[(r,s)] / r.a[s]))
                self.q_ub[(r,s)] = self.rmp.add_constraint(self.rmp.q[(r,s)] <= self.rmp.y[(r,s)] * math.exp(-self.mu_ub[(r,s)] / r.a[s]))
        
        
        for a in self.network.links:
            self.rmp.add_constraint(self.rmp.beta[a] >= a.t_ff * self.rmp.x[a], ctname="beta_"+str(a)+"0")
        
        if self.useCG:
            
            for a in self.network.links:
                self.link_consGE[a] = self.rmp.add_constraint(self.rmp.x[a] - sum(self.rmp.h[p] for p in self.getPaths() if a in p.links) >= 0, ctname="link consGE "+str(a))
                self.link_consLE[a] = self.rmp.add_constraint(-self.rmp.x[a] + sum(self.rmp.h[p] for p in self.getPaths() if a in p.links) >= 0, ctname="link consLE "+str(a))
        
            for r in self.network.origins:
                for s in self.network.zones:
                    if r.getDemand(s) > 0:
                        self.dem_consGE[(r,s)] = self.rmp.add_constraint(sum(self.rmp.h[p] for p in self.paths[r][s]) -self.rmp.q[(r,s)] >= 0, ctname="dem consGE"+str((r,s)))
                        self.dem_consLE[(r,s)] = self.rmp.add_constraint(-sum(self.rmp.h[p] for p in self.paths[r][s]) +self.rmp.q[(r,s)] >= 0, ctname="dem consLE"+str((r,s)))
            
        else:
            for a in self.network.links:
                self.rmp.add_constraint(sum(self.rmp.xc[(a,r)] for r in self.network.origins) == self.rmp.x[a], ctname="xc-eq_"+str(a)+"_"+str(r))


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
            print("\t", a, x[a], self.x_target[a], a.getTravelTime(x[a], "UE"))
           
        print("demand")
        for r in self.network.origins:
            for s in r.getDests():
                print("\t", r, s, q[(r,s)], self.q_target[(r,s)], y[(r,s)])
 
    def checkTargetValidity(self):

        for r in self.network.origins:
            for s in r.getDests():   
                
                if self.y_base[(r,s)] < self.rmp.y[(r,s)].lb or self.y_base[(r,s)] > self.rmp.y[(r,s)].ub:
                    print("INVALID", (r,s), self.y_base[(r,s)], self.rmp.y[(r,s)].lb, self.rmp.y[(r,s)].ub)
                    return False
        return True
        
    def validateFeasible(self):
        print("VALIDATE FEASIBLE")
        valid = self.checkTargetValidity()
        
        if not valid:
            return True
        
        solution = self.rmp.solution
      
        target_solution = self.rmp.new_solution()
        
        for a in self.network.links:
            target_solution.add_var_value(self.rmp.x[a], self.x_base[a]) 
            target_solution.add_var_value(self.rmp.beta[a], a.getPrimitiveTravelTime(self.x_base[a])+0.001)
            target_solution.add_var_value(self.rmp.mu_a[a], (self.x_base[a] - self.x_target[a]) ** 2)
            
            
            for r in self.network.origins:
                target_solution.add_var_value(self.rmp.xc[(a,r)], self.xc_base[(a, r)]) 
                

                    

        
        for r in self.network.origins:
            for s in r.getDests():
                target_solution.add_var_value(self.rmp.q[(r,s)], self.q_base[(r,s)])
                
                target_solution.add_var_value(self.rmp.y[(r,s)], self.y_base[(r,s)])
                target_solution.add_var_value(self.rmp.rho[(r,s)], -r.intDinv(s, self.q_base[(r,s)], self.y_base[(r,s)]))
                target_solution.add_var_value(self.rmp.mu_w[(r,s)], (self.q_base[(r,s)] - self.q_target[(r,s)]) ** 2)
                target_solution.add_var_value(self.rmp.theta[(r,s)], (self.y_base[(r,s)] - self.rmp.y[(r,s)].lb) / (self.rmp.y[(r,s)].ub - self.rmp.y[(r,s)].lb))
        
                #print("TARGET-THETA", r, s, target_solution.get_value(self.rmp.theta[(r,s)]))
        output = target_solution.is_feasible_solution()

        
        print("TEST", output, "obj", target_solution.get_value(self.rmp.objfunc), "ll obj", self.calcLLobj(self.x_base, self.q_base, self.y_base))
        
        
        if output == False:
            unsatisfied = target_solution.find_unsatisfied_constraints(self.rmp)
        
            print("UNSATISFIED: ", len(unsatisfied))
            for i in unsatisfied:
                print(i)
                print("\t", target_solution.get_value(i.lhs), target_solution.get_value(i.rhs))
                #print("\t", solution.get_value(i.lhs), solution.get_value(i.rhs))
                
                
            for a in self.network.links:
                self.checkBounds(self.rmp.x[a], "x", target_solution)
                self.checkBounds(self.rmp.mu_a[a], "mu_a", target_solution)
                self.checkBounds(self.rmp.beta[a], "beta", target_solution)
                
                for r in self.network.origins:
                    self.checkBounds(self.rmp.xc[(a,r)], "xc", target_solution)
                    
            for r in self.network.origins:
                for s in r.destSet:
                    self.checkBounds(self.rmp.q[(r,s)], "q", target_solution)
                    self.checkBounds(self.rmp.y[(r,s)], "y", target_solution)
                    self.checkBounds(self.rmp.rho[(r,s)], "rho", target_solution)
                    self.checkBounds(self.rmp.mu_w[(r,s)], "mu_w", target_solution)
                    self.checkBounds(self.rmp.theta[(r,s)], "theta", target_solution)
            
        if output == False:    
            #self.checkVFcut()
            exit()
                
        return output
        
    def checkBounds(self, var, varname, target_solution):
        if target_solution.get_value(var) < var.lb or target_solution.get_value(var) > var.ub:
            print("exceeds bounds var", varname, target_solution.get_value(var), var.lb, var.ub)
        
    def CG(self):
        conv = False
        nCG = 0

        CG_status = ""
        x_l = None
        q_l = None
        y_l = None
        obj_l = 1e15
        
        while conv == False:
        
            '''
            if not self.validateFeasible():
                print("TARGET INFEAS")
                exit()
            '''
               
            
            RMP_status, x_l, q_l, y_l, obj_l = self.solveRMP()
            #print("solved?", RMP_status)
            
            
            
            
            if RMP_status == 'infeasible':
                CG_status = 'infeasible'
                
                print("CG failed")
                
                return CG_status, dict(), dict(), dict(), 1e15
            else:
                CG_status = "solved"
                
            link_duals = {a: self.link_consGE[a].dual_value - self.link_consLE[a].dual_value for a in self.network.links}
            
            for a in self.network.links:
                print(a, link_duals[a], self.link_consGE[a].dual_value, self.link_consLE[a].dual_value)
                
            dem_duals = {(r,s): self.dem_consGE[(r,s)].dual_value - self.dem_consLE[(r,s)].dual_value for r in self.network.origins for s in r.getDests()}
                
            minrc = self.pricing(link_duals, dem_duals)

            if self.params.PRINT_BB_INFO:
                npaths = len(self.getPaths())
                print('CG: %d\t%d\t%.3f\t%.2f' % (nCG,npaths,obj_l,minrc))
                
            if minrc > -abs(self.params.CG_tol):
                conv = True
                
            nCG += 1
        

 
        
        return CG_status, x_l, q_l, y_l, obj_l
        
    def pricing(self, link_duals, dem_duals):

        t0_pricing = time.time()
        
        new = 0
        minrc = 1e15
        
        
        #print("len check", len(link_duals))
        for a in self.network.links:
            a.dual = link_duals[a]
        
        for r in self.network.origins:
            self.network.dijkstras(r,'RC')
            
            for s in r.getDests():
                

                rc = - dem_duals[(r,s)] + s.cost                    

                if rc < - abs(self.params.CG_tol):
                    p = self.network.trace(r,s)
                    self.paths[r][s].append(p) #---is it needed to store paths if directly adding to RMP?

                    #---add new path var to RMP                                                
                    self.rmp.h[p] = self.rmp.continuous_var(lb=0)

                    #---update RMP constraints
                    self.dem_consGE[(r, s)].lhs.add_term(self.rmp.h[p], 1)
                    self.dem_consLE[(r, s)].lhs.add_term(self.rmp.h[p], -1)

                    for a in p.links:
                        self.link_consGE[a].lhs.add_term(self.rmp.h[p], -1)
                        self.link_consLE[a].lhs.add_term(self.rmp.h[p], 1)

                    new += 1

                if rc < minrc:
                    minrc = rc
                
        self.rt_pricing += (time.time() - t0_pricing)        
        return minrc
        
    def getPaths(self):
        all_paths = []
        for r in self.network.origins:
            for s in self.paths[r].keys():
                for p in self.paths[r][s]:
                    all_paths.append(p)
        return all_paths  
            
            
        
                