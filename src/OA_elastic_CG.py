import time
from src import Params
from src import BB_node
from src import Network
from src import Zone
from src import Link
from docplex.mp.model import Model
import math


class OA_elastic_CG:
    
    def __init__(self, network, useCG, obj_weight, scenario):
        self.network = network
        
        self.useCG = useCG
        
        self.obj_weight = obj_weight
        
        for a in self.network.links2:
            a.y = 0
            a.add_cap = 0
            
        self.varlinks = self.network.links
        
        self.x_target = dict()
        self.q_target = dict()
        
    
        
        
        
        self.params = network.params
        

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



        self.x_base = dict()
        self.q_base = dict()
        self.xc_base = dict()
        
        scenario = "0"
        
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

       
        
        self.network.params.equilibrate_demand = True
        
        

        

        
        self.q_lb = dict()
        self.q_ub = dict()
        self.mu_lb = dict()
        self.mu_ub = dict()
        
        self.oacut_x = dict()
        self.oacut_eta = dict()
        self.oacut_q = dict()
        
        for a in self.network.links:
            self.oacut_x[a] = list()
            self.oacut_eta[a] = list()
            
        for r in self.network.origins:
            for s in r.destSet:
                self.oacut_q[(r,s)] = list()
        
        
        self.paths = {r:{s:[] for s in self.network.zones} for r in self.network.origins}
        
        
        
        self.ub = 1e15
        self.best_q = dict()
        self.best_x = dict()
        
        self.vfcut_q_ub = dict()
        self.vfcut_pi_ub = dict()
        
        
        
    
    def solve(self):
        timelimit = 14400
        starttime = time.time()
        
        
        self.initRMPFix()
        self.initRMP()
        
        

        
        self.ub = 1e15
        lb = 0
        
        bb_nodes = []
        
        
        
        root = BB_node.BB_node(self.q_lb.copy(), self.q_ub.copy(), lb, 1e15, 1e15, self.mu_lb, self.mu_ub)
        
        
        bb_nodes.append(root)
        
        max_node_iter = 20
        min_gap = 1e-2

        global_lb = 0
        global_ll_gap = 1e15

        max_iter = 100
        iter = 0
        
        print("iter", "global_lb", "best ub", "local_lb", "local_ub", "gap", "elapsed_time", "local ll gap", "global ll gap")
        
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
                            diff = max(diff, bb_nodes[n].q_ub[(r,s)] - bb_nodes[n].q_lb[(r,s)])

                    
                    if bb_nodes[n].lb < best_lb or (bb_nodes[n].lb == best_lb and diff > largest_diff):
                        best_lb = bb_nodes[n].lb
                        largest_diff = diff
                        idx = len(new_bb_nodes)-1
                
                
                    
            bb_nodes = new_bb_nodes
            
            
            if self.network.params.PRINT_BB_INFO:
                print("\tavail nodes")
                for n in bb_nodes:
                    print("\t\t", n.lb, n.ub)

                    #print("\t\t\t", n.q_lb)
                    #print("\t\t\t", n.q_ub)
                    
                    
            bb_node = bb_nodes.pop(idx)
            
            
            iter += 1
            
            if bb_node.lb > self.ub:
                if self.network.params.PRINT_BB_BASIC:
                    print(iter, "FATHOM")
                continue
            
            if self.network.params.PRINT_BB_INFO:
                print("--------------------")
            
                print("solving node", bb_node.lb)
                
            status, local_lb, local_ub, ll_gap = self.solveNode(bb_node, max_node_iter, timelimit, starttime)


            if status == "infeasible":
                if self.network.params.PRINT_BB_INFO:
                    print(iter, "solved node --- fathom infeasible")
                continue
                 
            global_lb = local_lb
            
            global_ll_gap = ll_gap

            if len(bb_nodes) > 0:
                for n in bb_nodes:
                    global_lb = min(global_lb, n.lb)
                    global_ll_gap = max(global_ll_gap, n.ll_gap)

            global_lb = max(0, global_lb) # numerical errors
            gap = min(1, self.ub)
            #gap = 1
            if global_lb > 0:
                gap = (self.ub - global_lb) / global_lb


            elapsed_time = time.time() - starttime

            print(iter, f"{global_lb:.3f}", f"{self.ub:.3f}", f"{local_lb:.3f}", f"{local_ub:.3f}", f"{gap:.3f}", f"{elapsed_time:.2f}", ll_gap, global_ll_gap)
            
            

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
                        gap = self.rmp.vf_ub[(r,s)].solution_value - self.rmp.q[(r,s)].solution_value * self.rmp.pi[(r,s)].solution_value

                        #print("best rs check", (r,s), gap)
                        total_gap += gap

                        if gap > worst_gap:
                            worst_gap = gap
                            worst = (r,s)


                if worst is None:
                    worst_gap = 0
                    
                    for r in self.network.origins:
                        for s in r.getDests():
                            gap = self.rmp.q[(r,s)].ub - self.rmp.q[(r,s)].lb 

                            if gap > worst_gap:
                                worst_gap = gap
                                worst = (r,s)
                    
                branch = True

                #if status == "end-ll":
                #    branch = False
               

                if branch:
                    
                    
                    q_lb_1 = bb_node.q_lb.copy()
                    q_lb_2 = bb_node.q_lb.copy()
                    q_ub_1 = bb_node.q_ub.copy()
                    q_ub_2 = bb_node.q_ub.copy()
                    
                    
                    mid = (q_ub_1[worst] + q_lb_1[worst]) / 2
                        
                    if self.params.PRINT_BB_INFO:
                        print("branching", worst, mid)
                    q_lb_2[worst] = mid
                    q_ub_1[worst] = mid

                    mu_ub_1 = self.calcMuBounds(q_ub_1)
                    mu_lb_2 = self.calcMuBounds(q_lb_2)

                    left = BB_node.BB_node(q_lb_1, q_ub_1, local_lb, ll_gap, local_ub, bb_node.mu_lb, mu_ub_1)
                    right = BB_node.BB_node(q_lb_2, q_ub_2, local_lb, ll_gap, local_ub, mu_lb_2, bb_node.mu_ub)




                    bb_nodes.append(left)

                    if self.network.params.PRINT_BB_INFO:
                        print("\tleft node", local_lb)
                        #print("\t\t", left.q_lb)
                        #print("\t\t", left.q_ub)

                    bb_nodes.append(right)

                    if self.network.params.PRINT_BB_INFO:
                        print("\tright node", local_lb)
                        #print("\t\t", right.q_lb)
                        #print("\t\t", right.q_ub)

                max_node_iter = 1

        if self.network.params.PRINT_BB_INFO:
            print("best obj", self.calcOFV(self.best_x, self.best_q))
            
        if self.params.PRINT_SOL:
            self.printSolution(self.best_x, self.best_q);
        
  
    def solveNode(self, bbnode, max_iter, timelimit, starttime):

       
        lb = bbnode.lb
        
        
        self.mu_lb = bbnode.mu_lb
        self.mu_ub = bbnode.mu_ub
       
        min_gap = 1e-2
            

        self.updateQbounds(bbnode.q_lb, bbnode.q_ub, bbnode.mu_lb, bbnode.mu_ub)
        
        best_ub = 1e15
        

        iteration = 0
        
        
        gap = 1
        cutoff = 0.01
        
        last_x_f = {a:0 for a in self.network.links}
        last_x_l = {a:0 for a in self.network.links}
        
        last_q_l = {(r,s):0 for r in self.network.origins for s in r.getDests()}
        status = "infeasible"
        x_l = None
        q_l = None

        obj_l = 1e15
        ll_gap = 1e15
        
        last_gap = 1e14
        
        node_ub = 1e15
        
        while gap > cutoff and iteration < max_iter:
            iteration += 1
            
            
            if self.params.VALIDATE_BASE and not self.validateFeasible():
                print("TARGET INFEAS")
                exit()
            
            
            # solve RMP -> y, LB
            
            
            #if self.params.PRINT_BB_INFO:
                #print("\tsolving RMP")
            status, x_l, q_l, obj_l = self.solveRMP()
            
            
            
            if status == 'infeasible':
                if self.network.params.PRINT_BB_BASIC:
                    print(status)

                return status, None, None, 1e15


            lb = obj_l

            
            ll_l = self.calcLLobj(x_l)
            #print("calc ll l", ll_l)
            


            x_f, obj_f = self.TAP(q_l)
            
            
            '''
            for r in self.network.origins:
                for s in r.destSet:
                    print("q", r, s, q_l[(r,s)])
            
            for a in self.network.links:
                print("x", a, x_l[a], x_f[a])
            '''
            
            ll_f = self.calcLLobj(x_f)
            #print("calc ll f", ll_f)
            
            
            
            if self.params.PRINT_BB_INFO:
                print("rmp obj ", obj_l, self.calcOFV(x_l, q_l), self.calcOFV(x_f, q_l), "from x", self.calcOFV_x(x_f), "from q", self.calcOFV_q(q_l))
                print("ll obj from rmp ", self.getRMP_ll_obj(), "rmp ll ub", self.getRMP_ll_ub(), "actual rmp ll", ll_l, "TAP", ll_f)
                #self.network.checkDualBeckmann()
            
            #self.printSolution(x_l, q_l, y_l)

            eta_l = {a : self.rmp.eta[a].solution_value for a in self.network.links}
            self.addCuts(x_l, x_f, q_l, eta_l)

            ll_gap = (ll_l - ll_f) / ll_f
            
            if ll_gap < self.params.ll_tol:
                print("\n\nSOLVE FIXED\n\n")
                ll_l2, q_l2 = self.solveFixed(x_f, bbnode.q_lb, bbnode.q_ub)
                x_f2, obj_f2 = self.TAP(q_l2)
                
                obj_l2 = self.calcOFV(x_f, q_l2)
                
                print("test 4 solve", obj_l, obj_l2, obj_f, obj_f2)
                # solve rmp with fixed l -> demand
                # solve tap with new demand
                
                for r in self.network.origins:
                    for s in r.getDests():
                        print((r,s), q_l[(r,s)], q_l2[(r,s)])
            
                exit(1)
            
            
            node_ub = min(node_ub, obj_f)
            
            if self.ub > obj_f:
                self.ub = obj_f
                self.best_x = x_f
                self.best_q = q_l
                
                
            if best_ub > obj_f:
                best_ub = obj_f
            
            elapsed = time.time() - starttime
            if lb > 0:
                gap = (self.ub - lb)/lb
            elif self.ub == 0:
                gap = 0
            else:
                gap = 1
                
            
            
            
            #print(obj_f)
            if self.network.params.PRINT_BB_BASIC:
                print("\tBB node", iteration, lb, obj_f, self.ub, f"{gap:.2f}", f"{elapsed:.2f}", f"{ll_l:.2f}", f"{ll_f:.2f}", ll_gap)
                
                if self.params.PRINT_BB_INFO:
                    #print("\t", q_l)
                    print("\t\tl", self.calcLLobj(x_l), self.calcLLobj(x_l))
                    #print("\t\tf", q_f)
                    print("\t\tdem gap", self.calcGap(x_l, q_l))
                    
                

            '''
            if ll_gap < self.params.ll_tol:
                if self.params.PRINT_BB_INFO:
                    print("end by low ll gap", gap, self.ub, lb)

                return "end-ll", lb, node_ub, ll_gap
            '''
 
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
            if self.isSameSolution(last_x_l, x_l, last_q_l, q_l):
                if self.params.PRINT_BB_INFO:
                    print("end due to same sol")
                break
                
            if elapsed > timelimit:
                break
                
          

            last_x_l = x_l
            last_q_l = q_l
            last_gap = gap
            
        
        return "solved", lb, node_ub, ll_gap
        
        
        
    def isSameSolution(self, last_x, x, last_q, q):
        for a in self.network.links:
            if abs(last_x[a] - x[a]) > self.network.params.SOL_TOL_X:
                return False
        
        for r in self.network.origins:
            for s in r.getDests():
                if abs(last_q[(r,s)] - q[(r,s)]) > self.network.params.SOL_TOL_Q:
                    return False
                    
        return True
        
    
        
    def getRMP_ll_ub(self):
        
        #print("eta gap", self.getEtaGap())
        
        total = 0
        
        gap = 0
        
        for r in self.network.origins:
            for s in r.destSet:
                
                
                #self.vfcut_q_ub[(r,s)].rhs = self.rmp.q[(r,s)].ub * self.rmp.pi[(r,s)] + self.rmp.q[(r,s)] * self.mu_lb[(r,s)] - self.rmp.q[(r,s)].ub * self.mu_lb[(r,s)]
                #self.vfcut_pi_ub[(r,s)].rhs = self.rmp.q[(r,s)] * self.mu_ub[(r,s)] + self.rmp.q[(r,s)].lb * self.rmp.pi[(r,s)] - self.rmp.q[(r,s)].lb * self.mu_ub[(r,s)]
                
                mccormick1 = self.rmp.q[(r,s)].ub * self.rmp.pi[(r,s)].solution_value + self.rmp.q[(r,s)].solution_value * self.mu_lb[(r,s)] - self.rmp.q[(r,s)].ub * self.mu_lb[(r,s)]
                mccormick2 = self.rmp.q[(r,s)].solution_value * self.mu_ub[(r,s)] + self.rmp.q[(r,s)].lb * self.rmp.pi[(r,s)].solution_value - self.rmp.q[(r,s)].lb * self.mu_ub[(r,s)]
                #print("\t", "mccormick ub", (r,s), self.rmp.q[(r,s)].lb, self.rmp.q[(r,s)].ub, self.mu_lb[(r,s)], self.mu_ub[(r,s)])
                #print("\t\tbound 1", mccormick1, self.vfcut_q_ub[(r,s)].rhs.solution_value)
                #print("\t\tbound 2", mccormick2, self.vfcut_pi_ub[(r,s)].rhs.solution_value)
                
                
                #print("\tmccormick bound vs value", min(mccormick1, mccormick2), self.rmp.q[(r,s)].solution_value * self.rmp.pi[(r,s)].solution_value)
        
                #total += self.rmp.q[(r,s)].solution_value * self.rmp.pi[(r,s)].solution_value
                total += min(mccormick1, mccormick2)
                
                gap += min(mccormick1, mccormick2) - self.rmp.q[(r,s)].solution_value * self.rmp.pi[(r,s)].solution_value

        
        subtract = sum(self.rmp.eta_oa[a].solution_value for a in self.network.links)
        actual = 0
        
        for a in self.network.links:
            g = a.getConst()
            
            p = a.beta
            ge = pow(g, 1/p)
            eta = self.rmp.eta[a].solution_value
            eta_term = p / ((p+1) * ge) * pow(eta, (p+1)/p)
            actual += eta_term
            
            #print(a, eta_term, self.rmp.eta_oa[a].solution_value)
        
        
        print("sum eta", subtract, actual)
        print("mccormick gap", gap)
        
        return (1+self.params.ub_eps) * total - subtract
        
          
    def getRMP_ll_obj(self):
        return sum(self.rmp.beta[a].solution_value for a in self.network.links)
        
   
        
    def calcGap(self, x, q):
    
        beck_diff = 0
        
        for a in self.network.links:
            a.x = x[a]
            beck_diff += abs(a.getPrimitiveTravelTime(x[a]) - self.rmp.beta[a].solution_value)
        

        return beck_diff
        

    def calcOFV_x(self, x):
        return (1-self.obj_weight) * sum((x[a] - self.x_target[a]) ** 2 for a in self.network.links if a in self.x_target) 
        
    def calcOFV_q(self, q):
        return self.obj_weight * sum( (q[(r,s)] - self.q_target[(r,s)]) ** 2 for r in self.network.origins for s in r.getDests() if (r,s) in self.q_target)


    
    def calcOFV(self, x, q):
        return self.calcOFV_x(x) + self.calcOFV_q(q)

    def calcLLobj(self, xhat):
        total = 0
        
        total = sum(a.getPrimitiveTravelTimeC(xhat[a], 0) for a in self.network.links)
        
        return total
    
    def addOApoint_x(self, a, x_l):
        return self.addOApoint_helper(self.oacut_x[a], x_l, self.network.params.OA_TOL_X)
        
    def addOApoint_eta(self, a, eta_l):
        return self.addOApoint_helper(self.oacut_eta[a], eta_l, self.network.params.OA_TOL_X)
    
    def addOApoint_q(self, r, s, q_l):
        return self.addOApoint_helper(self.oacut_q[(r,s)], q_l, self.network.params.OA_TOL_Q)
        
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
    
    def addObjCut_x(self, a, xhat, idx):
        self.rmp.add_constraint(self.rmp.mu_a[a] >= (xhat-self.x_target[a]) ** 2 + (self.rmp.x[a] - xhat)* 2*(xhat-self.x_target[a]), ctname="oa_mu_x_"+str(idx)+"-"+str(a))
        
    
    def addObjCut_q(self, r, s, qhat, idx):
        self.rmp.add_constraint(self.rmp.mu_w[(r,s)] >= (qhat-self.q_target[(r,s)]) * (qhat-self.q_target[(r,s)]) + (self.rmp.q[(r,s)] - qhat)* 2*(qhat-self.q_target[(r,s)]), ctname="oa_mu_q_"+str(idx)+"-"+str(r)+"_"+str(s))
        self.rmp_fix.add_constraint(self.rmp_fix.mu_w[(r,s)] >= (qhat-self.q_target[(r,s)]) * (qhat-self.q_target[(r,s)]) + (self.rmp_fix.q[(r,s)] - qhat)* 2*(qhat-self.q_target[(r,s)]), ctname="oa_mu_q_"+str(idx)+"-"+str(r)+"_"+str(s))
     
    def addOACut_x(self, a, x_l, idx):
        self.rmp.add_constraint(self.rmp.beta[a] >= a.getPrimitiveTravelTime(x_l) + (self.rmp.x[a] - x_l) * a.getTravelTime(x_l, self.network.type), ctname="oa_beta_"+str(idx)+"-"+str(a))
        
    def addOACut_eta(self, a, eta_p):
        g = a.getConst()
        p = a.beta
        ge = pow(g, 1/p)
        self.rmp.add_constraint(self.rmp.eta_oa[a] >= p / ((p+1) * ge) * pow(eta_p, (p+1)/p) + pow(eta_p, (p+1)/p - 1) / ge * (self.rmp.eta[a] - eta_p))
        
    def addCuts(self, x_l, x_f, q_l, eta_l):
        for a in self.network.links:
            
           
            idx = self.addOApoint_x(a, x_l[a])
            if idx >= 0:
                self.addOACut_x(a, x_l[a], idx)
                
                
                if a in self.x_target:
                    self.addObjCut_x(a, x_l[a], idx)
        
        for r in self.network.origins:
            for s in r.getDests():
                if (r,s) in self.q_target:
                    idx = self.addOApoint_q(r, s, q_l[(r,s)])
                    if idx >= 0:
                        self.addObjCut_q(r, s, q_l[(r,s)], idx)

        for a in self.network.links:
            idx = self.addOApoint_eta(a, eta_l[a])
            
            if idx >= 0:
                self.addOACut_eta(a, eta_l[a])
        
    def calcMuBounds(self, q):
        output = dict()
        
        for r in self.network.origins:
            for s in r.getDests():
                r.demand[s] = q[(r,s)]
            
        self.network.tapas("UE", None)
        
        for r in self.network.origins:
            self.network.dijkstras(r, "UE")
            for s in r.destSet:
                path = self.network.trace(r, s)
                
                cost = 0
                for a in path.links:
                    cost += a.getTravelTime(a.x, "UE")
                output[(r,s)] = cost
                
        return output
        
    def calcTTs(self, q):
        output = dict()
        
        for r in self.network.origins:
            for s in r.getDests():
                r.demand[s] = q[(r,s)]
            
        self.network.tapas("UE", None)
        
        for r in self.network.origins:
            self.network.dijkstras(r, "UE")
            for s in r.destSet:
                output[(r,s)] = s.cost
                
        return output
        
    def updateQbounds(self, q_lb, q_ub, mu_lb, mu_ub):
        
        self.mu_lb = mu_lb
        self.mu_ub = mu_ub
        self.q_lb = q_lb
        self.q_ub = q_ub
        
        
        for r in self.network.origins:
            for s in r.destSet:
                self.rmp.q[(r,s)].lb = q_lb[(r,s)]
                self.rmp.q[(r,s)].ub = q_ub[(r,s)]
                
                mu_ub = self.mu_ub[(r,s)]
                mu_lb = self.mu_lb[(r,s)]

                #print("\t", "mccormick", (r,s), self.rmp.q[(r,s)].lb, self.rmp.q[(r,s)].ub, self.mu_lb[(r,s)], self.mu_ub[(r,s)])
                self.vfcut_q_ub[(r,s)].rhs = self.rmp.q[(r,s)].ub * self.rmp.pi[(r,s)] + self.rmp.q[(r,s)] * mu_lb - self.rmp.q[(r,s)].ub * mu_lb
                self.vfcut_pi_ub[(r,s)].rhs = self.rmp.q[(r,s)].lb * self.rmp.pi[(r,s)] + self.rmp.q[(r,s)] * mu_ub - self.rmp.q[(r,s)].lb * mu_ub
        
        
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
        

        for r in self.network.origins:
            for s in r.destSet:
                if (r,s) in self.q_target:
                    self.q_lb[(r,s)] = self.q_target[(r,s)]*0.8
                    self.q_ub[(r,s)] = self.q_target[(r,s)]*1.2
        
        self.rmp.parameters.read.scale = -1

        
        
        self.mu_lb = self.calcMuBounds(self.q_lb)
        self.mu_ub = self.calcMuBounds(self.q_ub)
                

        
        
        self.rmp.mu_a = {a:self.rmp.continuous_var(lb=0) for a in self.network.links}
        self.rmp.mu_w = {(r,s):self.rmp.continuous_var(lb=0) for r in self.network.origins for s in r.getDests()}

        self.rmp.x = {a:self.rmp.continuous_var(lb=0, name="x_"+str(a)) for a in self.network.links}
        
        self.rmp.q = {(r,s): self.rmp.continuous_var(lb=self.q_lb[(r,s)], ub=self.q_ub[(r,s)]) for r in self.network.origins for s in r.getDests()}
        
        self.rmp.xc = {(a,r):self.rmp.continuous_var(lb=0, name="xc_"+str(a)+"_"+str(r)) for a in self.network.links for r in self.network.origins}
        
        self.rmp.eta_oa = {a: self.rmp.continuous_var(lb = 0) for a in self.network.links}
        
        self.rmp.beta = {a:self.rmp.continuous_var(lb=0) for a in self.network.links}
        
        self.rmp.vf_ub = {(r,s): self.rmp.continuous_var(lb = 0) for r in self.network.origins for s in r.destSet}
        
        self.rmp.eta = {a : self.rmp.continuous_var(lb = 0) for a in self.network.links}
        self.rmp.pi = {(r,i): self.rmp.continuous_var(lb = 0) for r in self.network.origins for i in self.network.nodes}
        
        for r in self.network.origins:
            self.rmp.pi[(r,r)].ub = 0
        
        print("CG", self.useCG)
        
        
        
        # eta from dual
        
        for r in self.network.origins:
            self.network.dijkstras(r, "UE")
            for a in self.network.links:
                if not(a.end.pred is None and a.start.pred is None):
                    self.rmp.add_constraint(self.rmp.eta[a] >= self.rmp.pi[(r,a.end)] - self.rmp.pi[(r, a.start)] - a.t_ff)
        
        
        # OA for eta
        
        for a in self.network.links:
            intervals = 10
            for pct in range (0, intervals):
                x_pct = a.C * pct / (intervals/2) # 0 to 200% of capacity
                #print(a, pct, x_pct, a.C)
                eta_pct = a.getTravelTime(x_pct, "UE") - a.t_ff
                
                self.addOACut_eta(a, eta_pct)
                
        
        
        # OFV for ll
        self.rmp.add_constraint(sum(self.rmp.beta[a] for a in self.network.links) + sum(self.rmp.eta_oa[a] for a in self.network.links) <= (1+self.params.ub_eps) *sum(self.rmp.vf_ub[(r,s)] for r in self.network.origins for s in r.destSet))
        
        
        for r in self.network.origins:
            for s in r.destSet:
                #bounds are set in updateqbounds()
                self.vfcut_q_ub[(r,s)] = self.rmp.add_constraint(self.rmp.vf_ub[(r,s)] <= 1e15)
                self.vfcut_pi_ub[(r,s)] = self.rmp.add_constraint(self.rmp.vf_ub[(r,s)] <= 1e15)
        
        
        # avoid the x=0 solution
        for a in self.network.links:
            if a in self.x_target and self.addOApoint_x(a, 0) >= 0:
                self.addObjCut_x(a, 0, 0)
        
        for r in self.network.origins:
            for s in r.getDests():
                if (r,s) in self.q_target and self.addOApoint_q(r, s, 0) >= 0:
                    self.addObjCut_q(r, s, 0, 0)
        
        
        for a in self.network.links:
            self.rmp.add_constraint(self.rmp.beta[a] >= a.t_ff * self.rmp.x[a], ctname="beta_"+str(a)+"0")
        
       
        for a in self.network.links:
            self.rmp.add_constraint(sum(self.rmp.xc[(a,r)] for r in self.network.origins) == self.rmp.x[a], ctname="xc-eq_"+str(a)+"_r")


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

        
        
        self.rmp.objfunc_x = (1-self.obj_weight) * sum(self.rmp.mu_a[a] for a in self.network.links)
        self.rmp.objfunc_q = self.obj_weight * sum(self.rmp.mu_w[(r,s)] for r in self.network.origins for s in r.getDests())
        self.rmp.objfunc = self.rmp.objfunc_x + self.rmp.objfunc_q
        self.rmp.minimize(self.rmp.objfunc)
        
    def solveRMP(self):
    
        t_solve = time.time()
        self.rmp.solve(log_output=False)
        t_solve = time.time() - t_solve
        
        if self.rmp.solve_details.status == 'infeasible' or self.rmp.solve_details.status == 'integer infeasible':
            return 'infeasible', dict(), dict(), 1e15
        
        if self.rmp.solve_details.status == 'optimal with unscaled infeasibilities':
            print(self.rmp.solve.details.status)
            return 'infeasible', dict(), dict(), 1e15
            
            
        

        
        x_l = {a:self.rmp.x[a].solution_value for a in self.network.links}
        q_l = {(r,s): self.rmp.q[(r,s)].solution_value for r in self.network.origins for s in r.getDests()}
        obj_l = self.rmp.objective_value
        
        '''
        print("pi")
        for r in self.network.origins:
            for s in r.destSet:
                print((r,s), self.rmp.pi[(r,s)].solution_value, self.mu_ub[(r,s)])
        
        print("eta")
        for a in self.network.links:
            g = a.alpha / pow(a.C, a.beta)
            p = a.beta
            print(a, self.rmp.eta[a].solution_value, p / (g * (p+1)) * self.rmp.eta[a].solution_value, self.rmp.eta_oa[a].solution_value)
         '''   
            
        return "solved", x_l, q_l, obj_l
        
    def getEtaGap(self):
        total = 0
        
        for a in self.network.links:
            g = a.alpha / pow(a.C, a.beta)
            p = a.beta
            #print(a, self.rmp.eta[a].solution_value, p / (g * (p+1)) * self.rmp.eta[a].solution_value , self.rmp.eta_oa[a].solution_value)
            total +=  max(0, p / (g * (p+1)) * self.rmp.eta[a].solution_value - self.rmp.eta_oa[a].solution_value)
            
        return total
        
    def TAP(self, q):
    
        for r in self.network.origins:
            for s in r.getDests():
                r.demand[s] = q[(r,s)]
                
                
            
        self.network.tapas("UE", None)
        xhat = {a:a.x for a in self.network.links}
        obj_f = self.calcOFV(xhat, q)
        
        return xhat, obj_f
        
    def printSolution(self, x, q):
        print("link flows")
        total = 0
        
        for a in self.network.links:
            if a in self.x_target:
                val = (1-self.obj_weight) * (x[a] - self.x_target[a]) **2
                print("\t", a, round(x[a], 2), round(self.x_target[a], 2), round(val, 2))
                total += val
            else:
                print("\t", a, round(x[a], 2))
           
        print("demand")
        for r in self.network.origins:
            for s in r.getDests():
                val = self.obj_weight * (q[(r,s)] - self.q_target[(r,s)]) **2
                print("\t", r, s, round(q[(r,s)], 2), round(self.q_target[(r,s)], 2), round(val, 2))
                total += val
        
        print(total)
        
        
    def getPaths(self):
        all_paths = []
        for r in self.network.origins:
            for s in self.paths[r].keys():
                for p in self.paths[r][s]:
                    all_paths.append(p)
        return all_paths  
        
    
    def solveFixed(self, xfix, q_lb, q_ub):
        
        
        for a in self.network.links:
            self.fix_link_cuts[a] = self.rmp_fix.add_constraint(self.rmp_fix.x[a] == xfix[a])
        
        for r in self.network.origins:
            for s in r.getDests():
                self.fix_q_bounds_lb = self.rmp_fix.add_constraint(self.rmp_fix.q[(r,s)] >= q_lb[(r,s)])
                self.fix_q_bounds_ub = self.rmp_fix.add_constraint(self.rmp_fix.q[(r,s)] <= q_ub[(r,s)])

                
        #self.addTAPDualFixed(xfix)
        
        
        self.rmp_fix.solve(log_output=False)
        
        print(self.rmp.solve_details.status)
            
            
        q_l = {(r,s): self.rmp.q[(r,s)].solution_value for r in self.network.origins for s in r.getDests()}
        obj_l = self.rmp.objective_value    
        
        for a in self.network.links:
            self.rmp_fix.remove_constraint(self.fix_link_cuts[a])
            
        for r in self.network.origins:
            for s in r.getDests():
                self.rmp_fix.remove_constraint(self.fix_q_bounds_lb)
                self.rmp_fix.remove_constraint(self.fix_q_bounds_ub)

        
        #self.rmp_fix.remove_constraint(self.rmp_fix.tap_dual)
        
        return obj_l, q_l
        

    def initRMPFix(self):   
        self.rmp_fix = Model()
        

        self.rmp_fix.parameters.read.scale = -1

        self.fix_link_cuts = dict()
        self.fix_q_bounds_lb = dict()
        self.fix_q_bounds_ub = dict()
        
        
        
        self.rmp_fix.mu_w = {(r,s):self.rmp_fix.continuous_var(lb=0) for r in self.network.origins for s in r.getDests()}

        self.rmp_fix.x = {a:self.rmp_fix.continuous_var(lb=0, name="x_"+str(a)) for a in self.network.links}
        
        self.rmp_fix.q = {(r,s): self.rmp_fix.continuous_var() for r in self.network.origins for s in r.getDests()}
        
        self.rmp_fix.xc = {(a,r):self.rmp_fix.continuous_var(lb=0, name="xc_"+str(a)+"_"+str(r)) for a in self.network.links for r in self.network.origins}


        for a in self.network.links:
            self.rmp_fix.add_constraint(sum(self.rmp_fix.xc[(a,r)] for r in self.network.origins) == self.rmp_fix.x[a], ctname="xc-eq_"+str(a)+"_r")


            for i in self.network.nodes:                    
                for r in self.network.origins:            

                    if i.id == r.id:
                        dem = - sum(self.rmp_fix.q[(r,s)] for s in r.getDests())                
                    elif isinstance(i, type(r)) == True and i in r.getDests():
                        dem = self.rmp_fix.q[(r,i)]
                    else:
                        dem = 0

                    #print(i, r, dem)

                    self.rmp_fix.add_constraint(sum(self.rmp_fix.xc[(a,r)] for a in i.incoming) - sum(self.rmp_fix.xc[(a,r)] for a in i.outgoing) == dem, ctname="cons_"+str(i)+"_"+str(r))


        self.rmp_fix.objfunc_q = self.obj_weight * sum(self.rmp_fix.mu_w[(r,s)] for r in self.network.origins for s in r.getDests())
        self.rmp_fix.objfunc = self.rmp_fix.objfunc_q
        self.rmp_fix.minimize(self.rmp_fix.objfunc)
        
    def addTAPDualFixed(self, x):
        eta = dict()
        pi = dict()
        
        for a in self.network.links:
            eta[a] = 0
            a.x = x[a] # for travel time calculation to be correct
            
        for r in self.network.origins:
            
            self.network.dijkstras(r, "UE")
            
            for i in self.network.nodes:
                pi[(r,i)] = i.cost
        
        
        #self.rmp.add_constraint(self.rmp.eta[a] >= self.rmp.pi[(r,a.end)] - self.rmp.pi[(r, a.start)] - a.t_ff)
        for a in self.network.links:
            eta[a] = max(eta[a], pi[r, a.end] - pi[r,a.start] - a.t_ff)
            
            
        Beckmann = sum(a.getPrimitiveTravelTimeC(a.x, 0) for a in self.network.links)
        
        demand = sum(pi[(r,s)] * self.rmp_fix.q[(r,s)] for r in self.network.origins for s in r.getDests())
        
        eta_term = 0
        for a in self.network.links:
            g = a.getConst()
            p = a.beta
            ge = pow(g, 1/p)
            eta_term += p / ((p+1) * ge) * pow(eta[a], (p+1)/p)
            
        self.rmp_fix.tap_dual = self.rmp_fix.add_constraint(Beckmann <= demand - eta_term)
                