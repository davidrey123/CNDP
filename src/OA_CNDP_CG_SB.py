import time
from src import Params
from docplex.mp.model import Model
from src import BB_node

class OA_CNDP_CG_SB:
    
    def __init__(self, network, inflate_costs, useCG=True):
        self.network = network
        self.CG_tol = 1e-4
        self.inf = 1e+9
        
        self.useCG = useCG
        self.solveSO_only = False
        
        
        self.gap = 100
        self.tstt = 10000000000
        self.last_xhat = None
        self.last_obj_f = 1000000000
        
        self.sameycuts = []
        self.sameyvars = []
        
        self.g = {a:a.cost * inflate_costs for a in self.network.links}
        
        for a in self.network.links2:
            a.y = 0
            a.add_cap = 0
            
        self.varlinks = []
        
        for a in self.network.links:
            #print(a, a.cost)
            if a.cost > 1e-6:
                
                self.varlinks.append(a)
                
        print("yvars", len(self.varlinks))
        
        self.paths = {r:{s:[] for s in self.network.zones} for r in self.network.origins}
        self.link_cons = dict()
        self.dem_cons = dict()
        
        self.rt_pricing = 0.0
        self.params = Params.Params()
        
        
        self.xl_points = dict()
        self.yl_points = dict()
        self.xf_points = list()
        self.vf_rhs_cut = list()
        
        for a in self.network.links:
            self.xl_points[a] = []
            
        for a in self.varlinks:
            self.yl_points[a] = []
            
       
        
        
    
    def solve(self):
       
        self.timelimit = 3600
        iteration = 0
        starttime = time.time()
        ub = 1e100
        lb = 0
        
        last_lb = 0
        
        gap = 1
        cutoff = self.params.min_CNDP_gap
        self.tap_time = 0
        
        elapsed = 0
        
        
        self.best_y = {a:-1 for a in self.network.links}
        self.best_x = {a: -1 for a in self.network.links}
        self.best_ub = 1e15
        
        
        self.add_vf_cut = False
                
        y_l = None
        x_f = None
        x_l = None
        
        first_node_iter = 20
        other_node_iter = 1
        node_iter = first_node_iter
        
        
        max_iter = 500
        
        B_f = 10000000
        obj_f = 100000000
        
        self.initRMP()
        
        bbnodes = []
        
        backtrack_nodes = []
        
        y_lb = {a:0 for a in self.network.links}
        y_ub = {a:0 for a in self.network.links}
        
        for a in self.varlinks:
            y_ub[a] = a.max_add_cap
        
        bbnodes.append(BB_node.BB_node(y_lb, y_ub, 0, 1e15))
        
        while iteration < max_iter and len(bbnodes) > 0 and gap > cutoff:
            iteration += 1
            
            # solve RMP -> y, LB
            
            lowest_lb = 1e15
            largest_diff = 0
            best_idx = -1
            idx = 0
            best_node = None
            
            next_bbnodes = []
            
            for n in bbnodes:
                if n.lb < self.best_ub:
                    next_bbnodes.append(n)
                    idx += 1
                    
                    if n.lb < lowest_lb or (n.lb == lowest_lb and n.largest_diff > largest_diff):
                        largest_diff = n.largest_diff
                        lowest_lb = n.lb
                        best_idx = idx-1
                        best_node = n
                    
            bbnodes = next_bbnodes
            
            
            if self.params.PRINT_BB_NODE:
                print("avail nodes", len(bbnodes))
                for n in bbnodes:
                    self.printNode(n)
            
            del bbnodes[best_idx]
            
            if self.params.PRINT_BB_INFO:
                print("evaluate node")
                self.printNode(best_node)
            
            status, local_lb, local_ub = self.solveNode(best_node, node_iter, starttime)
            
            best_node.lb = local_lb
            
            if best_node.lb < self.best_ub:
                backtrack_nodes.append(best_node)
            
            min_lb = local_lb
                
            for n in bbnodes:
                min_lb = min(min_lb, n.lb)
                
            gap = (self.best_ub - min_lb)/min_lb
            
            elapsed = time.time() - starttime
            
            print(iteration, min_lb, self.best_ub, local_ub, gap, elapsed)
            
            if iteration > 1 and gap < cutoff:
                print("ending due to low gap", iteration)
                break
                
            if elapsed > self.timelimit:
                break
            
            if status == 'solved' and lb < self.best_ub:
                # branch
                y_lb_1 = best_node.y_lb.copy()
                y_lb_2 = best_node.y_lb.copy()
                y_ub_1 = best_node.y_ub.copy()
                y_ub_2 = best_node.y_ub.copy()
                
                
                worst = self.findWorstGap()
                
                # this could happen if y=y_lb or y=y_ub
                if worst is None:
                    continue
                    
                split = self.findBranchSplit(worst)
                
                
                
                if split is not None:
                
                    if self.params.PRINT_BB_INFO:
                        print("branch on ", worst, self.rmp.y[worst].lb, self.rmp.y[worst].ub, split)


                    y_ub_1[worst] = split
                    y_lb_2[worst] = split

                    left = BB_node.BB_node(y_lb_1, y_ub_1, local_lb, local_ub)
                    right = BB_node.BB_node(y_lb_2, y_ub_2, local_lb, local_ub)

                    bbnodes.append(left)
                    bbnodes.append(right)

            node_iter = other_node_iter
            
            
 
        
        self.rmp.end()
        self.tstt = self.best_ub
        self.gap = gap
        
        elapsed = time.time() - starttime
            
        print(self.best_y)
        return ub, elapsed, self.tap_time, iteration
    
    def printNode(self, n):
        print("\t", n.lb, n.ub)
        
        if self.params.PRINT_BB_NODE:
            for a in self.varlinks:
                print("\t\t", a, n.y_lb[a], n.y_ub[a])

            
    def findBranchSplit(self, a):
    
        possible = []
 
        y_ub = self.rmp.y[a].ub
        y_lb = self.rmp.y[a].lb
        
        for xf_vec in self.xf_points:
            xf = xf_vec[a]
            
            rho_one = (a.getPrimitiveTravelTimeC(xf, y_ub) - a.getPrimitiveTravelTimeC(xf, y_lb)) / ((y_ub - y_lb) * a.t_ff * a.alpha * (-a.beta-1) )
            
            rho_two = xf / pow(rho_one * (a.beta+1), 1.0/(a.beta+1))
            theta = (rho_two - a.C - y_lb) / (y_ub-y_lb)
            y = theta * y_ub + (1-theta) * y_lb
            
            if y >= y_lb and y <= y_ub:
                possible.append(y)
            else:
                possible.append( (y_lb+y_ub)/2)
                #print(a, y_lb, y, y_ub)
            
            
           
            
            
        best_gap = 0
        best_y = None
        
  
        
        for y in possible:
            gap = self.calcVFgap(a, y)
            
            
            if gap > best_gap:
                best_gap = gap
                best_y = y
                
 
        return best_y
            
            
    
    def findWorstGap(self):
        worst = None
        worst_gap = 0

        total_gap = 0

        for a in self.varlinks:
        
            gap = 0
            
            yl = self.rmp.y[a].solution_value
            if not (yl == self.rmp.y[a].lb or yl == self.rmp.y[a].ub):
                gap = self.calcVFgap(a, self.rmp.y[a].solution_value)

            total_gap += gap
            

            if gap > worst_gap:
                worst_gap = gap
                worst = a
                
        '''        
        if worst is None:
            for a in self.varlinks:
                gap = self.calcVFgap(a, self.rmp.y[a].solution_value)
                print(a, gap, self.rmp.y[a].lb, self.rmp.y[a].ub, self.rmp.y[a].solution_value)
        '''
        
                
        return worst
        
        
    def solveNode(self, bbnode, max_iter, starttime):
    
        last_x_f = {a:-1 for a in self.network.links}
        last_x_l = {a:-1 for a in self.network.links}
        last_y_l = {a:-1 for a in self.network.links}
        last_lb = 1e15
        local_ub = 1e15
        
        self.updateYbounds(bbnode.y_lb, bbnode.y_ub)
        
        cutoff = self.params.min_CNDP_gap
        
        lb = 1e15
        iter = 0
        
        run_TAP = True
        add_cut = False
        
        if self.params.PRINT_BB_INFO:
            print("solve node")
        
        
        while iter < max_iter:
            iter += 1
            
            if self.useCG:
                SP_status, obj_l, x_l, y_l = self.CG()
            else:
                SP_status, obj_l, x_l, y_l = self.solveRMP()
            
            if SP_status == "infeasible":
                elapsed = time.time() - starttime
                if self.params.PRINT_BB_INFO:
                    print("end node due to infeasible")
                
                return "infeasible", 1e15, 1e15
            
            
            elapsed = time.time() - starttime
            
            if elapsed >= self.params.BB_timelimit:
                if self.params.PRINT_BB_INFO:
                    print("end due to elapsed time")
                break
            
            lb = obj_l
            B_l = self.calcBeckmann(x_l, y_l)
            # solve TAP -> x, UB

            run_TAP, add_cut = self.isYDifferent(y_l, last_y_l)

            if run_TAP:
                if self.params.PRINT_BB_INFO:
                    print("\t\tSolving TAP")
                t1 = time.time()



                x_f, obj_f = self.TAP(y_l)


                t1 = time.time()-t1
                self.tap_time += t1
                B_f = self.calcBeckmann(x_f, y_l)


                if local_ub > obj_f:
                    local_ub = obj_f
                    
                if self.best_ub > obj_f:
                    self.best_ub = obj_f
                    self.best_y = y_l
                    self.best_x = x_f

                if self.add_vf_cut:
                    self.addVF_RHScut(x_f, y_l, bbnode.y_lb, bbnode.y_ub)
                else:
                    self.add_vf_cut = True
                    
                    
            else:
                if self.params.PRINT_BB_INFO:
                    print("\tSkipping TAP")
                x_f = last_x_f

                # y may be similar enough to skip TAP but different enough to add new cut
              
              

            
            
            
                
                
            #self.checkVFCut(x_l, yhat, x_l, xhat, yhat)
            # add TSTT cut
            
            #self.addTSTTCut(x_l, y_l)
            #self.addBeckmannOACut(x_l, y_l)
            
            self.addXlCuts(x_l, y_l)
            
            if self.params.PRINT_BB_INFO:
                print("num cuts", sum(len(self.xl_points[a]) for a in self.network.links), "OA cuts", len(self.xf_points), "VF cuts")

            
            elapsed = time.time() - starttime
            
            if elapsed >= self.params.BB_timelimit:
                if self.params.PRINT_BB_INFO:
                    print("end due to elapsed time")
                break
                

            if lb > 0:
                gap = (self.best_ub - lb)/lb
            else:
                gap = 1
            
            if self.params.PRINT_BB_BASIC:
                print("\t", iter, lb, self.best_ub, obj_f, gap, elapsed, self.tap_time)
            
            if iter > 1 and gap < cutoff:
                print("end node due to low gap", iter)
                break
                
                
            #for a in self.varlinks:
                #print(a, self.rmp.theta[a].solution_value, self.rmp.theta[a].solution_value * self.rmp.y[a].ub + (1-self.rmp.theta[a].solution_value)*self.rmp.y[a].lb, self.rmp.y[a].solution_value)
            
            
            if self.params.PRINT_BB_INFO:
                B_approx = sum(self.rmp.beta[a].solution_value for a in self.network.links)
                print("\t\tBeckmann", B_l, B_f, B_l-B_f, B_approx, B_approx-B_f)
                
            #for a in self.network.links:
                #print("\t", a, xhat[a], x_l[a], yhat[a])
                #print("\t\t", x_l[a] * a.getTravelTimeC(x_l[a], yhat[a], "UE"), self.rmp.mu[a].solution_value)
            
            
            if iter > 1 and (lb - last_lb) / lb < 1e-3:
                if self.params.PRINT_BB_INFO:
                    print("end solve node due to small change")
                break
 

                
            
            

                
            if not run_TAP and last_lb == lb:
                break
                
            last_xhat = x_f
            last_yhat = y_l
            last_x_l = x_l
            last_lb = lb
            
        
        return "solved", lb, local_ub
    
    def calcVFgap(self, a, y):
        y_lb = self.rmp.y[a].lb
        y_ub = self.rmp.y[a].ub
        
        if y_lb == y_ub:
            return 0
        
        theta = (y - y_lb) / (y_ub - y_lb)
        
        output = 1e15
        
        for x_f in self.xf_points:
            actual = a.getPrimitiveTravelTimeC(x_f[a], y)
            approx = theta * a.getPrimitiveTravelTimeC(x_f[a], y_ub) + (1-theta)* a.getPrimitiveTravelTimeC(x_f[a], y_lb)
            diff = approx - actual
            
            #print("\t", actual, approx, theta, y, y_lb, y_ub)
            #print("\t\tgap", r, s, actual, approx, theta, self.rmp.y[(r,s)].solution_value, self.rmp.y[(r,s)].lb, self.rmp.y[(r,s)].ub)
            output = min(output, diff)
            
        return output        
     
   
    
    def addVF_RHScut(self, x_f, y_l, y_lb, y_ub):


        if self.isDifferentX(x_f, self.xf_points):
            self.vf_rhs_cut.append(self.rmp.add_constraint(sum(self.rmp.beta[a] for a in self.network.links) <= sum(self.calcVF_RHScut_val(a, x_f[a], y_lb[a], y_ub[a]) for a in self.network.links)))

            self.xf_points.append(x_f)
                
    def updateYbounds(self, y_lb, y_ub):      
        for a in self.varlinks:
            self.rmp.y[a].lb = y_lb[a]
            self.rmp.y[a].ub = y_ub[a]
            
            
            self.rmp.remove_constraint(self.theta_constraint[a])
            self.theta_constraint[a] = self.rmp.add_constraint(self.rmp.y[a] == self.rmp.theta[a] * y_ub[a] + (1-self.rmp.theta[a]) * y_lb[a])
        
        for n in range(0, len(self.vf_rhs_cut)):
            self.rmp.remove_constraint(self.vf_rhs_cut[n])
            
            self.vf_rhs_cut[n] = self.rmp.add_constraint(sum(self.rmp.beta[a] for a in self.network.links) <= sum(self.calcVF_RHScut_val(a, self.xf_points[n][a], y_lb[a], y_ub[a]) for a in self.network.links))
        
           
                
    def isDifferentX(self, x, x_list):
    

        for xhat in x_list:
            
            for a in self.network.links:
                #print("compare", xhat[a], x[a], a.C, abs(xhat[a] - x[a]) / a.C )
                if abs(xhat[a] - x[a]) / a.C > self.params.VF_TOL_X:
                    return True
        return len(x_list) == 0
        
    def calcVF_RHScut_val(self, a, x_f, y_lb, y_ub):
        output = 0
        
        if a in self.varlinks:
            output += self.rmp.theta[a] * a.getPrimitiveTravelTimeC(x_f, y_ub) + (1-self.rmp.theta[a])* a.getPrimitiveTravelTimeC(x_f, y_lb)
        else:
            output += a.getPrimitiveTravelTimeC(x_f, 0)

        return output

    def calcTSTT(self, x, y):
        output = 0
        for a in self.network.links:
            y_ext = 0
            
            if a in self.varlinks:
                y_ext = y[a]
            
            output += x[a] * a.getTravelTimeC(x[a], y_ext, "UE")
            
            
        return output
    
    def getAvgLinkCost(self):
        return sum(self.g[a] for a in self.varlinks) / len(self.varlinks)
        
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
        
    

    def calcBeckmann(self, xhat, yhat):
        total = 0
        
        for a in self.network.links:
            if a in self.varlinks:
                total += a.getPrimitiveTravelTimeC(xhat[a], yhat[a])
            else:
                total += a.getPrimitiveTravelTimeC(xhat[a], 0)
                
        return total
            
    
        
        
    def addBeckmannOACut(self, xl, yl):
    
        #print("\tadding Beckmann OA cut")
        
        for a in self.network.links:
            self.addBeckmannOACutLink(a, xl[a], yl[a])
    
    def addBeckmannOACutLink(self, a, xl, yl):
        y_ext = 0
        yterm = 0

        if a in self.varlinks:
            y_ext = yl
            yterm = (self.rmp.y[a] - yl) * a.intdtdy(xl, yl) 

        xterm = (self.rmp.x[a] - xl) * a.getTravelTimeC(xl, y_ext, "UE")

        B = a.getPrimitiveTravelTimeC(xl, y_ext)

        self.rmp.add_constraint(self.rmp.beta[a] >= B + yterm + xterm)
        #print("\tBcut", xl[a], lastx[a], a.getTravelTimeC(xl[a], y_ext, "UE"))
        
    def addXlCuts(self, x_l, y_l):
    
        #print("\tadding Beckmann OA cut")
        
        for a in self.network.links:
            
            y_ext = 0
            
            if a in self.varlinks:
                y_ext = y_l[a]
                
            if self.addLeadPoint(a, x_l[a], y_ext):
                self.addBeckmannOACutLink(a, x_l[a], y_ext)
                self.addTSTTCutLink(a, x_l[a], y_ext)
        
        
    def checkVFCut(self, x, y, x_l, xhat, yhat):
        for a in self.network.links:
            if a in self.varlinks:
                print("\tbeta", a, self.rmp.beta[a].solution_value, a.getPrimitiveTravelTimeC(x[a], y[a]), self.rmp.mu[a].solution_value)
            else:
                print("\tbeta", a, self.rmp.beta[a].solution_value, a.getPrimitiveTravelTimeC(x[a], 0), self.rmp.mu[a].solution_value)
    
    def addTSTTCut(self, xhat, yhat):
        for a in self.network.links:
            self.addTSTTCutLink(a, xhat[a], yhat[a])
                
    def addTSTTCutLink(self, a, xhat, yhat):
        if a in self.varlinks:
            firstterm = xhat * a.getTravelTimeC(xhat, yhat, "UE")
            secondterm = a.getTravelTimeC(xhat, yhat, "UE") + xhat * a.getDerivativeTravelTimeCx(xhat, yhat)
            thirdterm = xhat * a.getDerivativeTravelTimeCy(xhat, yhat)

            self.rmp.add_constraint(self.rmp.mu[a] >= firstterm + secondterm * (self.rmp.x[a] - xhat) + thirdterm * (self.rmp.y[a] - yhat))
        else:
            firstterm = xhat * a.getTravelTimeC(xhat, 0, "UE")
            secondterm = a.getTravelTimeC(xhat, 0, "UE") + xhat * a.getDerivativeTravelTimeCx(xhat, 0)

            self.rmp.add_constraint(self.rmp.mu[a] >= firstterm + secondterm * (self.rmp.x[a] - xhat))
        
    
      
    def initRMP(self):   
    
        # init paths
        
        if self.useCG:
            for r in self.network.origins:
                if r.totaldemand > 0:
                    self.network.dijkstras(r, "UE")

                    self.paths[r] = dict()

                    for s in self.network.zones:
                        if r.getDemand(s) > 0:
                            p = self.network.trace(r, s)
                            self.paths[r][s] = list()
                            self.paths[r][s].append(p)
                    
                
        self.rmp = Model()
        
        
        self.rmp.parameters.read.scale = -1
        
        
        self.rmp.mu = {a:self.rmp.continuous_var(lb=0,ub=1e10) for a in self.network.links}
        self.rmp.beta = {a:self.rmp.continuous_var(lb=0) for a in self.network.links}
        
        self.rmp.y = {a:self.rmp.continuous_var(lb=0, ub=a.max_add_cap) for a in self.varlinks}
        
        
        
        self.rmp.theta = {a:self.rmp.continuous_var(lb=0, ub=1) for a in self.varlinks}
        self.theta_constraint = dict()
        
        for a in self.varlinks:
            self.theta_constraint[a] = self.rmp.add_constraint(self.rmp.y[a] == self.rmp.theta[a] * a.max_add_cap)
        
        self.rmp.x = {a:self.rmp.continuous_var(lb=0, ub=self.network.TD) for a in self.network.links}
        
        if self.useCG:
            self.rmp.h = {p:self.rmp.continuous_var(lb=0) for p in self.getPaths()}
            
            for a in self.network.links:
                self.link_cons[a] = self.rmp.add_constraint(self.rmp.x[a] - sum(self.rmp.h[p] for p in self.getPaths() if a in p.links) >= 0)

            for r in self.network.origins:
                for s in self.network.zones:
                    if r.getDemand(s) > 0:
                        self.dem_cons[(r,s)] = self.rmp.add_constraint(sum(self.rmp.h[p] for p in self.paths[r][s]) >= r.getDemand(s))
        else:
            self.rmp.xc = {(a,r):self.rmp.continuous_var(lb=0, ub=r.totaldemand) for a in self.network.links for r in self.network.origins}
            
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
        
        for a in self.network.links:
            self.rmp.add_constraint(self.rmp.mu[a] >= self.rmp.x[a] * a.t_ff)
        
        
        
        
        
        self.rmp.minimize(sum(self.rmp.mu[a] for a in self.network.links) + sum(self.g[a] * self.rmp.y[a] for a in self.varlinks))
    
    
    
    
    
    def pricing(self, link_duals, dem_duals):
        
        t0_pricing = time.time()
        
        new = 0
        minrc = self.inf
        
        #print("len check", len(link_duals))
        for a in self.network.links:
            a.dual = link_duals[a]
        
        for r in self.network.origins:
            self.network.dijkstras(r,'RC')
            
            for s in self.network.zones:
                
                if r.getDemand(s) > 0:
                
                    rc = - dem_duals[(r,s)] + s.cost                    
                    
                    if rc < - self.CG_tol:
                        p = self.network.trace(r,s)
                        self.paths[r][s].append(p) #---is it needed to store paths if directly adding to RMP?
                        
                        #---add new path var to RMP                                                
                        self.rmp.h[p] = self.rmp.continuous_var(lb=0)
                        
                        #---update RMP constraints
                        self.dem_cons[(r, s)].lhs.add_term(self.rmp.h[p], 1)
                        
                        for a in p.links:
                            self.link_cons[a].lhs.add_term(self.rmp.h[p], -1)
                        
                        new += 1
                    
                        if rc < minrc:
                            minrc = rc
                
        self.rt_pricing += (time.time() - t0_pricing)        
        return minrc
        
    def CG(self):
        conv = False
        nCG = 0
        OFV = self.inf
        
        starttime = time.time()
        
        
        while conv == False:
            RMP_status, OFV, link_duals, dem_duals = self.solveRMP()
            #print("solved?", RMP_status)
            
            
            if RMP_status == 'infeasible':
                CG_status = 'infeasible'
                
                return CG_status, self.inf, dict(), dict()
            else:
                CG_status = "solved"
                
            minrc = self.pricing(link_duals, dem_duals)

            elapsed = time.time() - starttime
            
            if self.params.PRINT_CG_INFO:
                npaths = len(self.getPaths())
                print('CG: %d\t%d\t%.1f\t%.2f\t%.2f' % (nCG,npaths,OFV,minrc, elapsed))
                
            if minrc > -self.CG_tol:
                conv = True
                
            nCG += 1
        
        yhat = {a:self.rmp.y[a].solution_value for a in self.varlinks}
        x_l = {a:self.rmp.x[a].solution_value for a in self.network.links}
 
        
        return CG_status, OFV, x_l, yhat
        
        
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
        
        
        if self.params.PRINT_CG_INFO:
            print("\t\t", RMP_status)
        
        OFV = self.rmp.objective_value
        
        if self.useCG:
            link_duals = {a: self.link_cons[a].dual_value for a in self.network.links}

            dem_duals = dict()
            for r in self.network.origins:
                for s in self.network.zones:
                    if r.getDemand(s) > 0:
                        dem_duals[(r,s)] = self.dem_cons[(r,s)].dual_value
        
        
        
            return RMP_status, OFV, link_duals, dem_duals
        else:
            yhat = {a:self.rmp.y[a].solution_value for a in self.varlinks}
            x_l = {a:self.rmp.x[a].solution_value for a in self.network.links}
            return RMP_status, OFV, x_l, yhat
        
    def isYDifferent(self, y, lasty):
        output1 = False
        for a in self.varlinks:
            if abs(y[a] - lasty[a]) > 1e-9:
                return True, True
            elif abs(y[a] - lasty[a]) > 1e-9:
                output1 = True
        return output1, False
        
    def TAP(self, y):
    
        for a in self.varlinks:
            a.add_cap = y[a]
            
        self.network.tapas("UE", None)
        xhat = {a:a.x for a in self.network.links}
        obj_f = self.calcOFV(xhat, y)

        return xhat, obj_f
        
        
    def getPaths(self):
        all_paths = []
        for r in self.network.origins:
            for s in self.paths[r].keys():
                for p in self.paths[r][s]:
                    all_paths.append(p)
        return all_paths   
        
   
    def addLeadPoint(self, a, xl, yl):
        for idx in range(0, len(self.xl_points[a])):
            px = self.xl_points[a][idx]
            
            if a in self.varlinks:
                py = self.yl_points[a][idx]

                if abs(px-xl) / a.C < self.network.params.OA_TOL_X and abs(py-yl) / a.C < self.network.params.OA_TOL_Y: # was tol
                    return False
            else:
                if abs(px-xl) / a.C < self.network.params.OA_TOL_X:
                    return False
        
        self.xl_points[a].append(xl)
        
        if a in self.varlinks:
            self.yl_points[a].append(yl)
        return True
        
    def addOApoint_helper(self, list, pt, tol):
        # decide whether to add this point and return true or false
        similar = False
        
        # this is inefficient, I should change it!
        for p in list:
            if abs(pt - p) < tol: # was tol
                similar = True
                break
                
        if not similar:
            list.append(pt)
            return len(list)-1
        else:
            return -1  
 