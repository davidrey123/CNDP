import time
from src import Params
from docplex.mp.model import Model

class OA_CNDP_CG:
    
    def __init__(self, network, inflate_costs, useCG=True, useLinkVF=True):
        self.network = network
        self.CG_tol = 1e-4
        self.inf = 1e+9
        
        self.useCG = useCG
        self.useLinkVF = useLinkVF
        self.solveSO_only = False
        
        print("CG", useCG, "link VF", useLinkVF)
        
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
        best_x = {a: -1 for a in self.network.links}
        
        
        
                
                
        yhat = None
        xhat = None
        x_l = None
        
        run_TAP = True
        add_cut = False
        
        B_f = 10000000
        obj_f = 100000000
        
        self.initRMP()
        
        while gap > cutoff:
            iteration += 1
            
            # solve RMP -> y, LB

            if self.useCG:
                SP_status, obj_l, x_l, yhat = self.CG()
            else:
                SP_status, obj_l, x_l, yhat = self.solveRMP()
            
            if SP_status == "infeasible":
                elapsed = time.time() - starttime
                exit()
                break
            
            
            elapsed = time.time() - starttime
            
            if elapsed >= self.params.BB_timelimit:
                break
            
            lb = obj_l
            B_l = self.calcBeckmann(x_l, yhat)
            # solve TAP -> x, UB

            if not self.solveSO_only: 
                # add VF cut
                if run_TAP:
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
                        best_x = xhat

                    if iteration > 1:
                        if self.useLinkVF:

                            self.addVFCut2(x_l, xhat, yhat)
                        else:
                            self.addVFCut(x_l, xhat, yhat)
                else:
                    if self.params.PRINT_BB_INFO:
                        print("\tSkipping TAP")
                    xhat = last_xhat
                    #self.addVFCutSameY(x_l, xhat, yhat)
                
                    # y may be similar enough to skip TAP but different enough to add new cut
                    if self.useLinkVF and add_cut and iteration > 1:
                    
                        self.addVFCut2(x_l, xhat, yhat)
                        
                if self.useLinkVF:
                    self.addBeckmannOACut(x_l, yhat)
            
                
                
            #self.checkVFCut(x_l, yhat, x_l, xhat, yhat)
            # add TSTT cut
                self.addTSTTCut(x_l, yhat)
            else:
                self.addTSTTCut(x_l, yhat)
            
            
            
            elapsed = time.time() - starttime
            
            if elapsed >= self.params.BB_timelimit:
                break
                
            if self.solveSO_only:
                gap = (lb - last_lb)/lb
            else:
                if lb > 0:
                    gap = (ub - lb)/lb
                else:
                    gap = 1
            
            print(iteration, lb, ub, obj_f, gap, elapsed, tap_time)
            
            
            
            if self.params.PRINT_BB_INFO:
                B_approx = sum(self.rmp.beta[a].solution_value for a in self.network.links)
                print("Beckmann", B_l, B_f, B_l-B_f, B_approx, B_approx-B_f)
                '''
                for a in self.network.links:
                    y_ext = 0
                    
                    if a in self.varlinks:
                        y_ext = yhat[a]
                    print("\t", a, a.getPrimitiveTravelTimeC(x_l[a], y_ext), self.rmp.beta[a].solution_value)
                print("obj", self.calcOFV(x_l, yhat), obj_l, self.calcOFV(xhat, yhat))
                '''
            #for a in self.network.links:
                #print("\t", a, xhat[a], x_l[a], yhat[a])
                #print("\t\t", x_l[a] * a.getTravelTimeC(x_l[a], yhat[a], "UE"), self.rmp.mu[a].solution_value)
            
            
 
            if elapsed > timelimit:
                break
                
            run_TAP, add_cut = self.isYDifferent(yhat, last_yhat)
            
            if gap < 0.2:
                self.network.params.min_gap = 1E-3
                run_TAP = True
                
            last_xhat = xhat
            last_yhat = yhat
            last_x_l = x_l
            last_lb = lb
            
        '''   
        for a in self.network.links:
            y_ext = 0
            
            if a in self.varlinks:
                y_ext = best_y[a]
            print("[", a.start, ",", a.end, ",", best_x[a], ",", y_ext, "], ")
        '''
        
        self.rmp.end()
        self.tstt = self.calcTSTT(best_x, best_y)
        self.gap = gap
            
        return ub, elapsed, tap_time, iteration

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
            
    def addVFCut(self, x_l, xhat, yhat):

        # doesn't work
        '''
        vba = {a:self.rmp.continuous_var(lb=-1e300, ub=1e300) for a in self.network.links}
        
        sum_b = 0
        for a in self.network.links:
            B1 = 0
            B2 = 0
            
            y_ext = 0
            
            if a in self.varlinks:
                y_ext = yhat[a]
            else:
                y_ext = 0
                
            B1 = a.getPrimitiveTravelTimeC(x_l[a], y_ext)
            B2 = a.getPrimitiveTravelTimeC(xhat[a], y_ext)

            #sum( (self.rmp.x[a] - x_l[a]) * a.getTravelTimeC(x_l[a], yhat_ext[a], "UE") for a in self.network.links) 
            xterm = (self.rmp.x[a] - x_l[a]) * a.getTravelTimeC(x_l[a], y_ext, "UE")
            
            yterm = 0
            
            if a in self.varlinks:
                yterm = (self.rmp.y[a] - yhat[a]) * (a.intdtdy(x_l[a], yhat[a]) - a.intdtdy(xhat[a], yhat[a]))
            #sum( (self.rmp.y[a] - yhat[a]) * (a.intdtdy(x_l[a], yhat[a]) - a.intdtdy(xhat[a], yhat[a])) for a in self.varlinks)
            
            #print(a, B1-B2)
            self.rmp.add_constraint(self.rmp.eta[a] >= B1-B2 + xterm + yterm)
        
            sum_b += B1-B2+xterm+yterm
        #self.rmp.add_constraint(sum_b <= 0)

        #self.rmp.add_constraint(sum(vba[a] for a in self.network.links) <= 0)
        '''    
        
        if self.params.PRINT_BB_INFO:
            print("adding VF cut Delta B")
           
        B1 = self.calcBeckmann(x_l, yhat)
        B2 = self.calcBeckmann(xhat, yhat)
        
        yhat_ext = dict()
        for a in self.network.links:
            if a in self.varlinks:
                yhat_ext[a] = yhat[a]
            else:
                yhat_ext[a] = 0
                
        self.rmp.add_constraint(B1-B2 + sum( (self.rmp.x[a] - x_l[a]) * a.getTravelTimeC(x_l[a], yhat_ext[a], "UE") for a in self.network.links) 
        + sum( (self.rmp.y[a] - yhat[a]) * (a.intdtdy(x_l[a], yhat[a]) - a.intdtdy(xhat[a], yhat[a])) for a in self.varlinks) <= 0)
    
    def addVFCut2(self, xl, xf, yl):
        lhs = sum(self.rmp.beta[a] for a in self.network.links) -1e-4
        
        lhs -= self.calcBeckmann(xf, yl)
        
        lhs += sum( (self.rmp.y[a] - yl[a]) * a.intdtdy(xf[a], yl[a]) for a in self.varlinks)
        
        self.rmp.add_constraint(lhs <= 0)
        
        #self.checkVFCut(xl, yl, xl, xf, yl)
        
        
    def addBeckmannOACut(self, xl, yl):
    
        #print("\tadding Beckmann OA cut")
        
        for a in self.network.links:
            y_ext = 0
            yterm = 0
            
            if a in self.varlinks:
                y_ext = yl[a]
                yterm = (self.rmp.y[a] - yl[a]) * a.intdtdy(xl[a], yl[a]) 
                
            xterm = (self.rmp.x[a] - xl[a]) * a.getTravelTimeC(xl[a], y_ext, "UE")
                
            B = a.getPrimitiveTravelTimeC(xl[a], y_ext)
            
            self.rmp.add_constraint(self.rmp.beta[a] >= B + yterm + xterm)
            #print("\tBcut", xl[a], lastx[a], a.getTravelTimeC(xl[a], y_ext, "UE"))
        
    def addVFCutSameY(self, xl, xf, yl):
    
        if self.params.PRINT_BB_INFO:
            print("adding VF cut same y")
       
            
        
        yl_ext = dict()
        for a in self.network.links:
            if a in self.varlinks:
                yl_ext[a] = yl[a]
            else:
                yl_ext[a] = 0

        ydiff = None
        yl_old = None
        
        for i in range(0, len(self.sameyvars)):
            if not self.isYDifferent(self.sameyvars[i], yl):
                ydiff = self.sameycuts[i]
                yl_old = self.sameyvars[i]
                break
            
        if ydiff == None:
            ydiff = {a:self.rmp.continuous_var(lb=0, ub=a.C) for a in self.varlinks}
            yl_old = yl
            
            for a in self.varlinks:
                self.rmp.add_constraint(ydiff[a] >= yl[a] - self.rmp.y[a])
                
            
            self.sameyvars.append(yl)
            self.sameycuts.append(ydiff)
            
            rhs = sum(a.getPrimitiveTravelTimeC(xf[a], yl_ext[a]) for a in self.network.links)  
            rhs += sum(-a.getDerivativeTravelTimeCy(xf[a], 0) * ydiff[a] for a in self.varlinks)
            
            #self.rmp.add_constraint(sum(self.rmp.beta[a] for a in self.network.links) <= rhs)
          
        #print("\tcheck1", sum(self.rmp.beta[a].solution_value for a in self.network.links), sum(a.getPrimitiveTravelTimeC(xf[a], yl_ext[a]) for a in self.network.links)  )

        #for a in self.varlinks:
        #    print("\t\tydiff", a, ydiff[a].solution_value, yl[a]-yl_old[a], a.getDerivativeTravelTimeCy(xf[a], 0))
        #self.checkVFCut(xl, yl, xl, xf, yl)

        
    def checkVFCut(self, x, y, x_l, xhat, yhat):
        for a in self.network.links:
            if a in self.varlinks:
                print("\tbeta", a, self.rmp.beta[a].solution_value, a.getPrimitiveTravelTimeC(x[a], y[a]), self.rmp.mu[a].solution_value)
            else:
                print("\tbeta", a, self.rmp.beta[a].solution_value, a.getPrimitiveTravelTimeC(x[a], 0), self.rmp.mu[a].solution_value)
    
    def addTSTTCut(self, xhat, yhat):
        for a in self.network.links:
            if a in self.varlinks:
                firstterm = xhat[a] * a.getTravelTimeC(xhat[a], yhat[a], "UE")
                secondterm = a.getTravelTimeC(xhat[a], yhat[a], "UE") + xhat[a] * a.getDerivativeTravelTimeCx(xhat[a], yhat[a])
                thirdterm = xhat[a] * a.getDerivativeTravelTimeCy(xhat[a], yhat[a])
            
                self.rmp.add_constraint(self.rmp.mu[a] >= firstterm + secondterm * (self.rmp.x[a] - xhat[a]) + thirdterm * (self.rmp.y[a] - yhat[a]))
            else:
                firstterm = xhat[a] * a.getTravelTimeC(xhat[a], 0, "UE")
                secondterm = a.getTravelTimeC(xhat[a], 0, "UE") + xhat[a] * a.getDerivativeTravelTimeCx(xhat[a], 0)

                self.rmp.add_constraint(self.rmp.mu[a] >= firstterm + secondterm * (self.rmp.x[a] - xhat[a]))
        
    
      
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
        
        while conv == False:
            RMP_status, OFV, link_duals, dem_duals = self.solveRMP()
            #print("solved?", RMP_status)
            
            
            if RMP_status == 'infeasible':
                CG_status = 'infeasible'
                
                return CG_status, self.inf, dict(), dict()
            else:
                CG_status = "solved"
                
            minrc = self.pricing(link_duals, dem_duals)

            if self.params.PRINT_BB_INFO:
                npaths = len(self.getPaths())
                print('CG: %d\t%d\t%.1f\t%.2f' % (nCG,npaths,OFV,minrc))
                
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
        
        
        if self.params.PRINT_BB_INFO:
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
        
    def TAP(self, y, lasty):
    
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
 