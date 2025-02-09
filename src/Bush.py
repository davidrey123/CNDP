# Created on : Mar 27, 2024, 5:02:24 PM
# Author     : michaellevin
import sys
#print(sys.executable)
import llist
import contextlib
from collections import deque
import heapq

from src import Params
from src import PASList
from src import Branch
from src import Zone


import llist
#print('llist module is available.')

from llist import sllist
from src import PAS
from src import NodeReturn
from src import Node

class Bush:
    def __init__(self, network, origin):
        self.network = network
        self.flow = {}
        self.relevantPAS = []
        self.sorted = list()
        self.branches = []
        self.origin = origin 
        
        self.demand = dict()
        
        for s in origin.destSet:
            self.demand[s] = 0
        
        #print(origin)
        self.relevantPAS = PASList.PASList()
        
        #for l in self.network.links:
        #    self.flow[l] = 0
            
        self.equilibrateDemand()
        #self.loadDemand()
        
    def contains(self, l):
        return l in self.flow and self.flow[l] > self.network.params.flow_epsilon
        
    def getFlow(self, l):
        return self.flow.get(l, 0)
        
    #def containsNode(self, n):
        #if n == self.origin:
            #print(n)
            #return True


        #for l in n.incoming:
            #if self.contains(l):
                #return True

        #return False
     
    def loadDemand(self):
        #with open('result130.txt', 'w') as file, contextlib.redirect_stdout(file):
            self.network.dijkstras(self.origin, self.network.type)
            #with open('result37.txt', 'w') as file, contextlib.redirect_stdout(file):
            for s in self.network.zones:
                d = self.origin.getDemand(s)
                #print(d)
                if d > 0:
                    curr = s
                    
                    while curr != self.origin:
                        #print(f"self origin is{self.origin}")
                        uv = curr.pred
                        #print(f"uv is {uv}")
                        self.addFlow(uv, d)
                        #print(self.addFlow(uv, d))
                        curr = uv.start
                        #print(f"curr is {curr}")
                        #print(curr)

                self.demand[s] = d
            self.topologicalSort()
        #print(self.topologicalSort)
        
        
    


    def topologicalSort(self):
        
        for n in self.network.nodes:
            n.in_degree = len(n.getBushIncoming(self))

        # Queue to store vertices with indegree 0
        q = deque()
        for n in self.network.nodes:
            if n.in_degree == 0:
                q.append(n)
        self.sorted = []
        
        idx = 0
        while q:
            i = q.popleft()
            self.sorted.append(i)
            i.top_order = idx
            idx += 1
            # Decrease indegree of adjacent vertices as the current node is in topological order
            for ij in i.getBushOutgoing(self):
                j = ij.end
                j.in_degree -= 1
                # If indegree becomes 0, push it to the queue
                if j.in_degree == 0:
                    q.append(j)
    
    def testTopologicalSort(self):
    
        self.topologicalSort()
        
        for l in self.flow:
            if self.contains(l) and l.start.top_order > l.end.top_order:
                return False

        return True
    
    
    def removeLinks(self, removedlinks):

        minPathTree = self.network.getSPTree(self.origin)
        
        for l in removedlinks:
            if self.contains(l):
                if l.end == self.origin:
                    continue
                elif l.start == self.origin:
                    # need a PAS shift here
                    pas = self.createOriginLinkPAS(minPathTree, l)
                    pas.zeroBackwardFlow(self)
                else:
                    branch = self.createBranch(l)
                    if branch is not None:
                        #print("branch shift ", self.origin, l, self.getFlow(l))
                        branch.init()
                        branch.flowShift(self.network.type)
                        #print(self.getFlow(l))
                    #else:
                        #print("branch is null ", self.origin, l)
            
    def bushSPTT(self):
        sptt = 0
        self.network.dijkstras(self.origin, self.network.type)
        
        for s in self.network.zones:
            if self.origin.getDemand(s) > 0:
                sptt += self.origin.getDemand(s) * s.cost
        return sptt
        
    def bushTSTT(self):
        tstt = 0
        
        for ij in self.network.links:
            if self.getFlow(ij) > 0:
                tstt += self.getFlow(ij) * ij.getTravelTime(ij.x, self.network.type)
        return tstt
        
    def checkEquilibrium(self):
        sptt = 0
        totflow = 0
        demand = 0
        self.network.dijkstras(self.origin, self.network.type)
        
        cost = {}
        
        for s in self.network.zones:
            if self.origin.getDemand(s) > 0:
                demand += self.origin.getDemand(s)
                sptt += self.origin.getDemand(s) * s.cost
            cost[s] = s.cost
        
        tstt = 0
        
        for ij in self.origin.outgoing:
            totflow += self.getFlow(ij)
        
        for ij in self.network.links:
            if self.getFlow(ij) > 0:
                tstt += self.getFlow(ij) * ij.getTravelTime(ij.x, self.network.type)
                
        self.minUsedTree()
        
        gap = (tstt-sptt)/tstt
        
        if gap > 1e-4:     
            print("check equilibrium ", self.origin, sptt, tstt, gap, demand - totflow, self.validateFlowConservation())
            return False
        return True
        
            
    def validateFlowConservation(self):

        for n in self.network.nodes:
            #print(n)
            inflow = 0
            outflow = 0
            
            for l in n.incoming:
                inflow += self.getFlow(l)
            

            for l in n.outgoing:
                outflow += self.getFlow(l)
            
            if n == self.origin:
                inflow += self.origin.getProductions()
            
            
            outflow += self.origin.getDemand(n)
            

            if abs(inflow-outflow) > self.network.params.flow_epsilon:
                return False

        return True
        
    def addLinks(self, newlinks):  
        # does this create a new shortest path?
        
        included = set()
        
        minPathTree = self.network.getSPTree(self.origin)
        
        for l in newlinks:
            if l not in included and l.end != self.origin and l.hasHighReducedCost(self.network.type, self.network.params.bush_gap):
                self.network.params.good_bush_gap += 1
                #print(f"included is {included}")
                #print(f"l.end is {l.end} and l.start is {l.start}")
                #print(f"self.flow of l is {self.flow}")
                #print(f"self.network.params.flow_epsilon is {self.network.params.flow_epsilon}")
                #print(f"the last thing {l.hasHighReducedCost(self.network.type, self.network.params.pas_cost_mu)}")
                #print(l)
                #print(f"selforigin is {self.origin}")
                #System.out.println(l+" "+l.getDest().cost+" "+l.getSource().cost+" "+(l.getDest().cost-l.getSource().cost)+" "+l.getTT());
                # we need a PAS!
                if not self.hasRelevantPAS(l):

                    # should check if we can borrow one from network
                    fromNetwork = self.network.findPAS(l, self)
                    if fromNetwork is None:

                        if self.network.params.PRINT_PAS_INFO:
                            print("\n\tCreate PAS for " + str(l) + " for origin " + str(self.origin))

                        newPAS = None
                        
                        if l.start != self.origin:
                            newPAS = self.createPAS(minPathTree, l, self.getFlow(l) * self.network.params.pas_flow_mu)
                        #else:
                        #    newPAS = self.createOriginLinkPAS(minPathTree, l)
                        #print(l)
                        
                        if newPAS != None:
                            self.savePAS(newPAS)
                            for ij in newPAS.backwardlinks:
                                included.add(ij)
                            newPAS.flowShift(self.network.type, self.network.params)
                            
                    else:
                        if self.network.params.PRINT_PAS_INFO:
                            print("Take PAS for " + str(l))

                        fromNetwork.addRelevantOrigin(self.origin)
                        #print(self.origin)
   
    
    #def checkReducedCosts(self):
    
        #self.minPath()
        
        #output = True
        #for l in self.network.links:
        
            #if self.flow[l] > self.network.params.flow_epsilon and l.getReducedCost() < -self.network.params.bush_gap:
                #print("Negative reduced cost origin "+str(self.origin.id))
                #print(str(l)+"\t"+str(flow[l])+"\t"+str(l.getSource().cost)+"\t"+str(l.getTT())+"\t"+str(l.getDest().cost) + "\t"+ str(l.getDest().cost - (l.getSource().cost + l.getTT())))
                #output = False

        #return output
    
    
    
    def algBShift(self):
    
        self.topologicalSort()
        
        minTree = self.minUsedTree()
        maxTree = self.maxUsedTree()
        
        #for n in self.network.nodes:
        #    print(n, n.cost, n.pred, self.origin.getDemand(n))
        
        #for ij in self.flow:
        #    print(ij, self.getFlow(ij))
            
        for s in self.network.zones:
            if self.origin.getDemand(s) < self.network.params.flow_epsilon:
                continue
            
            # find common node
            minNodes = set()
            curr = s
            
            #print("iter", s, self.origin, s.pred)
            while curr.pred != None:
                #print("\t", curr, s, self.origin)
                curr = curr.pred.start
                minNodes.add(curr)
                
            curr = s
            
            while curr.pred2 != None:
                #print("\tmax", curr, s, self.origin)
                if curr in minNodes:
                    break
                curr = curr.pred2.start
            
            common = curr
            
            # trace back minpath and maxpath
            minPath = self.tracePath(common, s)
            maxPath = self.tracePath2(common, s)
            
            #print("alg b paths", common, s, minPath, maxPath)
            
            # line search to equilibrate minpath and maxpath
            if minPath != None and maxPath != None:
                self.flowShift(minPath, maxPath)
    
    def getTT(self, add, path):
        output = 0
        for ij in path:
            output += ij.getTravelTime(ij.x + add, self.network.type)
        
        return output
        
    def flowShift(self, minPath, maxPath):
        maxshift = 1.0e9
        
        lowPath = minPath
        highPath = maxPath
        
        lowTT = self.getTT(0, lowPath)
        highTT = self.getTT(0, highPath)
        
        #print("starting alg b shift", lowTT, highTT, lowPath, highPath)
        
        if highTT - lowTT < self.network.params.bush_gap:
            return
        
        '''
        minPathTT = self.getTT(0, minPath)
        maxPathTT = self.getTT(0, maxPath)

        if minPathTT < maxPathTT:
            lowPath = minPath
            highPath = maxPath
        else:
            lowPath = maxPath
            highPath = minPath
        '''
        
        
        
        for ij in highPath:
            maxshift = min(maxshift, self.getFlow(ij))
            
        if maxshift < self.network.params.flow_epsilon:
            return
        
        bot = 0
        top = maxshift
        #with open('flowShift3.txt', 'a') as file, contextlib.redirect_stdout(file):
        while top - bot > self.network.params.line_search_gap:
            #print(line_search_gap)
            mid = (top + bot)/2
            
            check = self.getTT(mid, lowPath) - self.getTT(-mid, highPath)
            #print("\t", check, bot, top, mid)
            
            if check < 0:
                bot = mid
            
            else:
                top = mid
        
        for ij in lowPath:
            self.addFlow(ij, bot)
        
        for ij in highPath:
            self.addFlow(ij, -bot)
            
        #print("after algB shift:", bot, maxshift, self.getTT(0, highPath) - self.getTT(0, lowPath))
        
    def printFlows(self):
        self.network.dijkstras(self.origin, self.network.type)
        
        cost = {}
        
        for ij in self.network.links:
            cost[ij] = ij.end.cost
        
        self.minUsedTree()
        self.maxUsedTree()
        
        if not self.checkEquilibrium():
            for i in self.sorted:
                for ij in i.outgoing:
                    if self.contains(ij):
                        costdiff = ij.end.cost-cost[ij]
                        if costdiff > 0.1:
                            print("\t", ij, "%.2f, %.2f, %.2f, %.2f" % (self.getFlow(ij), ij.getTravelTime(ij.x, self.network.type), costdiff, ij.end.cost2 - ij.end.cost))

    
    def improveBush(self):
        self.minUsedTree()
        
        for ij in self.network.links:
            if ij.getReducedCost() < 0:
                # add link to bush
                self.flow[ij] = 0
        
        self.removeCycles()
            
    
    def minUsedTree(self):
        for u in self.sorted:
            u.cost = Params.INFTY*1000
            u.pred = None

        self.origin.cost = 0
        
        for u in self.sorted:
            for uv in u.outgoing:
                if self.getFlow(uv) < self.network.params.flow_epsilon:
                    continue    
                
                v = uv.end
                temp = uv.getTravelTime(uv.x, self.network.type) + u.cost
                
                if temp < v.cost:
                    v.cost = temp
                    v.pred = uv
                    #print(v.pred)

        
        output = {}
        
        for n in self.sorted:
        
            output[n] = n.pred
        
        
        return output
    
    def maxUsedTree(self):
        for u in self.sorted:
            u.cost2 = 0
            u.pred2 = None

        self.origin.cost2 = 0
        
        for u in self.sorted:
            for uv in u.outgoing:
                if self.getFlow(uv) < self.network.params.flow_epsilon:
                    continue    
                
                v = uv.end
                temp = uv.getTravelTime(uv.x, self.network.type) + u.cost2
                
                if temp > v.cost2:
                    v.cost2 = temp
                    v.pred2 = uv
                    #print(v.pred2)

        
        output = {}
        
        for n in self.sorted:
        
            output[n] = n.pred2
        
        
        return output

    

    
    def tracePath(self, i, j):
        output = []
        
        curr = j
        #print("tracing", i, j)
        
        while curr != i:
            if curr.pred == None:
                return None
            output.append(curr.pred)
            #print(curr, curr.pred, i, j, self.origin, output)
            curr = curr.pred.start
        
        return output
    
    def tracePath2(self, i, j):
        output = []
        
        curr = j
        
        while curr != i:
            if curr.pred2 == None:
                return None
            #print(str(curr)+" "+str(curr.pred2)+" "+str(i)+" "+str(j))
            output.append(curr.pred2)
            curr = curr.pred2.start
        
        return output
    
    def tracePathSet(self, i, j):
        output = []
        
        curr = j
        
        while curr != i:

            output.append(curr.pred)
            #print(curr.pred)
            curr = curr.pred.start
        
        return output
    
   
    
    def minPath(self):
        for u in self.network.nodes:
            u.cost = Params.INFTY
            u.pred = None
        
        self.origin.cost = 0
        
        
        for u in self.sorted:
            
            for uv in u.getBushOutgoing(self):
                v = uv.end
                
                temp = uv.getTravelTime(uv.x, self.network.type) + u.cost
                
                if temp < v.cost:
                    v.cost = temp
                    v.pred = uv
               

        output = {}
        
        for n in self.sorted:
            output[n] = n.pred
        
        
        return output
   
    
    def branchShifts(self):

        if self.network.params.PRINT_TAPAS_INFO and len(self.branches) > 0:
            print("starting branch shifts", self.origin.id)
        
        for b in self.branches:
            if self.network.params.PRINT_BRANCH_INFO:
                print("\tshift branch with flow",self.getFlow(b.endlink))
                
            if self.getFlow(b.endlink) > self.network.params.flow_epsilon:
                if b.init():
                    b.flowShift(self.network.type)
                #print(b)
                #return(b.flowShift())
            
        
        self.branches = []
    
  
    
    
    #def validateDemand(self):

        #for s in self.sorted:
            #print(s)
            #if s != self.origin and isinstance(s, self.Dest):
                #d = self.origin.getDemand(s)
                
                #actual = 0
                
                #for id in s.incoming:
                    #actual += self.flow[id]
                
                
                #for dj in sorted(s.outgoing, key=lambda edge: edge.end.id):
                    #actual -= self.flow[dj]
                
                
                #if abs(d - actual) > Params.flow_epsilon:
                    #print("Origin "+str(self.origin)+": demand to "+str(s)+" is "+str(d)+" but flow is "+str(actual))
                    #return False
        #return True
    
    
    def addFlow(self, l, x):
        curr_flow = self.getFlow(l)
        #with open('result121.txt', 'a') as file, contextlib.redirect_stdout(file):
        if curr_flow + x < -self.network.params.flow_epsilon:
            raise("attempt to add flow to link "+str(l.id)+" from bush "+str(self.origin)+" of flow "+str(self.flow[l])+" "+str(x))

            #print(f"the l.x is {l.x}")
        l.x += x
            #print(f"The l.x after adding x is {l.x} and x is {x}")
        curr_flow += x
        
        if curr_flow > self.network.params.flow_epsilon:
            self.flow[l] = curr_flow
        elif l in self.flow:
            del self.flow[l]
            #print(f"The self.flow of l is {self.flow[l]}")

        
    def checkPAS(self):
        
        minPathTree = self.network.getSPTree(self.origin)
        #print(self.origin)
        
        # only create 1 PAS per link per origin per iteration
        included = set()
        
        # look for all used links not part of the tree of least cost routes
        # search in backwards topological order
        
        
        #for l in self.flow:
            #print(l, self.origin, self.getFlow(l), l.getReducedCost(self.network.type), l.hasHighReducedCost(self.network.type, self.network.params.bush_gap))
            
        #with open('result171.txt', 'a') as file, contextlib.redirect_stdout(file):
        for n in self.sorted[::-1]:
            #print(n)

            #for l in n.incoming:
            #print(f"n.incoming is {n.incoming}")
            #for l in sorted(n.incoming, key=lambda l: (l.start.id, l.end.id)): 
            for l in n.incoming:
                #print(sorted(n.incoming, key=lambda l: (l.start.id, l.end.id)))
                #print(f"l before {l}")
                #print(n.incoming)
                
                #if l.end != self.origin and self.getFlow(l) > self.network.params.flow_epsilon and l.hasHighReducedCost(self.network.type, self.network.params.bush_gap):
                #    if self.origin.id == 14:
                #        print(l in included, included)
                
                #if self.getFlow(l) > self.network.params.flow_epsilon:
                #    print("\tlink", l, self.origin, self.getFlow(l), l.getReducedCost(self.network.type), l.hasHighReducedCost(self.network.type, self.network.params.bush_gap))

                # check for links with high reduced cost and positive flow, not just links not on the shortest path
                if l not in included and l.end != self.origin and self.getFlow(l) > self.network.params.bush_gap * self.origin.totaldemand and l.hasHighReducedCost(self.network.type, self.network.params.bush_gap):
                    
                    #print(f"included is {included}")
                    #print(f"l.end is {l.end} and l.start is {l.start}")
                    #print(f"self.flow of l is {self.flow}")
                    #print(f"self.network.params.flow_epsilon is {self.network.params.flow_epsilon}")
                    #print(f"the last thing {l.hasHighReducedCost(self.network.type, self.network.params.pas_cost_mu)}")
                    #print(l)
                    #print(f"selforigin is {self.origin}")
                    #System.out.println(l+" "+l.getDest().cost+" "+l.getSource().cost+" "+(l.getDest().cost-l.getSource().cost)+" "+l.getTT());
                    # we need a PAS!

                    if self.network.params.PRINT_PAS_INFO:
                        print("\tneed PAS for link", l, self.origin, self.getFlow(l), l.getReducedCost(self.network.type), l.hasHighReducedCost(self.network.type, self.network.params.bush_gap))

 

                    if not self.hasRelevantPAS(l):

                        # should check if we can borrow one from network
                        fromNetwork = self.network.findPAS(l, self)
                        if fromNetwork is None:
                            self.network.params.good_bush_gap += 1
                            
                            if self.network.params.PRINT_PAS_INFO:
                                print("\n\tCreate PAS for " + str(l) + " for origin " + str(self.origin))

                            newPAS = None

                            if l.start != self.origin:
                                newPAS = self.createPAS(minPathTree, l, max(self.network.params.flow_epsilon, self.getFlow(l) * self.network.params.pas_flow_mu))

                            # this was causing a crash??
                            else:
                                newPAS = self.createOriginLinkPAS(minPathTree, l)

                            #print(l)
                            
                            

                            #if newPAS != None:
                                #print(l.getReducedCost(self.network.type))
                                #print(newPAS.isFlowEffective(self.network.params.pas_flow_mu, self.network.params.pas_flow_mu * self.origin.getProductions(), self.network.type), newPAS.isCostEffective(self.network.type, self.network.params.pas_cost_mu))
                                #print(self.getFlow(l))
                                
                            #newPAS = None
                            
                            if l.start != self.origin and (newPAS is None or not newPAS.isEffective(self.network.type, self, self.network.params.pas_cost_mu, self.network.params.pas_flow_mu, self.network.params)):
                                
                                # branch shift
                                if self.network.params.PRINT_BRANCH_INFO:
                                    print("\t\tbranch shift!")

                                branch = self.createBranch(l)
                                if branch is not None:
                                    self.branches.append(branch)
                                else:
                                    if self.network.params.PRINT_BRANCH_INFO:
                                        print("branch is null")

                            else:
                                
                                self.savePAS(newPAS)
                                for ij in newPAS.backwardlinks:
                                    included.add(ij)
                        else:
                            if self.network.params.PRINT_PAS_INFO:
                                print("\t\tSAVE PAS")

                            fromNetwork.addRelevantOrigin(self.origin)
                            self.relevantPAS.add(fromNetwork)
                            #print(self.origin)
                    else:
                        if self.network.params.PRINT_PAS_INFO:
                        #if self.origin.id == 21:
                            print("\thas relevant PAS already")


    

    def createBranch(self, endlink):
        minpath = self.tracePathSet(self.origin, endlink.end)
        
        output = Branch.Branch(self, endlink, minpath)
        
        return output
   
    
    def hasRelevantPAS(self, a):
        #with open('hasrelevantPAS4.txt', 'a') as file, contextlib.redirect_stdout(file):
            #print(a)
            
            
                        
            if a in self.relevantPAS.forward:
                #print(f"a is {a}")
                for p in self.relevantPAS.forward[a]:
                    #if p.isCostEffective(self.network.type, self.network.params.pas_cost_mu) and p.maxForwardBushFlowShift(self) > self.network.params.pas_flow_mu * self.origin.getProductions():
                    if p.isEffective(self.network.type, self, self.network.params.pas_cost_mu, self.network.params.pas_flow_mu * self.origin.getProductions(), self.network.params): 
                        #if self.origin.id == 21:
                            #print(p, "maxshift", str(p.relevant), p.maxForwardBushFlowShift(self), "check with", self.network.params.pas_flow_mu * self.origin.getProductions(), p.isCostEffective(self.network.type, self.network.params.pas_cost_mu))
                        p.relevant.add(self.origin)
                        return True
            
            if a in self.relevantPAS.backward:
                #print(a)     
                for p in self.relevantPAS.backward[a]:
                    #if p.isCostEffective(self.network.type, self.network.params.pas_cost_mu) and p.maxBackwardBushFlowShift(self) > self.network.params.pas_flow_mu * self.origin.getProductions():
                    if p.isEffective(self.network.type, self, self.network.params.pas_cost_mu, self.network.params.pas_flow_mu * self.origin.getProductions(), self.network.params):    
                        #if self.origin.id == 21:
                        #    print(p, "maxshift", str(p.relevant), p.maxBackwardBushFlowShift(self), "check with", self.network.params.pas_flow_mu * self.origin.getProductions(), p.isCostEffective(self.network.type, self.network.params.pas_cost_mu))
                        p.relevant.add(self.origin)    
                        return True
                    
            return False
            
    
            
    def removeCycles(self):

        
        # right now this restarts the entire loop when a cycle is detected. I think we don't need to restart everything...
        cycleDetected = True
        while cycleDetected:
            for n in self.network.nodes:
                n.visited = False
                n.pred2 = None
                n.top_order = -1
                
                for l in n.incoming:
                    l.visit_order = -1
            
            cycleDetected = False
            
            sorted = list()
            
            idx = len(self.network.nodes)-1
            
            #unvisited = [] # this is a stack
            #unvisited.append(self.origin)
            unvisited = [self.origin]
            #print(self.origin)

            
            while len(unvisited) > 0:
                #unvisited = sorted(unvisited, key=lambda n: n.node.id)
                n = unvisited.pop()
                
                
                
                if isinstance(n, Node.Node):

                    if n.top_order >= 0:
                        continue
                    elif n.visited:
                        # remove the cycle
                        
                        #print("cycle check "+str(n)+" for origin "+str(self.origin))
                        
                        #for nb in self.network.nodes:
                        #    if nb.pred2 is not None:
                        #        print("\t"+str(nb)+" "+str(nb.pred2)+" "+str(self.flow[nb.pred2]))
                        
                        self.removeCycleAtNode(n)

                        cycleDetected = True
                        break
                    else:
                        n.visited = True

                        unvisited.append(NodeReturn.NodeReturn(n))

                        #outgoing_edges_list = list(n.outgoing)
                        #outgoing_edges_list.sort(key=lambda e: (e.start.id, e.end.id))
                        #for l in n.outgoing:
                        #for l in sorted(n.outgoing, key=lambda edge: edge.end.id):
                        for l in n.outgoing:
                            if self.contains(l):
                                j = l.end
                                j.pred2 = l
                                unvisited.append(j)
                else:
                    node = n.node
                    if node.top_order < 0:
                        sorted.append(node)
                        node.top_order = idx
                        idx -= 1


    def removeCycleAtNode(self, n):
        # n is the root node of the cycle
        
        #for l in self.flow:
        #    print(l, l.start.pred2)
            
        list = []
        
        curr = n
        
        
        start_idx = 0
        
        idx = 0
        #print(self.origin, self.flow)
        loopbreak = True
        
        while True:
            
            pred = curr.pred2
            
            if pred.visit_order >= 0:
                start_idx = pred.visit_order
                break
            else:
                pred.visit_order = idx
                idx += 1
                #if pred is None:
                    #print(f"Error: Attempt to access 'start' attribute of None at node {curr}")
                    #break  #
                list.append(pred)

                curr = pred.start
           
            
        
        n.pred2 = None
        
        maxflow = Params.INFTY
        #print(maxflow)
        
        #with open('result119.txt', 'a') as file, contextlib.redirect_stdout(file):
        for l in list[start_idx:]:
            maxflow = min(maxflow, self.getFlow(l))
            #print(f"maxflow is {maxflow}")
            #print(f"self.flow[l] is {self.flow[l]}")
        
        
        #print("removing "+str(maxflow)+" from "+str(list))
        #with open('result116.txt', 'a') as file, contextlib.redirect_stdout(file):
        for l in list[start_idx:]:
            self.addFlow(l, -maxflow)
            #print(-maxflow)

    def getPathAsNodeSet(self, minPathTree, i, j):
        #with open('result53.txt', 'a') as file, contextlib.redirect_stdout(file):
            output = set()
            
            output.add(j)
            
            curr = j
            
            while curr != i:
                curr = curr.pred.start
                output.add(curr)
                #print(curr)
                
            return output
        
    def createOriginLinkPAS(self, minPathTree, a):
        output = PAS.PAS()

            
        output.backwardlinks.append(a)
        
        for l in self.network.traceTree(minPathTree, a.start, a.end):
            #print(f'tree entry {minPathTree}, r {firstSimilar.start}, s {a.end}')
            output.forwardlinks.append(l)
            #print(l)
            
        output.start = a.start
        output.end = a.end

        if self.network.params.PRINT_PAS_INFO:
            print("\t\tPAS is "+str(output))
                
        return output
    
    # create a PAS for link a
    def createPAS(self, minPathTree, a, minflow):
        #with open('result149.txt', 'a') as file, contextlib.redirect_stdout(file):
            output = PAS.PAS()
            #print(a)
            # min path to a.dest
            minPath = self.getPathAsNodeSet(minPathTree, self.origin, a.end)
            #print(self.origin)
            if self.network.params.PRINT_PAS_INFO:
                print("\t\tminPath is "+str(minPath)+" for "+str(a.end))
            
            #with open('result147.txt', 'a') as file, contextlib.redirect_stdout(file):
        # store trace to avoid repeating breadth first search
            trace = {}
            unvisited = deque([a])
            #print(f"a is {a}")
        
            #unvisited = sllist()
            #print(unvisited)
            
            #unvisited.append(a)
            #print(f"unvisited is {unvisited}")
            
            firstSimilar = None
            
            #with open('result145.txt', 'a') as file, contextlib.redirect_stdout(file):
            while len(unvisited) > 0:
                jk = unvisited.popleft()
                #print(f"jk is {jk}")
                #print("consider "+str(jk)+" "+str(unvisited))
                j = jk.start
                
                if j in minPath:
                    firstSimilar = jk
                    #print(jk)
                    break
                
                for ij in j.incoming:
                    #print(str(ij.start)+" "+str(ij.end)+" "+str(self.flow[ij])+" "+str(minflow))
                    if self.getFlow(ij) > minflow and ij.start.top_order < j.top_order:
                    #if self.flow[ij] > minflow:
                        unvisited.append(ij)
                        trace[ij] = jk
            
            
            
            
            
            #if self.network.params.PRINT_PAS_INFO:
                #print("\tfirstSimilar is "+str(firstSimilar))
                #print(str(trace)+" "+str(minflow)+" "+str(self.getFlow(a))+" "+str(firstSimilar))
            
            if firstSimilar is None:
                return None
            

            
            # trace firstSimilar to a in used flow bush: this is the backward side of the PAS
            # if this loop breaks it's probably because we shifted flow and tried to find a PAS after shifting flow before re-running shortest path
            #with open('result24.txt', 'a') as file, contextlib.redirect_stdout(file):
            curr = firstSimilar
            #print(curr)
            output.backwardlinks.append(firstSimilar)
            while curr != a:
                
                #System.out.println("\t"+curr);
                curr = trace[curr]
                output.backwardlinks.append(curr)
                
            
            
            # trace firstSimilar to a in min path tree: this is the forward side of the PAS
            for l in self.network.traceTree(minPathTree, firstSimilar.start, a.end):
                #print(f'tree entry {minPathTree}, r {firstSimilar.start}, s {a.end}')
                output.forwardlinks.append(l)
                #print(l)
                    
            if len(output.forwardlinks) == 0:
                return None
            
            output.start = firstSimilar.end
            output.end = a.end
            
            if self.network.params.PRINT_PAS_INFO:
                print("\t\tPAS is "+str(output))
            
            output.addRelevantOrigin(self.origin)
            #print(self.origin)
            return output
            
    def savePAS(self, pas):
        self.network.allPAS.add(pas)
        self.relevantPAS.add(pas)
        
           
    
    def equilibrateDemand(self):
        #with open('result130.txt', 'w') as file, contextlib.redirect_stdout(file):
            self.network.dijkstras(self.origin, self.network.type)
            #with open('result37.txt', 'w') as file, contextlib.redirect_stdout(file):
            for s in self.origin.destSet:
            
                path = self.tracePath(self.origin, s)
                dem = self.demand[s]
                
                tt = self.calcPathTT(path, 0)
                pathflow = self.calcPathFlow(path)

                bot = -min(pathflow, dem)
                top = self.origin.demandFunc(s, tt)
                #print("starting with ", top, tt)



                while top -bot > self.network.params.line_search_gap:
                    mid = (bot+top)/2

                    tt_path = self.calcPathTT(path, mid)
                    tt_dem = self.origin.demandFuncInv(s, dem + mid)

                    #print("\t", dest, tt_path-tt_dem, bot, top)

                    # if path TT > dem TT, we should reduce dem to equilibrate
                    # if path TT < dem TT, we could support more flow
                    if tt_path - tt_dem < -self.network.params.line_search_gap:
                        bot = mid
                    elif tt_path - tt_dem > self.network.params.line_search_gap:
                        top = mid
                    else:
                        bot = mid
                        break


                
                #d = self.origin.getDemand(s)
                d = bot
                
                    
                #print("\t", self.origin, "add", d)   
                    
                for ij in path:
                    self.addFlow(ij, d)
                    
                self.demand[s] += d
            self.topologicalSort()
        #print(self.topologicalSort)        
                
        
    def calcPathFlow(self, path):
        output = 1e10
        
        for ij in path:
            output = min(output, self.getFlow(ij))
        
        return output
        
    def calcPathTT(self, path, shift):
        output = 0
        
        for ij in path:
            output += ij.getTravelTime(ij.x + shift, self.network.type)
        
        return output
            
    def __str__(self):
        return "bush "+str(self.origin.id)