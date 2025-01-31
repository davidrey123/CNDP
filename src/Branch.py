# Created on : Mar 28, 2024, 2:08:29 PM
# Author     : michaellevin

class Branch:

    # this is just a placeholder that says this branch has a specific endlink and a minpath. The actual links in the branch will be determined later.
    def __init__(self, bush, endlink, minpath):
        self.bush = bush
        self.linkflows = {}
        
        self.endlink = endlink
        self.minpath = minpath
        self.maxflow = 0
        
        
    def init(self):
        for n in self.bush.network.nodes:
            n.visited = False
        
        
        unvisited = []
        
        branchlinks = []
        
        unvisited.append(self.endlink.start)
        
        while len(unvisited) > 0:
            j = unvisited.pop()
            #print(j)
            for ij in j.incoming:
 
                #print("\t"+str(ij.start)+" "+str(ij.end)+" "+str(self.bush.flow[ij])+" "+str(ij.x))
                if self.bush.contains(ij):
                    i = ij.start
                    
                    branchlinks.append(ij)
                    
                    
                    if not i.visited:
                        unvisited.append(i)
                        i.visited = True

        self.maxflow = self.bush.flow[self.endlink]
        #print(self.maxflow)
        
        #print("maxflow = "+str(self.maxflow)+" "+str(branchlinks)+" "+str(len(branchlinks))+" "+str(self.endlink))
            
        # now do Ford-Fulkerson to figure out branch flow on each link
        # the "capacities" are the bush flow on each link
        # due to conservation of flow I don't need to add flow in reverse. DFS will be sufficient.
        
        
        for l in branchlinks:
            self.linkflows[l] = 0.0
        
        
        start = self.bush.origin
        end = self.endlink.start
        self.linkflows[self.endlink] = self.maxflow
        #print(self.maxflow)
        
        
        assignedFlow = 0
        

        
        # while there is flow left to assign
        # use flow epsilon to avoid numerical error causing infinite loop
        while self.maxflow - assignedFlow > self.bush.network.params.flow_epsilon:
            #print(self.maxflow - assignedFlow)
            #print(str(self.maxflow)+" "+str(assignedFlow))
            
            # DFS find path
            unvisited = []
            unvisited.append(start)
            
            for n in self.bush.network.nodes:
                n.visited = False
                n.pred2 = None
            
            start.visited = True
            
            while len(unvisited) > 0:
                i = unvisited.pop()
                
                #print("evaluate "+str(i))
                # once DFS finds a path, stop and add flow. That path will become unusable
                if i == end:
                    break
                
                expanded = []
                for ij in i.outgoing:
                    # only expand links with positive bush flow - temporary branch flow
                    if ij in branchlinks and not ij.end.visited and self.bush.flow[ij] - self.linkflows[ij] > self.bush.network.params.flow_epsilon:
                        expanded.append(ij)
                    #if ij in branchlinks:
                        #print("\t"+str(ij.start)+" "+str(ij.end)+" "+str(self.bush.flow[ij] - self.linkflows[ij]))
                
                # sort in order of decreasing flow
                expanded.sort(key = lambda x: self.bush.flow[x] - self.linkflows[x])
                
                #print("expanded "+str(expanded)+" "+str(branchlinks)+" "+str(len(branchlinks)))
                
                #Collections.sort(expanded, new Comparator<Link>(){
                #    public int compare(Link i, Link j){
                #        double flowi = bush.getFlow(i) - linkflows.get(i);
                #        double flowj = bush.getFlow(j) - linkflows.get(j);
                #        return (int)Math.ceil(flowj - flowi);
                
                for ij in expanded:
                    j = ij.end
                    j.pred2 = ij
                    j.visited = True
                    unvisited.append(j)
            
            # trace path and label flows
            augmentedPath = self.bush.tracePath2(start, end)
            
            if augmentedPath == None:
                return False
                
                
            sendFlow = self.maxflow - assignedFlow
            
            for l in augmentedPath:
                sendFlow = min(sendFlow, self.bush.getFlow(l) - self.linkflows[l])
            
            
            for l in augmentedPath:
                self.linkflows[l] = self.linkflows[l] + sendFlow
            
            assignedFlow += sendFlow
     
        return True
           
    def flowShift(self, type):
        avgTT = self.getAvgTT(0, type)
        minTT = self.getMinTT(0, type)
        

        '''
        print("shortest path", self.bush.origin)
        self.bush.network.getSPTree(self.bush.origin)
        for n in self.bush.network.nodes:
            print(n, n.cost, type)
            for ij in n.outgoing:
                print("\t", ij, ij.x, ij.start.cost, ij.end.cost, ij.getTravelTime(ij.x, type))
        '''
        # difference is too small to be worth shifting
      
        if avgTT - minTT < minTT * self.bush.network.params.pas_cost_epsilon:
            
            #print("hi1", avgTT, minTT, avgTT-minTT, self.endlink.end.cost, self.bush.origin)
            #for l in self.minpath:
            #    print("\t", l, l.end.cost, l.getTravelTime(l.x, type))
            
            return 0
        
        bot = 0
        top = self.maxflow
        #print(top)
        
        while top - bot > self.maxflow * self.bush.network.params.line_search_gap:
            mid = (bot+top)/2
            
            newTTDiff = self.getAvgTT(mid, type) - self.getMinTT(mid, type)
            
            #System.out.println(bot+" "+mid+" "+top+" "+getAvgTT(mid)+" "+getMinTT(mid)+" "+newTTDiff);
            
            if newTTDiff > 0:
                # shift more
                bot = mid
            else:
                # shift less
                top = mid
                
        
        flowshift = self.propAddFlow(bot)
        
        #if self.bush.network.params.PRINT_BRANCH_INFO:
        #print("after branch shift is "+str(self.getAvgTT(0, type))+" "+str(self.getMinTT(0, type))+" "+str(flowshift))

        return bot
    
    def getMinTT(self, shift, type):
        if self.maxflow == 0:
        # Return a default value or handle the zero case appropriately
            return 0 
        output = 0
        
        # this is used in case the min links are also on the branch
        prop = shift/self.maxflow
        #print(prop)
        
        for l in self.minpath:
            newflow = l.x
            #print(newflow)
            
            # add flow to the minpath
            newflow += shift
            #print(newflow)
            
            # if link is in branch, some flow will be shifted:
            if l in self.linkflows:
                newflow -= prop * self.linkflows[l]
            
            output += l.getTravelTime(newflow, type)
        
        return output
        
    # consider shifting flow from branch to minpath
    
    def getAvgTT(self, shift, type):
        if self.maxflow == 0:
        # Return a default value or handle the zero case appropriately
            return 0 


        output = 0
        
        # proportional shift based on how much of the branch flow is on each link
        prop = shift/self.maxflow
  
        for l in self.linkflows:
            
             #subtract the flowshift
            
            flowchange = self.linkflows[l]*prop
            #print(l.x)
            newflow = l.x - flowchange
            #print(newflow)
            
            # if link is on minpath, then add the entire shift to it
            
            if l in self.minpath:
                newflow += shift

            output += (self.linkflows[l] - flowchange) * l.getTravelTime(newflow, type)
        
        return output / (self.maxflow - shift)
    
    def propAddFlow(self, shift):
        if self.maxflow == 0:
        # Return a default value or handle the zero case appropriately
            return 0 
        prop = shift/self.maxflow
        #print(self.maxflow)
        
        for l in self.linkflows:
            self.bush.addFlow(l, -self.linkflows[l]*prop)
            
            # this isn't needed: the branch will be discarded after equilibrating.
            if self.bush.network.params.DEBUG_CHECKS:
                self.linkflows[l] = self.linkflows[l] * (1-prop)

        for l in self.minpath:
            self.bush.addFlow(l, shift)
            
            # this isn't needed: the branch will be discarded after equilibrating.
            #if self.bush.network.params.DEBUG_CHECKS:
                #self.linkflows[l] = self.linkflows[l] + shift

        self.maxflow -= shift
        return shift