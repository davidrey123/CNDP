from src import Node
from src import Link
from src import Path
from src import Zone
from src import Params
from src import Heap
import random

class Network:

    # construct this Network with the name; read files associated with network name
    def __init__(self,name, num_candidate, B_):
        self.nodes = [] 
        self.links = []
        self.zones = []
        self.origins = []
        
        random.seed(10)
        
        self.name = name
        
        self.params = Params.Params()
        
    

            
        self.readNetwork("data/"+name+"/net.txt",1,1)
        self.readTrips("data/"+name+"/trips.txt",1,1,1)
        
        
        for r in self.origins:
            r.setDests()
        
        self.B = B_
        
        self.candidates = []
        
        for a in self.links:
            a.enabled = False
        
        
        for i in range (0, num_candidate):
            while True:
                link = self.links[random.randint(0, len(self.links))]
                if link not in self.candidates:
                    self.candidates.append(link)
                    break

        for a in self.candidates:
            a.enabled = True
        
    # read file "/net.txt"
    def readNetwork(self,netFile,scal_time,scal_flow):
        
        firstThruNode = 1
        numZones = 0
        numNodes = 0
        numLinks = 0
        newLinks = 0
        file = open(netFile, "r")

        line = ""
        
        while line.strip() != "<END OF METADATA>":
            line = file.readline()
            if "<NUMBER OF ZONES>" in line:            
                numZones = int(line[line.index('>') + 1:].strip())            
            elif "<NUMBER OF NODES>" in line:
                numNodes = int(line[line.index('>') + 1:].strip())
            elif "<NUMBER OF LINKS>" in line:
                numLinks = int(line[line.index('>') + 1:].strip())
            elif "<NUMBER OF NEW LINKS>" in line:
                newLinks = int(line[line.index('>') + 1:].strip())
            elif "<FIRST THRU NODE>" in line:
                firstThruNode = int(line[line.index('>') + 1:].strip())

        for i in range(0, numZones):
            self.zones.append(Zone.Zone(i + 1))

        for i in range(0, numNodes):
            if i < numZones:
                self.nodes.append(self.zones[i])
                
                if i + 1 < firstThruNode:
                    self.zones[i].setThruNode(False)

            else:
                self.nodes.append(Node.Node(i + 1))

        line = ""
        id = 0
        
        
        while len(line) == 0:
            line = file.readline().strip()

        for i in range(0, numLinks + newLinks):
            line = file.readline().split()
            if len(line) == 0:
                continue
            start = self.nodes[int(line[0]) - 1]
            end = self.nodes[int(line[1]) - 1]
            C = float(line[2]) * scal_flow   
            
            length = float(line[3])


            t_ff = float(line[4]) * scal_time
            alpha = float(line[5])
            beta = float(line[6])
            
   
            try:
                cost = float(line[10])
            except ValueError:
                cost = 0
                
            #print(start, end, cost, line)

            foundLink = False
            for ij in start.outgoing:
                if ij.end == end:
                    foundLink = True
                    break
            
            if not foundLink:
                link = Link.Link(id, start ,end, length)
                id = id +1
                #print(start,end)
                self.links.append(link)
                
                link2 = Link.Link(id, end ,start, length)
                id = id +1
                #print(start,end)
                self.links.append(link2)
                #print(self.links)
            
            if i >= numLinks:
                self.links2.append(link)
            
        file.close()

        print("num links", len(self.links))


    def readTrips(self,tripsFile,scal_time,scal_flow,inflate_trips):
        
        file = open(tripsFile, "r")
        
        lines = file.readlines()
        
        line_idx = 0
        
        while lines[line_idx].strip() != "<END OF METADATA>":
            line_idx += 1
            
        line_idx += 1
        
        while lines[line_idx].strip() == "":
            line_idx += 1
            
        r = None
        
        idx = 0
        
        splitted = lines[line_idx].split()
        #print(splitted)
        
        while len(lines) < line_idx or idx < len(splitted):

            next = splitted[idx]

            if next == "Origin":
                
                idx += 1
                r = self.zones[int(splitted[idx]) - 1]

            else:
                s = self.zones[int(splitted[idx]) - 1]

                #print(s)
                idx += 2
                next = splitted[idx]
                d = float(next[0:len(next) - 1]) * scal_flow
                d = d * inflate_trips
                
                
                
                r.addDemand(s, d)


            idx += 1

            if idx >= len(splitted):
                line_idx += 1
                while line_idx < len(lines) and lines[line_idx].strip() == "":
                    line_idx += 1
                    
                if line_idx < len(lines):
                    line = lines[line_idx].strip()
                    splitted = line.split()
                    idx = 0
            
        file.close()
        
        for r in self.zones:
            if r.getProductions() > self.params.flow_epsilon:
                self.origins.append(r)

    def getLinks(self):
        return self.links
    
    def getNodes(self):
        return self.nodes
    
    def getZones(self):
        return self.zones

    # find the node with the given id
    def findNode(self, id):
        if id <= 0 or id > len(self.nodes):
            return None
        return self.nodes[id - 1]

    # find the link with the given start and end nodes
    def findLink(self, i_in, j_in):
        
        i = None
        j = None
        
        if isinstance(i_in, Node.Node):
            i = i_in
        elif isinstance(i_in, int):
            i = self.findNode(i_in)
            
        if isinstance(j_in, Node.Node):
            j = j_in
        elif isinstance(j_in, int):
            j = self.findNode(j_in)   
        
        if i is None or j is None:
            return None
            

        for link in i.outgoing:
            if link.end == j:
                return link

        return None
    

    

    def dijkstras(self, origin, max_cost, type):
        
            for n in self.nodes:
                n.cost = Params.INFTY
                n.pred = None

            origin.cost = 0.0

            Q = Heap.Heap()
            Q.insert(origin)
            


            while Q.size() > 0:

                u = Q.removeMin()
                
                #if type == 'RC':
                #    print('u',u)

                for uv in u.outgoing:
                    v = uv.end
                    tt = uv.getCost(type)

                    #if u.cost + tt < v.cost:
                    if u.cost + tt < v.cost and v.cost - u.cost - tt >= self.params.SP_tol and u.cost + tt <= max_cost:
                        v.cost = u.cost + tt
                        v.pred = uv
                        
                        #if type == 'RC':
                        #    print('v',v,v.pred,v.cost)

                        if v.isThruNode():
                            Q.insert(v)
            

    def trace(self, r, s):
        curr = s

        output = Path.Path()
        output.r = r
        output.s = s
        
        while curr != r and curr is not None:
            ij = curr.pred

            if ij is not None:
                output.add(ij)
                curr = curr.pred.start
              
        #print('trace',r,s,output)
              
        return output
        
    
