class Node:

    # construct this Node with the given id
    def __init__(self, id):
        # used for Dijkstra's implementation
        self.cost = 0.0
        self.cost2 = 0.0
        
        self.id = id
        self.outgoing = []
        self.incoming = []
        self.outgoing_edges_list = []
        self.top_order = -1
        self.visited = False
        self.in_degree = 0
        
        self.heap_idx = -1    

        # pred2 is used by bush when not finding shortest paths
        self.pred = None
        self.pred2 = None
        
    def __repr__(self):
        return str(self)
        
    def getBushOutgoing(self, b):
            output = []
            for l in sorted(self.outgoing, key=lambda edge: edge.end.id):
                if b.contains(l):
                    output.append(l)
            return output
    
    def getBushIncoming(self, b):    
        output = []        
        for l in self.incoming:
            if b.contains(l):
                output.append(l)        
        return output
    

        
    # returns a list of links containing the outgoing links of this node

    # returns True if this node is a thru node
    def isThruNode(self):
        return True
  
    # returns the id of this node
    def getId(self):
        return self.id
    
    def __str__(self):
        return str(self.id)

    def addOutgoingLink(self, ij):
        #with open('result39.txt', 'a') as file, contextlib.redirect_stdout(file):
            self.outgoing.append(ij)
            #print(ij)
    
    def addIncomingLink(self, ij):
        #with open('result102.txt', 'a') as file, contextlib.redirect_stdout(file):
            self.incoming.append(ij)
            #print(ij)

    def __lt__(self, other):
        return self.id < other.id