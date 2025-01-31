from src import Node


class Zone(Node.Node):


    def __init__(self, id):
        super().__init__(id)
        self.demand = {}
        self.totaldemand = 0
        self.thruNode = True
        self.bush = None
    

    # adds the specified demand to an internal data structure for the demand from this node to the destination
    def addDemand(self, dest, dem):
        if dest in self.demand.keys():
            self.demand[dest] = self.demand[dest] + dem 
        else:
            self.demand[dest] = dem
            
        self.totaldemand += dem
        
    def getProductions(self):
        return self.totaldemand
        
    # returns the number of trips from this node to the destination
    def getDemand(self, dest):
        if dest in self.demand.keys():
            return self.demand[dest]
        else:
            return 0
    


    # returns aboolean indicating whether this node is a thru node
    def isThruNode(self):
        return self.thruNode
    
    # set a boolean indicating whether this node is a thru node
    def setThruNode(self, thru):
        self.thruNode = thru