from src import Node
import math

class Zone(Node.Node):


    def __init__(self, id):
        super().__init__(id)
        self.demand = {}
        self.totaldemand = 0
        self.thruNode = True
        self.bush = None
        self.destSet = set()
        
        self.y = dict()
        self.a = dict()
        self.calTT = dict()
    
    
    def setY(self, ynew):
        for s in self.destSet:
            self.y[s] = ynew[(self,s)]
        
    def demandFunc(self, dest, tt):
        #print(self.y[dest], tt, self.a[dest], math.exp(-tt/self.a[dest]), self.getDemand(dest))
        
        output = self.y[dest] * math.exp(-tt/self.a[dest]) 
        #check = self.a[dest] * math.log(self.y[dest]/self.getDemand(dest)) 
        #print("check calc 2 ", output, check, tt, self.getDemand(dest))
        return output
        
    def demandFuncInv(self, dest, dem):
        return self.a[dest] * math.log(self.y[dest]/dem) 
        
    def Dinv(self, dest, q, y):
        return self.a[dest] * math.log(y / q)
        
    def intDinv(self, dest, q, y):
        if q == 0 or y == 0:
            #print(self.id, dest.id, self.getDemand(dest), q, y)
            return 0
        return self.a[dest] * q * (math.log(y/q) + 1)
        
    def intDerivDinv(self, dest, q, y):
        return self.a[dest] * q / y
    
    def demandFuncY(self, dest, dem, tt):
        y = math.exp(tt/self.a[dest]) * dem
        self.calTT[dest] = tt
        
        #print("check calc ",  dem * math.exp(-tt/self.a[dest]), self.a[dest], dem, tt, tt/self.a[dest] )
        #print("\t", y * math.exp(-tt*1.5/self.a[dest]), dem, tt*1.5)
        return y
    
        
    def getDests(self):
        return self.destSet
        
    def setDests(self):
        for dest in self.demand.keys():
            if self.demand[dest] > 0:
                self.destSet.add(dest)
  
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