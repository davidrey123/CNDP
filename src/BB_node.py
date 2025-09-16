from src import Params
from src import Node

class BB_node:
    def __init__(self, y_lb, y_ub, local_lb, local_ub, parent, iter):
        self.y_lb = y_lb
        self.y_ub = y_ub
        self.lb = local_lb
        self.ub = local_ub
        
        self.iter = iter

        self.level = 0
        if parent is not None:
            self.level = parent.level+1
            
        self.left = None
        self.right = None
        self.parent = parent
        
        self.updated = False
        
        self.largest_diff = 0
        
        for a in y_lb.keys():
            self.largest_diff = max(self.largest_diff, y_ub[a] - y_lb[a])
            
    def updateLb(self, new_lb, tol):
        
        # this could happen from recursion
        if new_lb - self.lb > tol:
        
            self.updated = True
            #print("update", self.lb, new_lb)

            self.lb = new_lb

            if self.left is not None:
                self.left.updateLb(new_lb, tol)

            if self.right is not None:
                self.right.updateLb(new_lb, tol)
            
    def delChild(self, child):
        if self.left == child:
            self.left = None
        if self.right == child:
            self.right = None
   