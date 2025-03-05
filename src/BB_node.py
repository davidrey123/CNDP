from src import Params
from src import Node

class BB_node:
    def __init__(self, y_lb, y_ub, local_lb, local_ub):
        self.y_lb = y_lb
        self.y_ub = y_ub
        self.lb = local_lb
        self.ub = local_ub

        
        self.largest_diff = 0
        
        for a in y_lb.keys():
            self.largest_diff = max(self.largest_diff, y_ub[a] - y_lb[a])
   