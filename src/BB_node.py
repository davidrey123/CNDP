from src import Params
from src import Node

class BB_node:
    def __init__(self, y_lb, y_ub, lb, ll_gap, ub, mu_lb, mu_ub):
        self.y_lb = y_lb
        self.y_ub = y_ub
        self.lb = lb
        self.ub = ub
        self.mu_lb = mu_lb
        self.mu_ub = mu_ub
        self.ll_gap = ll_gap
        
        self.largest_diff = 0
        
        for w in y_lb.keys():
            self.largest_diff = max(self.largest_diff, y_ub[w] - y_lb[w])
   
