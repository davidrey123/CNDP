from src import Params
from src import Node

class BB_node:
    def __init__(self, q_lb, q_ub, lb, ub, mu_lb, mu_ub):
        self.q_lb = q_lb
        self.q_ub = q_ub
        self.lb = lb
        self.ub = ub
        self.mu_lb = mu_lb
        self.mu_ub = mu_ub
        
        self.largest_diff = 0
        
        for w in q_lb.keys():
            self.largest_diff = max(self.largest_diff, q_ub[w] - q_lb[w])
   
