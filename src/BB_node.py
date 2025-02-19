from src import Params
from src import Node

class BB_node:
    def __init__(self, y_lb, y_ub, lb, ub, mu_lb, mu_ub):
        self.y_lb = y_lb
        self.y_ub = y_ub
        self.lb = lb
        self.ub = ub
        self.mu_lb = mu_lb
        self.mu_ub = mu_ub
   
