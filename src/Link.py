from src import Params

class Link:


    # construct this Link with the given parameters
    def __init__(self, id, start, end, t_ff):
        self.id = id
        self.start = start
        self.end = end
        self.t_ff = t_ff

        self.enabled = True


        self.y = 0

        
        if start is not None:
            start.addOutgoingLink(self)
            
        if end is not None:
            end.addIncomingLink(self)


    def __repr__(self):
        return str(self)

    def getCost(self, type):
        if self.enabled and (type == True or self.y == 1):
            return self.t_ff
        else:
            return 1e6
        
    def __str__(self):
        return "(" + str(self.start.getId()) + ", " + str(self.end.getId()) + ")"
        
    
    def __hash__(self):
        return hash(self.id)