import random

class Individual:
    def __init__(self, varlinks):
        
        self.data = {a: 0 for a in varlinks}
        self.fitness = 1e15
    
    def setData(self, data):
        self.data = data
        
    def binary_string(self, link, string_len):
        output = bin(self.data[link])[2:]
        
        while len(output) < string_len:
            output = "0"+output
                
        return output
    
    def calcY(self, varlinks, string_len):
        y = {a: self.data[a] / pow(2, string_len+1) for a in varlinks}
        
        return y
        
    def randomize(self, varlinks, string_len):
        # create random individual
        self.data = {a: random.randint(0, pow(2, string_len)-1) for a in varlinks}
        
