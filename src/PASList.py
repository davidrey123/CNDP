# Created on : Mar 27, 2024, 10:28:15 PM
# Author     : michael
import contextlib

class PASList:
    def __init__(self):
        self.forward = {}
        self.backward = {}
    
    def containsKey(self, ij):
        return ij in self.forward or ij in self.backward
        
    def size(self):
        return len(self.backward)
    
    def add(self, pas):
        #with open('result49.txt', 'a') as file, contextlib.redirect_stdout(file):
            ij = pas.getEndLinkBwd()
            #print(ij)
            if ij not in self.backward:
                self.backward[ij] = []
            #print(ij)
            #print(type(ij))
            self.backward[ij].append(pas)
            
            ij = pas.getEndLinkFwd()
            if ij not in self.forward:
                self.forward[ij] = []
            
            self.forward[ij].append(pas)
            #print(ij)


    def remove(self, pas):
        ij = pas.getEndLinkBwd()
        
        if ij in self.backward and pas in self.backward[ij]:
            self.backward[ij].remove(pas)
            
            if len(self.backward[ij]) == 0:
                del self.backward[ij]
        
        ij = pas.getEndLinkFwd()
        
        if ij in self.forward and pas in self.forward[ij]:
            self.forward[ij].remove(pas)
            
            if len(self.forward[ij]) == 0:
                del self.forward[ij]
    