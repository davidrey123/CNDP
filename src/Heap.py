from src import Params

class Heap:

    # Number of elements in the 
    # Max heap currently. 
    heapSize = 0
    maxSize = 0

    # Constructor function. 
    def __init__(self): 
        self.arr = []
        self.heapSize = 0
        self.maxSize = 0
        
    def compare(self, i, j):
        if abs(i.cost - j.cost) > Params.Params().SP_tol:
            return i.cost - j.cost
        else:
            return i.id - j.id

    # Heapifies a sub-tree taking the 
    # given index as the root. 
    def Heapify(self, i): 
        
        l = self.lChild(i) 
        r = self.rChild(i) 
        smallest = i 
        if l < self.heapSize and self.compare(self.arr[l], self.arr[i]) < 0: 
            smallest = l 
        if r < self.heapSize and self.compare(self.arr[r], self.arr[smallest]) < 0: 
            smallest = r 
        if smallest != i: 
            temp = self.arr[i] 
            self.arr[i] = self.arr[smallest] 
            self.arr[smallest] = temp 
            
            self.arr[i].heap_idx = i
            temp.heap_idx = smallest
            
            self.Heapify(smallest) 
        
        '''

        while True:
            l = self.lChild(i) 
            r = self.rChild(i) 
            
            smallest = i 
            if l < self.heapSize and self.arr[l].cost < self.arr[i].cost: 
                smallest = l 
            if r < self.heapSize and self.arr[r].cost < self.arr[smallest].cost: 
                smallest = r 
            if smallest != i: 
                temp = self.arr[i] 
                self.arr[i] = self.arr[smallest] 
                self.arr[smallest] = temp 

                self.arr[i].heap_ix = i
                temp.heap_idx = smallest
                
                i = smallest
            else:
                break
        '''   
    # Returns the index of the parent 
    # of the element at ith index. 
    def parent(self, i): 
        return (i - 1) // 2

    # Returns the index of the left child. 
    def lChild(self, i): 
        return (2 * i + 1) 

    # Returns the index of the 
    # right child. 
    def rChild(self, i): 
        return (2 * i + 2) 

    # Removes the root which in this 
    # case contains the maximum element. 
    def removeMin(self): 
        
        # Checking whether the heap array 
        # is empty or not. 
        if self.heapSize <= 0: 
            return None
        if self.heapSize == 1: 
            self.heapSize -= 1
            self.arr[0].heap_idx = -1
            return self.arr[0] 

        # Storing the maximum element 
        # to remove it. 
        root = self.arr[0] 
        
        # validate
        #for i in range(0, self.heapSize):
        #    if root.cost > self.arr[i].cost:
        #        raise Exception("heap failed "+str(root.cost)+" "+str(self.arr[i].cost))
        
        self.arr[0] = self.arr[self.heapSize - 1] 
        self.arr[0].heap_idx = 0
        self.heapSize -= 1
        
        # To restore the property 
        # of the Max heap. 
        self.Heapify(0) 
        
        root.heap_idx = -1
        

        return root 

    # Increases value of key at 
    # index 'i' to new_val. 
    def decreaseKey(self, node): 
        
        i = node.heap_idx
        
        if self.arr[i].id != node.id:
            raise Exception("bad order for "+str(node.id))
        

        while i != 0 and self.compare(self.arr[self.parent(i)], self.arr[i]) > 0: 

            parent_idx = self.parent(i)
            
            temp = self.arr[i] 
            self.arr[i] = self.arr[parent_idx] 
            self.arr[parent_idx] = temp 
            
            self.arr[i].heap_idx = i
            temp.heap_idx = parent_idx
            
            i = parent_idx

    # Returns the maximum key 
    # (key at root) from max heap. 
    def getMin(self): 
        return self.arr[0] 

    def size(self): 
        return self.heapSize 

    def printHeap(self):
        for i in range(0, self.heapSize):
            print("\t", self.arr[i].id, self.arr[i].heap_idx)


    # Inserts a new key 'x' in the Max Heap. 
    def insert(self, x): 

        # if x is already in heap, then decreaseKey instead of inserting a copy
        if x.heap_idx >= 0:
            self.decreaseKey(x)
        else: 
            # To check whether the key 
            # can be inserted or not. 

            # The new key is initially 
            # inserted at the end. 
            self.heapSize += 1
            i = self.heapSize - 1
            if self.heapSize <= self.maxSize:
               self.arr[i] = x  
            else:
                self.arr.append(x)
                self.maxSize += 1
            x.heap_idx = i

            # The max heap property is checked 
            # and if violation occurs, 
            # it is restored. 
            while i != 0 and self.compare(self.arr[self.parent(i)], self.arr[i]) > 0: 
                parent_idx = self.parent(i)
                temp = self.arr[i] 
                self.arr[i] = self.arr[parent_idx] 
                self.arr[parent_idx] = temp 
                self.arr[i].heap_idx = i
                temp.heap_idx = parent_idx

                i = parent_idx 


