# Created on : May 18, 2024, 10:39:11 AM
# Author     : michael

from src import Heap
from src import Node

def test():
    
    # Driver program to test above functions. 

    # Assuming the maximum size of the heap to be 15. 
    h = Heap.Heap() 

    # Asking the user to input the keys: 
    k, i, n = 6, 0, 6
    print("Entered 6 keys: 3, 10, 12, 8, 2, 14 \n") 
    n1 = Node.Node(1)
    n1.cost = 3
    n2 = Node.Node(2)
    n2.cost = 10
    n3 = Node.Node(3)
    n3.cost = 12
    n4 = Node.Node(4)
    n4.cost = 3
    n5 = Node.Node(5)
    n5.cost = 5
    n6 = Node.Node(6)
    n6.cost = 14

    h.insert(n1)
    h.insert(n2)
    h.insert(n3)
    h.insert(n4)
    h.insert(n5)
    h.insert(n6)
    
    
    n5.cost = 3
    h.insert(n5)

    # Printing the current size 
    # of the heap. 
    print("The current size of the heap is "
            + str(h.size()) + "\n") 

    # Printing the root element which is 
    # actually the maximum element. 
    print("The current minimum element is " + str(h.getMin().id)+" "+str(h.getMin().cost) 
            + "\n") 

    # Deleting key at index 2. 
    n = h.removeMin() 

    # Printing the size of the heap 
    # after deletion. 
    print("After deleting "+str(n.cost)+" the current size of the heap is "
            + str(h.size()) + "\n") 

    print("The current minimum element is " + str(h.getMin().id) +" "+str(h.getMin().cost) 
            + "\n") 

    n = h.removeMin() 
    
    print("After deleting "+str(n.cost)+" the current size of the heap is "
            + str(h.size()) + "\n") 

    print("The current minimum element is " + str(h.getMin().id) +" "+str(h.getMin().cost) 
            + "\n") 
            
    n6.cost = 1
    h.insert(n6)

    print("After changing "+str(n6.cost)+" the current size of the heap is "
            + str(h.size()) + "\n") 

    print("The current minimum element is " + str(h.getMin().cost) 
            + "\n") 

    n4.cost = 0
    h.insert(n4)

    print("After changing "+str(n4.cost)+" the current size of the heap is "
            + str(h.size()) + "\n") 

    print("The current minimum element is " + str(h.getMin().cost) 
            + "\n") 

    # Inserting 2 new keys into the heap.
    n7 = Node.Node(7)
    n7.cost = 15
    n8 = Node.Node(9)
    n8.cost = 5 
    h.insert(n7) 
    h.insert(n8) 
    print("After adding "+str(n7.cost)+" and "+str(n8.cost)+" the current size of the heap is "
            + str(h.size()) + "\n") 
    print("The current minimum element is " + str(h.getMin().cost) 
            + "\n")

    while h.size() > 0:
        print(h.removeMin().cost)


    raise Exception("terminate")
