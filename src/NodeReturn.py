

# Created on : Mar 28, 2024, 10:25:40 AM
# Author     : michaellevin

# this is a wrapper class around a Node to indicate that I'm returning to it in the list after visiting it before
# I can't use a flag variable because it doesn't indicate the location in the stack

class NodeReturn:

    def __init__(self, node):
        self.node = node