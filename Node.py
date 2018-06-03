class Node():
    """
    is a class that holdes the cities and their info.
    """
    def __init__(self):
        self.index = None
        self.x = None
        self.y = None
        self.demand = None

    def toString(self):
        print "Index: " + str(self.index) + "    x: " + str(self.x) + "     y: " + str(self.y) + "      demand: " + str(self.demand)
