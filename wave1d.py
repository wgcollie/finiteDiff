

class InitialCondition:

    def __init__(self,function,grid):

class GhostPoints:

    def __init__(self,grid):

class TimeIterator:

    def __init__(self,startTime,stopTime):

class ExplicitFiniteDifferenceScheme:

    def __init__(self,initialCondition,timeIterator):

    def setGhostPoints(self):

    def solver(self,tStart,tStop):
        timeIterator.setInterval(tStart,tStop)
        self.setGhostPoints
        for step in timeIterator:
             #Update all inner points
            up = 
             #Insert boundary conditions
            up[0]
             #Initiaalize for next step
            um = u; u = up
        
