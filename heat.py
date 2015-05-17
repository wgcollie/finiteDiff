import abc
from numpy import *

# One Dimensional Steady State Heat Conduction

#Parameters
problemInterval = [0,1]
nMeshPoints = 10
heatSource = "Exponential"
HeatEquationResultFile = "HeatEquationSolutionFile.dat"

#Constants
FIXED = 0
FLUX = 1

#Grid
class finiteDiff_1dGrid:

    def __init__(self,nMeshPoints,interval):
        print "1d Grid"
        self.nDim = 1
        self.interval = interval
        self.nMeshPoints = nMeshPoints
        self.meshSize = (interval[1]-interval[0])/nMeshPoints
        self.FIRST = 0
        self.LAST = nMeshPoints

    def setBoundaryProperties(self,pointID,Property):
        print "Set Boundary Properties"

    def __getitem__(self,k):
        x = self.interval[0] + k*self.meshSize

class finiteDiff_sources:

    def __init__(self,source):
        self.source = source

    def evaluateSourceContribution(self,x):
        print "HeatSource Contribution: %s"% self.source

#Solution
class finiteDiff_1dSolution:

    def __init__(self, nDim, nGridPoints):
        self.nDim = nDim
        self.nGridPoints = nGridPoints
        self.ss = zeros( nDim*nGridPoints )

#System of Equations
class finiteDiff_SystemOfEquations(object): 
    __metaclass__  = abc.ABCMeta

    @abc.abstractmethod
    def assembleSystemOfEquations(self,grid):
        pass

    @abc.abstractmethod
    def solveSystemOfEquations(self):
        pass

class finiteDiff_SystemOfEquations1dHeat(finiteDiff_SystemOfEquations): 

    def __init__(self,grid):
       print "System of Equations: 1d Heat"
       self.grid = grid
       self.ySolution = finiteDiff_1dSolution( grid.nDim, grid.nMeshPoints )

    def assembleSystemOfEquations(self,sources):
        print "Assemble System Of Equations"
        condition = True
        return condition

    def solveSystemOfEquations(self):
        print "Solve System Of Equations: %d"% len(self.ySolution.ss) 
        status = True
        return status

    def saveSolution(self, fileName):
        print "Save the solution in file %s"% fileName


#Main ===============
print "1D Heat Equation Solver"

xx = finiteDiff_1dGrid(nMeshPoints,problemInterval)

xx.setBoundaryProperties(xx.FIRST,FIXED)
xx.setBoundaryProperties(xx.LAST,FLUX)

hSource = finiteDiff_sources( heatSource )

AAbb = finiteDiff_SystemOfEquations1dHeat( xx )
fd_Condition = AAbb.assembleSystemOfEquations( hSource )
if fd_Condition:
    SOLVED = AAbb.solveSystemOfEquations()
    if SOLVED:
        AAbb.saveSolution( HeatEquationResultFile )
    else:
        print "ERROR: Solution Failed!"
else:
    print "ERROR: Assembly of system of equations Failed!"
