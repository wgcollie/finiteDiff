import abc
from numpy import *
from scipy import *
from scipy.sparse import *

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

    FIRST = 0
    LAST = 1
    UNDEFINED = 0
    FIXED = 1
    FLUX = 2


    def __init__(self,nMeshPoints,interval):
        print "1d Grid"
        self.nDim = 1
        self.interval = interval
        self.nMeshPoints = nMeshPoints
        self.meshSize = (interval[1]-interval[0])/nMeshPoints
        self.bdryProperties = {finiteDiff_1dGrid.FIRST:(finiteDiff_1dGrid.UNDEFINED,0.0), \
                               finiteDiff_1dGrid.LAST:(finiteDiff_1dGrid.UNDEFINED,0.0)}

    def setFirstBoundaryProperties(self,bdType,bdValue):
        print "Set Boundary Properties"
        self.bdryProperties[ finiteDiff_1dGrid.FIRST ] = (bdType,bdValue)

    def setLastBoundaryProperties(self,bdType,bdValue):
        print "Set Boundary Properties"
        self.bdryProperties[ finiteDiff_1dGrid.LAST ] = (bdType,bdValue)

    def nGridPoints(self):
        return self.nMeshPoints

    def isBdryTypeFixed( pointID ):
        if self.bdryProperties[ pointID ][0] == finiteDiff_1dGrid.FIXED:
            rval = True
        else:
            rval = False
        return rval

    def isBdryTypeFlux( pointID ):
        if self.bdryProperties[ pointID ][0] == finiteDiff_1dGrid.FLUX:
            rval = True
        else:
            rval = False
        return rval

    def __getitem__(self,k):
        x = self.interval[0] + k*self.meshSize

class finiteDiff_sources:

    def __init__(self):
        self.sourceList = []

    def addSource(self,source):
        self.sourceList.append( source )

    def evaluateSources(self,x):
        tt = 0.0
        for source in self.sourceList:
            tt += source.evaluateSourceContribution(x)
        return tt
            
class finiteDiff_sourceElement:

    def __init__(self,alpha):
        self.alpha = alpha
        if alpha == -1 or alpha == -2:
            print "ERROR: %d is FORBIDDEN for alpha"% alpha
            sys.exit()

    def evaluateSourceContribution(self,x):
        #print "HeatSource Contribution: %s"% self.source
        tt = -(alpha+1)*( x**alpha )
        return tt

#Solution
class finiteDiff_1dSolution:

    def __init__(self, nGridPoints):
        self.nDim = 1
        self.nGridPoints = nGridPoints
        self.ss = zeros( nGridPoints )

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
       self.ySolution = finiteDiff_1dSolution( grid.nMeshPoints )

    def assembleSystemOfEquations(self,sources):
        print "Assemble System Of Equations"
        nn = self.grid.nGridPoints()
        ii = 0
        dd1 = zeros( nn )
        dd2 = zeros( nn )
        dd3 = zeros( nn )
        for p in self.grid:
            dd1[ii] =  1.0
            dd2[ii] = -2.0
            dd3[ii] =  1.0
            if p.isBdryTypeFixed:
                self.vector[ii] = p.boundaryProperties.Value
            elif p.isBdryTypeFlux:
                self.vector[ii] = p.boundaryProperties.Value
                self.vector[ii] = sources.evaluateSourceContribution(xx)
            else:
                xx = p.gridPoint()
                self.vector[ii] = sources.evaluateSourceContribution(xx)
            ii += 1
        data = numpy.array( [dd1, dd2, dd3] )
        offsets = numpy.array( [-1,0,1])
        self.matrix = dia_matrix( (data,offsets), shape=(nn,nn))
        status = True
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

xx.setFirstBoundaryProperties(finiteDiff_1dGrid.FIXED,0.0)
xx.setLastBoundaryProperties(finiteDiff_1dGrid.FLUX,1.0)

heatSource = finiteDiff_sourceElement( 2.0 )
deqSources = finiteDiff_sources()
deqSources.addSource( heatSource )

AAbb = finiteDiff_SystemOfEquations1dHeat( xx )
fd_Staus = AAbb.assembleSystemOfEquations( deqSources )
if fd_Condition:
    SOLVED = AAbb.solveSystemOfEquations()
    if SOLVED:
        AAbb.saveSolution( HeatEquationResultFile )
    else:
        print "ERROR: Solution Failed!"
else:
    print "ERROR: Assembly of system of equations Failed!"
