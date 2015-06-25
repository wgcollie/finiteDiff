import abc
from numpy import *
from numpy.linalg import *

# One Dimensional Steady State Heat Conduction

#Parameters
problemInterval = [0,1]
nMeshPoints = 10
heatSource = "Exponential"
HeatEquationResultFile = "HeatEquationSolutionFile.dat"

#Constants
FIXED = 0
FLUX = 1

class GridParameters:

    def __init__(self,nMeshPoints,interval):
        print( "1d Grid" )
        self.nDim = 1
        self.interval = interval
        self.nMeshPoints = nMeshPoints
        self.meshSize = (interval[1]-interval[0])/nMeshPoints

class BoundaryProperties:

    FIRST = 0
    LAST = 1
    UNDEFINED = 0
    FIXED = 1
    FLUX = 2

    def __init__(self,lowerBdryPrprty,upperBdryPrprty):
        print( "Boundary Properties")
        self.bdryProperties = [(0,0.0),(0,0.0)]

#Grid
class gridPoint1d:

    def setStencils1d(ii,m_max):    
        # l:lower (_ X X)    u:upper (X X _)
       if ii == 0:
           self.stencilType = lowerLimStencil
       elif ii = m_max:
           self.stencilType = upperLimStencil
       else:
           self.stencilType = interiorStencil

    def setCoordinateValue1d(ii,baseValue,delta):
        x = baseValue + ii*delta
       return x
 
    def __init__(self,ii,gridParameters):
        m_max = 
        baseValue = 
        delta = 
        self.stencilType = gridPoint1d.setStencils1d(ii)       
        self.value = setCoordinateValue1d(ii,baseValue,delta)

    def x(self):
        return self.val

    def bdryValue(self):
        return self.bdryVal

    def bdryTypeCode(self):
        return self.bdryType


class FiniteDiff_Grid:

    def __init__(self,gridParams):
        self.GridParams = GridParams

    def __iter__(self):
        self.ii = 0
       return self
    
    def __next__(self):
       idx = self.ii
       if self.ii > Grid.m_max:
            raise StopIteration
        self.ii += 1
        return GridPoint(idx)


class finiteDiff_1dGrid:

    def __init__(self,nMeshPoints,interval):
        print( "1d Grid" )
        self.nDim = 1
        self.interval = interval
        self.nMeshPoints = nMeshPoints
        self.meshSize = (interval[1]-interval[0])/nMeshPoints
        self.bdryProperties = [(0,0.0),(0,0.0)]

    def setFirstBoundaryProperties(self,bdType,bdValue):
        print( "Set Boundary Properties", end=" " )
        self.bdryProperties[0] = (bdType,bdValue)
        print( "Type: %s, Value: %f"% (bdType, bdValue) )

    def setLastBoundaryProperties(self,bdType,bdValue):
        print( "Set Boundary Properties", end=" " )
        self.bdryProperties[1] = (bdType,bdValue)
        print( "Type: %s, Value: %f"% (bdType, bdValue) )


    def nGridPoints(self):
        return self.nMeshPoints

class FiniteDiff_Grid:

    def __init__(self,gridParams):
        self.GridParams = GridParams

    def __iter__(self):
        self.ii = 0
       return self
    
    def __next__(self):
       idx = self.ii
       if self.ii > Grid.m_max:
            raise StopIteration
        self.ii += 1
        return GridPoint(idx)


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
            print( "ERROR: %d is FORBIDDEN for alpha"% alpha )
            sys.exit()

    def evaluateSourceContribution(self,x):
        #print "HeatSource Contribution: %s"% self.source
        tt = -(self.alpha+1)*( x**self.alpha )
        return tt

#Solution
class finiteDiff_1dSolution:

    def __init__(self, nGridPoints):
        self.nDim = 1
        self.nGridPoints = nGridPoints
        self.ss = zeros( nGridPoints )

#System of Equations
class finiteDiff_SystemOfEquations(object):

    def __init__(self,appSpecificEqnMethds):
        self.appEqnMthds

    def assembleSystemOfEquations(self,grid,sources):
        print( "Assemble System Of Equations" )
        for p in self.grid:
            self.appEqnMthds.equationContribution(sources,p)
            
    def solveSystemOfEquations(self):
        pass

class finiteDiff_SystemOfEquations1dHeatMthds: 

    def __init__(self,grid):
       print( "System of Equations: 1d Heat" )
       self.grid = grid
       self.dd1 = []
       self.dd2 = []
       self.dd3 = []
       self.vector = []

    def equationContribution(self,sources,gridPoint):
        print( "equation contribution" )
        if gridPoint.bdryIdx() == BoundaryProperties.FIRST:
            if gridPointp.bdryTypeCode() == gridPoint.FIXED:
                self.dd2[gridPointp.ii] = 1.0
                self.dd3[gridPointp.ii] = 0.0
                self.vector[gridPointp.ii] = gridPointp.bdryValue()
            elif gridPointp.bdryTypeCode() == gridPoint.FLUX:
                self.dd2[gridPointp.ii] = -1.0
                self.dd3[gridPointp.ii] = 1.0
                self.vector[gridPointp.ii] = sources.evaluateSources(p.x())
        elif gridPointp.bdryIdx() == BoundaryProperties.LAST:
            if gridPointp.bdryTypeCode() == gridPoint.FIXED:
                self.dd1[gridPointp.ii] = 0.0
                self.dd2[gridPointp.ii] = 1.0
                self.vector[gridPointp.ii] = gridPointp.bdryValue()
            elif gridPointp.bdryTypeCode() == gridPoint.FLUX:
                self.dd1[gridPointp.ii] = -1.0
                self.dd2[gridPointp.ii] = 2.0
                self.vector[gridPointp.ii] = sources.evaluateSources(p.x())
        else:
            self.dd1[gridPointp.ii] =  1.0
            self.dd2[gridPointp.ii] = -2.0
            self.dd3[gridPointp.ii] =  1.0
            self.vector[gridPointp.ii] = sources.evaluateSources(p.x())
  

    def assembleEqns(self)
        data = numpy.array( [self.dd1, self.dd2, self.dd3] )
        offsets = numpy.array( [-1,0,1])
        self.matrix = dia_matrix( (data,offsets), shape=(self.nn,self.nn))
        status = True
        return status

    def solveSystemOfEquations(self):
        print( "Solve System Of Equations: %d"% len(self.ySolution.ss) )
        status = True
        return status

    def saveSolution(self, fileName):
        print( "Save the solution in file %s"% fileName )


#Main ===============
print( "1D Heat Equation Solver" )

xx = finiteDiff_1dGrid(nMeshPoints,problemInterval)

xx.setFirstBoundaryProperties('FIXED',0.0)
xx.setLastBoundaryProperties('FLUX',1.0)

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
        print( "ERROR: Solution Failed!" )
else:
    print( "ERROR: Assembly of system of equations Failed!" )
