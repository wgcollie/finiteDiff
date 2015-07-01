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

    def __init__(self,nGridPoints):
        print( "Boundary Properties")
        self.FIRST = 0
        self.LAST = nGridPoints
        self.bdryProperties = [(0,0.0),(0,0.0)]

    def setFirstBoundaryProperties(self,bdType,bdValue):
        print( "Set Boundary Properties", end=" " )
        if bdType == 'FIXED':
            bt = BoundaryProperties.FIXED
        elif bdTypd == 'FLUX':
            bt = BoundaryProperties.FLUX
        else:
            print( 'Invalid Boundary Type %s for first boundary point'%bdType )
        self.bdryProperties[0] = (bdType,bdValue)
        print( "Type: %s, Value: %f"% (bdType, bdValue) )

    def setLastBoundaryProperties(self,bdType,bdValue):
        print( "Set Boundary Properties", end=" " )
        if bdType == 'FIXED':
            bt = BoundaryProperties.FIXED
        elif bdType == 'FLUX':
            bt = BoundaryProperties.FLUX
        else:
            print( 'Invalid Boundary Type %s for last boundary point'%bdType )
        self.bdryProperties[1] = (bdType,bdValue)
        print( "Type: %s, Value: %f"% (bdType, bdValue) )

    def getFirstBdryPrprty(self):
        return self.bdryProperties[0][0]

    def getLastBdryPrprty(self):
        return self.bdryProperties[1][0]

    def getFirstBdryValue(self):
        return self.bdryProperties[0][1]

    def getLastBdryValue(self):
        return self.bdryProperties[1][1]

#Grid
class GridPoint1d:

    baseValue = 0
    delta = 0

    def setCoordinateValue1d(ii):
        x = gridPoint1d.baseValue + ii*gridPoint1d.delta
        return x
 
    def __init__(self,ii):
        self.ii = ii
        self.value = setCoordinateValue1d(ii)

    def x(self):
        return self.value

    def idx(self):
        return self.ii

class FiniteDiff_Grid:

    def __init__(self,gridParams):
        self.GridParams = gridParams

    def __iter__(self):
        self.ii = 0
        return self
    
    def __next__(self):
        idx = self.ii
        if self.ii > self.GridParams.nMeshPoints:
            raise StopIteration
        self.ii += 1
        return GridPoint1d(idx,self.GridParams)

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

    def __init__(self,grid,appSpecificEqnMethds):
        self.grid = grid
        self.appEqnMthds = appSpecificEqnMethds

    def assembleSystemOfEquations(self,sources):
        print( "Assemble System Of Equations" )
        for p in self.grid:
            self.appEqnMthds.equationContribution(sources,p)
        self.appEqnMthds.assembleEqns()
            
    def solveSystemOfEquations(self):
        pass

class finiteDiff_SystemOfEquations1dHeatMthds: 

    def __init__(self,bdryProperties):
       print( "System of Equations: 1d Heat" )
       self.bdryProperties = bdryProperties
       self.dd1 = []
       self.dd2 = []
       self.dd3 = []
       self.vector = []

    def equationContribution(self,sources,gridPoint):
        print( "equation contribution" )
        if gridPoint.bdryIdx() == BoundaryProperties.FIRST:
            if getFirstBdryType() == BoundaryProperties.FIXED:
                self.dd2[gridPointp.ii] = 1.0
                self.dd3[gridPointp.ii] = 0.0
                self.vector[gridPointp.ii] = self.bdryProperties.getFirstBdryValue()
            elif getFirstBdryType() == BoundaryProperties.FLUX:
                self.dd2[gridPointp.ii] = -1.0
                self.dd3[gridPointp.ii] = 1.0
                self.vector[gridPointp.ii] = sources.evaluateSources(p.x())
        elif gridPointp.bdryIdx() == BoundaryProperties.LAST:
            if getLastBdryType() == BoundaryProperties.FIXED:
                self.dd1[gridPointp.ii] = 0.0
                self.dd2[gridPointp.ii] = 1.0
                self.vector[gridPointp.ii] = self.bdryProperties.getLastBdryValue()
            elif getLastBdryType() == BoundaryProperties.FLUX:
                self.dd1[gridPointp.ii] = -1.0
                self.dd2[gridPointp.ii] = 1.0
                self.vector[gridPointp.ii] = sources.evaluateSources(p.x())
        else:
            self.dd1[gridPointp.ii] =  1.0
            self.dd2[gridPointp.ii] = -2.0
            self.dd3[gridPointp.ii] =  1.0
            self.vector[gridPointp.ii] = sources.evaluateSources(p.x())  

    def assembleEqns(self):
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

gridParams = GridParameters(10,[0.0,1.0])
xx = FiniteDiff_Grid(gridParams)

bdryPrp = BoundaryProperties(gridParams.nMeshPoints)
bdryPrp.setFirstBoundaryProperties('FIXED',0.0)
bdryPrp.setLastBoundaryProperties('FLUX',1.0)

deqSources = finiteDiff_sources()
deqSources.addSource( finiteDiff_sourceElement( 2.0 ) )

soeMethods = finiteDiff_SystemOfEquations1dHeatMthds(bdryPrp)
AAbb = finiteDiff_SystemOfEquations( xx,soeMethods )
fd_Staus = AAbb.assembleSystemOfEquations( deqSources )
if fd_Condition:
    SOLVED = AAbb.solveSystemOfEquations()
    if SOLVED:
        AAbb.saveSolution( HeatEquationResultFile )
    else:
        print( "ERROR: Solution Failed!" )
else:
    print( "ERROR: Assembly of system of equations Failed!" )
