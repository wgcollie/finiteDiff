import numpy
from scipy.sparse import dia_matrix
from scipy.sparse.linalg import spsolve
import sys

# One Dimensional Steady State Heat Conduction

#Parameters
problemInterval = [0,1]
nMeshPoints = 10
heatSource = "Exponential"
HeatEquationResultFile = "HeatEquationSolutionFile.dat"

#Constants

class GridParameters:

    def __init__(self,nPoints,interval):
        print( "1d Grid" )
        self.nDim = 1
        self.interval = interval
        self.nPoints = nPoints
        self.meshSize = (interval[1]-interval[0])/nPoints
        print( "Grid Parameters:" )
        print( "  interval = [%f,%f]"%(interval[0],interval[1]) )
        print( "  nPoints = %d,   meshSize = %f"%(nPoints,self.meshSize) )

class BoundaryProperties:

    UNDEFINED = 0
    FIXED = 1
    FLUX = 2

    def __init__(self,nGridPoints):
        print( "Boundary Properties")
        self.FIRST = 0
        self.LAST = nGridPoints
        self.bdryProperties = [(0,0.0),(0,0.0)]
        print( self.FIRST, self.LAST )

    def setFirstBoundaryProperties(self,bdType,bdValue):
        print( "Set Boundary Properties", end=" " )
        if bdType == 'FIXED':
            bt = BoundaryProperties.FIXED
        elif bdType == 'FLUX':
            bt = BoundaryProperties.FLUX
        else:
            print( 'Invalid Boundary Type %s for first boundary point'%bdType )
        self.bdryProperties[0] = (bt,bdValue)
        print( "Type: %s, Value: %f"% (bdType, bdValue) )

    def setLastBoundaryProperties(self,bdType,bdValue):
        print( "Set Boundary Properties", end=" " )
        if bdType == 'FIXED':
            bt = BoundaryProperties.FIXED
        elif bdType == 'FLUX':
            bt = BoundaryProperties.FLUX
        else:
            print( 'Invalid Boundary Type %s for last boundary point'%bdType )
        self.bdryProperties[1] = (bt,bdValue)
        print( "Type: %s, Value: %f"% (bdType, bdValue) )

    def getFirstBdryType(self):
        return self.bdryProperties[0][0]

    def getLastBdryType(self):
        return self.bdryProperties[1][0]

    def getFirstBdryValue(self):
        return self.bdryProperties[0][1]

    def getLastBdryValue(self):
        return self.bdryProperties[1][1]

#Grid
class GridPoint1d:

    baseValue = 0.0
    delta = 0.0

    def setCoordinateValue1d(ii):
        x = GridPoint1d.baseValue + ii*GridPoint1d.delta
        return x
 
    def __init__(self,ii):
        self.ii = ii
        self.value = GridPoint1d.setCoordinateValue1d(ii)

    def x(self):
        return self.value

    def idx(self):
        return self.ii

class FiniteDiff_Grid:

    def __init__(self,gridParams):
        self.GridParams = gridParams
        GridPoint1d.baseValue = gridParams.interval[0]
        GridPoint1d.delta = gridParams.meshSize

    def __iter__(self):
        self.ii = 0
        return self
    
    def __next__(self):
        idx = self.ii
        if self.ii > self.GridParams.nPoints:
            raise StopIteration
        self.ii += 1
        return GridPoint1d(idx)

class finiteDiff_sources:

    def __init__(self):
        self.sourceList = []

    def addSource(self,src):
        self.sourceList.append( src )

    def evaluateSources(self,x):
        tt = 0.0
        for src in self.sourceList:
            tt += src.evaluateSourceContribution(x)
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
        self.ss = numpy.zeros( nGridPoints )

#System of Equations
class finiteDiff_SystemOfEquations(object):

    def __init__(self,grid,appSpecificEqnMethds):
        self.grid = grid
        self.appEqnMthds = appSpecificEqnMethds

    def assembleSystemOfEquations(self,sources):
        print( "Assemble System Of Equations" )
        for p in self.grid:
            self.appEqnMthds.equationContribution(sources,p)
        self.matrix,self.vctr = self.appEqnMthds.assembleEqns()
        return True
            
    def solveSystemOfEquations(self):
        print( "Solve System Of Equations" )
        self.solution = spsolve(self.matrix.tocsr(),self.vctr)
        return True

    def saveSolution(self, fileName):
        print( "Save the solution in file %s"% fileName )
        fp = open( fileName, 'w' )
        for num in self.solution:
            fp.write('%f\n'%num)

class finiteDiff_SystemOfEquations1dHeatMthds: 

    def __init__(self,gridParams,bdryProperties):
       print( "System of Equations: 1d Heat" )
       self.gridParams = gridParams
       self.bdryProperties = bdryProperties
       neq = self.gridParams.nPoints+1
       self.dd1 = (neq)*[0.0]
       self.dd2 = (neq)*[0.0]
       self.dd3 = (neq)*[0.0]
       self.vctr = (neq)*[0.0]
       print( [self.dd1, self.dd2, self.dd3] )

    def equationContribution(self,sources,gridPoint):
        ii = gridPoint.idx()
        pp = gridPoint.x()
        h = self.gridParams.meshSize
        hSquared = h*h
        print( "equation contribution", ii, pp, end=" " )
        if ii == self.bdryProperties.FIRST:
            if self.bdryProperties.getFirstBdryType() == BoundaryProperties.FIXED:
                print( "First Fixed" )
                self.dd2[ii] = 1.0
                self.dd3[ii+1] = 0.0
                self.vctr[ii] = self.bdryProperties.getFirstBdryValue()
            elif self.bdryProperties.getFirstBdryType() == BoundaryProperties.FLUX:
                print( "First Flux" )
                self.dd2[ii] = -2.0
                self.dd3[ii+1] = 2.0
                self.vctr[gridPoint.ii] = -2*h - 2*hSquared*sources.evaluateSources(pp)
        elif ii == self.bdryProperties.LAST:
            if self.bdryProperties.getLastBdryType() == BoundaryProperties.FIXED:
                print( "Last Fixed" )
                self.dd2[ii] = 1.0
                self.vctr[ii] = self.bdryProperties.getLastBdryValue()
            elif self.bdryProperties.getLastBdryType() == BoundaryProperties.FLUX:
                print( "Last Flux" )
                self.dd1[ii-1] = 2.0
                self.dd2[ii] = -2.0
                self.vctr[ii] = -2*h - 2*hSquared*sources.evaluateSources(pp)
        else:
            self.dd1[ii-1] =  1.0
            self.dd2[ii] = -2.0
            self.dd3[ii+1] =  1.0
            self.vctr[ii] = -hSquared*sources.evaluateSources(pp)
        print( self.dd1[ii], self.dd2[ii], self.dd3[ii] )

    def assembleEqns(self):
        print( "Assemble Equation Components" )
        print( [self.dd1, self.dd2, self.dd3] )
        neq = self.gridParams.nPoints+1
        data = numpy.array( [self.dd1, self.dd2, self.dd3] )
        offsets = numpy.array( [-1,0,1])
        self.matrix = dia_matrix( (data,offsets), shape=(neq,neq))
        status = True
        print( self.matrix.todense() )
        print( self.vctr )
        return (self.matrix,self.vctr)

#Main ===============
print( "1D Heat Equation Solver" )

grdPrms = GridParameters(10,[0.0,1.0])
xx = FiniteDiff_Grid(grdPrms)

bdryPrp = BoundaryProperties(grdPrms.nPoints)
bdryPrp.setFirstBoundaryProperties('FIXED',0.0)
bdryPrp.setLastBoundaryProperties('FLUX',1.0)
print( "First boundary node: %d, Last boundary node: %d" \
       %(bdryPrp.FIRST,bdryPrp.LAST ) )

deqSources = finiteDiff_sources()
deqSources.addSource( finiteDiff_sourceElement( 2.0 ) )

soeMethods = finiteDiff_SystemOfEquations1dHeatMthds(grdPrms,bdryPrp)
AAbb = finiteDiff_SystemOfEquations( xx,soeMethods )
fd_Status = AAbb.assembleSystemOfEquations( deqSources )
print( fd_Status )
if fd_Status:
    SOLVED = AAbb.solveSystemOfEquations()
    if SOLVED:
        AAbb.saveSolution( HeatEquationResultFile )
    else:
        print( "ERROR: Solution Failed!" )
else:
    print( "ERROR: Assembly of system of equations Failed!" )
