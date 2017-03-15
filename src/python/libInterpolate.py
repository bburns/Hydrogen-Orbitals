
"""
libInterpolate
Various interpolation and related routines.
"""


# import python modules
import string
import time
from scipy import *


#verbose = True
verbose = False




def FindZero(f, xmin, xmax, xchangemax, niterations, args=None):
    
    """
    Find the zero for the given function f that lies between xmin and xmax, 
    ie try to find value for x for which f(x)=0. 
    xchangemax - max change between iterations
    niterations - max number of iterations
    args - can be any arguments needed by the function f 
    """
    
    x1 = xmin
    x2 = xmax
    
    # Get first two values for f
    f1 = f(x1, args)
    f2 = f(x2, args)
    if (verbose):
        print 0,x1,f1
        print 1,x2,f2
    
    # Make sure they have opposite signs, since we want a zero between them
    assert(sign(f1) != sign(f2))
    
    # Try to find value for x that will give f=0
    # ie iterate until f converges to zero
    x3old = 0.0
    for i in range(niterations):
        # find new x and f(x)
        deltaX = x2-x1
        deltaF = f2-f1
        slopeInverse = deltaX / deltaF
        x3 = x1-f1*slopeInverse
        f3 = f(x3, args)
        if (verbose): print i+2,x3,f3
        # replace one of the endpoints with out new values
        if (sign(f3)==sign(f1)):
            x1=x3
            f1=f3
        else:
            x2=x3
            f2=f3
        # exit when change in x is < xchangemax
        if (abs(x3-x3old) < xchangemax): break
        x3old = x3
        
    x = x3
    
    # Make sure x is within originally specified range
    assert (x > xmin and x < xmax)

    return x
    
    

def Interpolate(ax, ay, x, npoints):
    
    """
    Interpolate or extrapolate from the given data to find the value
        for y at the point x.
    Arbitrary spacing in the x array is okay, but should be sorted in ascending order!
    Uses Lagrangian (polynomial) interpolation with <npoints> data points:
       y(x) = y1*l1(x) + y2*l2(x) + ... yn*ln(x)
       li(x) = (x-x1)/(xi-x1) * (x-x2)/(xi-x2) * ... (x-xn)/(xi-xn) [excluding ith term]
       y(x) will be a polynomial of order <npoints-1>
    """

    assert(ax[1]>ax[0]) # test for ascending order, at least for first point
    
    if (verbose): 
        print 'interpolate/extrapolate to x=',x,', npoints=',npoints

    # Find best data points to use, based on which are closest to 
    # requested point x. Will find <npoints> (or fewer) best data points and 
    # return as an array.
    ibest = FindBest(ax,x,npoints)
    npoints = len(ibest) # make sure npoints is updated in case was reduced
    if (verbose): 
        print 'ibest',ibest

    # Build the polynomial y(x), evaluated at the point x.
    y = 0.0
    for i in range(npoints): # do i=0,npoints-1
        li = 1.0
        ni = ibest[i] # index to ith best point
        # build up li[x] term, evaluated at the point x
        for j in range(npoints): # do j=0,npoints-1
            if (i != j): # exclude j=i term
                nj = ibest[j] # index to jth best point
                li = li*(x-ax[nj])/(ax[ni]-ax[nj])
        y = y+ay[ni]*li
    
    return y



def FindBest(ax, x, npoints):
    
    """
    Find <npoints> data points closest to requested point x and 
    return an array with indices to the points.
    """

    # get size of data array (goes 0 to ndata-1)
    ndata = len(ax)

    # Find index of point closest to x
    iclosest = BSearch(ax,x)
    if (verbose): 
        print 'looking for closest to x=',x
        print 'ax',ax
        print 'closest point at ',ax[iclosest],'iclosest=',iclosest
    
    # Get npoints points in each direction, and find distance 
    # from x for each. 
    # This will handle cases where point is at start or end of 
    # data set, or where all closest points lie in one direction.
    imin = iclosest-npoints
    imax = iclosest+npoints
    # make sure imin and imax are in array range
    if (imin < 0): imin = 0
    if (imax >= ndata): imax = ndata-1
    ncandidates = imax-imin+1
    if (verbose):
        print 'imin,imax,ncandidates',imin,imax,ncandidates
        print 'candidate points:'
        print '    j,i,xdata[i],xdelta[j]:'
    xdelta = zeros(ncandidates,'d') # initialize array
    for i in range(imin,imax+1): # do i=imin,imax
        j = i-imin
        xdelta[j] = abs(ax[i]-x) # distance from x
        if (verbose): print '    ',j,i,ax[i],xdelta[j]
    
    # Sort points by xdelta, in ascending order
    ibest = IndexSort(xdelta)
    
    # Exclude closest point if it's actually the point we're searching for
    # (dr mayes requirement)
    npoints2 = npoints
    if (xdelta[ibest[0]] == 0.0):
        if (verbose): print 'excluding point with xdelta=0'
        # reduce number of available candidates by one
        ncandidates -=1
        # make sure we don't have more points than candidates
        npoints2 = ncandidates 
        # shift candidates down by one
        for i in range(ncandidates):  # do i=0,ncandidates-1
            ibest[i]=ibest[i+1]
            
    # trim the array down to the number of requested or available points
#    ibest.resize(npoints)  # having trouble with this sometimes
    ibest = ibest[:npoints2]

    # adjust ibest array to correct range
    # note: at this point the first <npoints> is all we need
#    for i in range(npoints):   # do i=0,npoints-1
    for i in range(npoints2):   # do i=0,npoints-1
        ibest[i]=ibest[i]+imin
    
    if (verbose):
        print 'best points (sorted by xdelta):'
        print '    i,ibest,xdata,xdelta'
#        for i in range(npoints):  # do i=0,npoints-1
        for i in range(npoints2):  # do i=0,npoints-1
            print '    ',i,ibest[i],ax[ibest[i]],abs(x-ax[ibest[i]])
    
    return ibest



def IndexSort(xdata):
    
    """
    Sorts array xdata into ascending numerical order, and returns results as
    an array with the indexes.
    Uses bubble sort, so keep number of points low.
    """

    ndata=len(xdata)

    # initialize iorder array
#    iorder = arange(0,ndata)  # 0 to ndata-1
    iorder = zeros(ndata)
    for i in range(ndata): # do i=0,ndata-1
        iorder[i]=i # initialize 
    
    # bubble sort
    for i in range(ndata): #do i=0,ndata-1
        ni=iorder[i]
        xi=xdata[ni]
        for j in range(i+1,ndata): #do j=i+1,ndata-1
            nj=iorder[j]
            if (ni != nj):
                xj=xdata[nj]
                if (xi > xj):
                    iorder[i]=nj
                    iorder[j]=ni
                    ni=nj
                    xi=xj
                    # xj will get reassigned next round so don't worry about it here
                    
    return iorder


def BSearch(xdata,x):
    
    """
    Do binary search to find the index of the point closest to 
    the given point x.
    Array xdata should be in ascending order.
    """

    # check that first two data points are in ascending order, at least
    assert(xdata[1]>xdata[0])
    
    imin = 0
    imax = len(xdata)-1
    while True: 
        i = (imin+imax)/2
        if (x < xdata[i]):
            imax = i
        else:
            imin = i
        if (imin+1 >= imax): # done
            # use imin or imax, depending on which is closer to x
            if (abs(x-xdata[imin]) < abs(x-xdata[imax])):
                i = imin
            else:
                i = imax
            # done, so exit
            break

    return i



#---------------------------------------------------------------------------

if __name__ == '__main__':

    # do some tests on this module
    
    print
    print 'libInterpolate tests running...'
    print 'date: ',time.asctime()
    print
    
    # read data from file
    filename='hw2data.txt'   # file with bessel x,y data
    print 'reading data from file'
    file = open(filename,'r')
    ndata = int(file.readline())
    ndata += 1  # correct number of points
    print '    ndata points:',ndata
    # initialize arrays
    xdata = zeros(ndata,'d')
    ydata = zeros(ndata,'d')
    for i in range(ndata):  # do i=0,ndata-1
    #  read (11, 100) xdata[i],ydata[i]
        line = file.readline()
        strvalues = line.strip().split()
        xdata[i] = string.atof(strvalues[0])
        ydata[i] = string.atof(strvalues[1])
    file.close()
    
    # do some interpolations
    verbose = True
    npoints = 7 # number of points to use from data arrays
    print 
    print 'test 1:'
    x = 14.5
    y = Interpolate(xdata,ydata,x,npoints)
    yactual = special.j1(x)
    yerror = yactual-y
    print 'x','yestimate','yactual','yerror'
    print x,y,yactual,yerror
    print 
    print 'test 2:'
    x = 15.0
    y = Interpolate(xdata,ydata,x,npoints)
    yactual = special.j1(x)
    yerror = yactual-y
    print 'x','yestimate','yactual','yerror'
    print x,y,yactual,yerror
    print 
    print 'test 3:'
    x = 16.0
    y = Interpolate(xdata,ydata,x,npoints)
    yactual = special.j1(x)
    yerror = yactual-y
    print 'x','yestimate','yactual','yerror'
    print x,y,yactual,yerror
    verbose = False
    
    print 
    print 
    print 'testing error vs. number of points:'
    x = 14.0
    print 'x=',x
    print 'npoints','yestimate','yactual','yerror'
    for npoints in range(1,21): # do i=1,20
        y = Interpolate(xdata,ydata,x,npoints)
        yactual = special.j1(x)
        yerror = yactual-y
        print npoints,y,yactual,yerror
    
    print 
    print 
    print 'doing interpolations and extrapolations:'
    print 'x','yestimate','yactual','yerror'
    npoints = 7 
    xmin = 14.0
    xmax = 16.0
    nfind = 40 # number of points to find
    h = (xmax-xmin)/nfind
    for i in range(nfind+1): # do i=0,nfind
        x = xmin+h*i
        y = Interpolate(xdata,ydata,x,npoints)
        yactual = special.j1(x)
        yerror = yactual-y
        print x,y,yactual,yerror
    
