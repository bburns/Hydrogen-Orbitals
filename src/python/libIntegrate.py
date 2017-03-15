

"""
libIntegrate
Various integration routines:
    Integrate1stOrder
    Solve2ndOrder
    Integrate2ndOrder
    IntegrateDiscrete
"""


# Import modules
import time
from scipy import *
from libGeneral import *


#verbose = True
verbose = False 



def Integrate1stOrder(f, xmin, xmax, y0=0, nsteps=10000, args=None):
    
    """
    Integrate the function f from xmin to xmax and return the result.
    Uses nsteps for number of steps,
    The function f should return y' at x,y. 
    y0 is the initial y value, usually zero.
    Uses 5th order Runge-Kutta.
    args is an arbitrary argument that will be passed to the function f.
    """
    
    h = (xmax - xmin) / nsteps    # step distance
    
    # initial point
    x = xmin
    y = y0
    
    for i in xrange(nsteps):
        k0 = h * f(x, y, args)
        k1 = h * f(x + h / 2.0, y + k0 / 2.0, args)
        k2 = h * f(x + h, y + k1, args)
        k3 = h * f(x + h, y + k2, args)
        y = y + (k0 + 2.0 * k1 + 2.0 * k2 + k3) / 6.0
        x = xmin + i * h
    return y




def Solve2ndOrder(g,xmin,xmax,y0=0,yprime0=0,nsteps=10000,
                            ndata=100,args=None): 
    """
    Integrate the function g from xmin to xmax and return arrays for x and y.
    Uses nsteps for number of steps, ndata+1 for number of data points in arrays. 
    y0 and yprime0 are initial values. 
    The function g should return y'' for x,y,y'
    args is an arbitrary argument that will be passed to the function g. 
    Uses 5th order Runge-Kutta.
    """
    
    if (verbose): print 'Solve2ndOrder(',xmin,xmax,y0,yprime0,nsteps,ndata,args,')...'
    
    # get step size
    h = (xmax - xmin) / nsteps
    
    # get number of points to skip between writing to arrays
    nskip = nsteps/ndata    
    assert(nsteps==nskip*ndata)  # nsteps must be evenly divisible by ndata
    
    # create arrays
    ax = zeros(ndata+1,'d')
    ay = zeros(ndata+1,'d')

    # initial point
    x = xmin
    y = y0
    yprime = yprime0
    
    if (verbose):
        print 'nsteps=',nsteps
        print 'ndata=',ndata
        print 'nskip=',nskip
        print 'h(step distance)=',h
        print 'args=',args
        print '    j,x,y:'
    
    # save first point to arrays
    j=0
    ax[j] = x
    ay[j] = y
    if (verbose): print '    ',j,x,y
    
    # loop over number of steps
    for i in xrange(1,nsteps+1):  # do i=1,nsteps
        
        k0 = h * yprime
        m0 = h * g(x, y, yprime, args)
        k1 = h * (yprime+m0/2)
        m1 = h * g(x+h/2, y+h*yprime/2, yprime+m0/2, args)
        k2 = h * (yprime+m1/2)
        m2 = h * g(x+h/2, y+h*yprime/2+h*m0/4, yprime+m1/2, args)
        k3 = h * (yprime+m2)
        m3 = h * g(x+h, y+h*yprime+h*m1/2, yprime+m2, args)
        
        # order of y and yprime calculation is important here!
        x = xmin + h*i
        y = y + h*yprime+h*(m0+m1+m2)/6
        yprime = yprime + (m0+2*m1+2*m2+m3)/6
        
        # save x,y to arrays every <nskip> data points, including last data point
        if (i % nskip == 0):
            j=j+1
            ax[j]=x
            ay[j]=y
            if (verbose): print '    ',j,x,y
    
    # done...
    return ax, ay
    


def Integrate2ndOrder(g,xmin,xmax,y0=0,yprime0=0,nsteps=10000,args=None):
    """
    Integrate a 2nd order equation by using our ode solver, and return the 
    last value for y.
    """
    ndata = 1 # save first and last points to array
    ax,ay = Solve2ndOrder(g,xmin,xmax,y0,yprime0,nsteps,ndata,args)
    y = ay[-1]  # get last value for y
    return y
    


def IntegrateDiscrete(ax,ay):
    """
    Find the area under the curve specified by the given x and y arrays.
    Equal spacing in x array is assumed. 
    """
    area = 0.0
    ndata = len(ax)
    h = abs(ax[1]-ax[0]) # assume equal spacing
    for i in range(ndata-1):  # do i=0,ndata-2
        a = h*0.5*(ay[i]+ay[i+1]) # simplest version
        area = area + a
    return area
    


#-------------------------------------------------------------------------------

    
def ftest(x, y, args):
    "Sample f function - calculates the derivative at some point x,y."
    # Handle bad behavior at x=0 and x=1
    if ((x == 0.0) or (x == 1.0)):
        yprime = 0.0
    else:
        yprime = math.log(x) * math.log(1.0 - x) # log = natural log!
    return yprime


def gtest(x, y, yprime, args):
    "Sample g function - calcs y'' at some x,y,y'."
    # eg Bessel equation: y'' = -y'/x - (x**2 - n**2) * y/x**2
    n = 1 # order of bessel function
    # Handle bad behavior at x=0
    if (x == 0.0):
        g = 0.0
    else:
        g = - yprime/x - (x**2-n**2)*y/x**2
    return g


if __name__ == "__main__":

    # test this module
    
    print 
    print 'module library - running tests...'
    print 'date:',now()
    print 

    verbose = True    # global variable, accessible to all functions
    passedTests = True
    
    print "test1:"

    # range to integrate over
    xmin = 0.0
    xmax = 1.0
    
    # initial y value
    y0 = 0.0
    
    # actual value of integral for comparison
    yexact = 2.0 - (pi ** 2.0) / 6.0
    
    print 'nsteps           error             time (sec)'

    # try different values of interval h - looking for minimum error
    for i in range(4):
        nsteps = 10**(i+1)
        t = time.time()
        y = Integrate1stOrder(ftest, xmin, xmax, y0, nsteps)
        t = time.time() - t
        yerror = yexact - y
        print nsteps,'    ',yerror,'     ',t
    print

    if (abs(yerror) > 1e-7): passedTests=False

    print "test2:"
    
    # range to integrate over
    xmin = 0.0
    xmax = 15.0
    
    # initial values (at xmin), calculated from power series for Bessel's equation
    y0 = 0.0
    yprime0 = 0.5
    
    # number of data points to output to file
    ndata = 10
    
    # actual value of integral for comparison
    # Bessel function of order 1
    yexact = special.j1(xmax) - special.j1(xmin) 
    
    # do integration
#    nsteps = 10**5
    nsteps = 10**4
    t = time.time()
    y = Integrate2ndOrder(gtest, xmin, xmax, y0, yprime0, nsteps)
    t = time.time() - t
    yerror = yexact - y
    
    print 'y=',y
    print 'yexact=',yexact
    print 'yerror=',yerror
    print 'time=',t,' sec'
    if (abs(yerror) > 1e-7): passedTests=False
    print
    
    
    print 'test3:'
    verbose = False
    ax,ay = Solve2ndOrder(gtest,xmin,xmax,y0,yprime0,nsteps)
    aj = special.j1(ax)  # vectorized function!
    for i in range(len(ax)):
        print ax[i],ay[i],aj[i],aj[i]-ay[i]
#    gplt.plot(ax,ay,ax,aj)
    
    
    
    if (not passedTests): print 'WARNING: tests not passed!!'
    print 'tests finished.'
        


    
