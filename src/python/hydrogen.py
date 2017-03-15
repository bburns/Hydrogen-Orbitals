#----------------------------------------------------------------------------
# Homework 5: Hydrogen Atom
# Brian Burns
# PHYS 4350 Computational Physics
# Dr. Mayes
# May 7, 2004
#----------------------------------------------------------------------------
# Find the eigenenergy and eigenfunction for the first excited state of Hydrogen. 
# Some variables: 
#    rho = distance from nucleus in Bohr radii = r/a0
#    epsilon = energy of electron in Rydbergs
#    Psi(r,angles) = R(r) * Y(angles)
#    Use change of variable u(r) = R(r) * r
# We'll solve for u(rho) to get R(rho), then plot R(rho) vs rho.
#----------------------------------------------------------------------------


# Import Python modules
from scipy import * # also includes arrays, which can be vectorized


# Import our routines into global namespace
from libConstants import *
from libGeneral import *
from libIntegrate import *
from libInterpolate import *


gifTerminal = 'gif large size 1024,760'
psTerminal = 'postscript enhanced color'


def upp(rho, u, uprime, args):  # ie g(x,y,yprime,args)

    """ 
    This is effectively the Shroedinger equation.
    Calculates u'' at some point rho, u, u'.
    From u'' = -(epsilon + 2z/rho -L(L+1)/rho^2)u, with z=1 for Hydrogen.
    Note: R(rho) = u(rho)/rho.
    args should contain epsilon and L
    """

    z = 1 # atomic number
    
    # Handle bad behavior at rho=0
    if (rho==0.0):
        upp = choose(args.L, (-2.0*z, 2.0, 1.0))
    else:
        upp = -(args.epsilon + 2.0*z/rho - args.L*(args.L + 1)/rho**2)*u
    
    return upp


def GetEnergy(epsilonMin,epsilonMax,L,rho0,nsteps):
    
    """ 
    Find the energy (eigenenergy) for which u(0)=0, by integrating 
    our Schroedinger equation from the starting point rho0 to the origin.
    If our energy is correct the SE should integrate to zero at the origin,
    since that is one of the boundary conditions. The other BC is at the 
    starting point rho0, for which we have approximate values for u and u'.
    The search will begin between the specified values for the energy. 
    nsteps is used in the integrations.
    """

    if (verbose): print 'GetEnergy(',epsilonMin,epsilonMax,L,rho0,nsteps,')...'

    niterations = 15
    xchangemax = 1e-12
    
    # Pack arguments needed by f into an empty object
    args = EmptyObject() 
    args.L = L
    args.rho0 = rho0
    args.nsteps = nsteps
    
    # Define a callback function for use by FindZero
    def f(epsilon,args): 
        return IntegrateSchroedingerToZero(epsilon,args.L,args.rho0,args.nsteps)

    # Call FindZero to find epsilon
    epsilon = FindZero(f,epsilonMin,epsilonMax,xchangemax,niterations,args)

    return epsilon



def IntegrateSchroedingerToZero(epsilon, L, rho0, nsteps):
    
    """ 
    Get value for u at 0 by integrating the Schroedinger equation
    from rho0 to 0 (going backwards). 
    Note: R(rho) = u(rho)/rho.
    """
    
    if (verbose): print 'IntegrateSchroedingerToZero(',epsilon,L,rho0,nsteps,')...'
        
    # Must have enough points here for discrete integral to work for normalization,
    # and nsteps must be evenly divisible by ndata to work properly!
    ndata = 200
    rhostart = rho0
    rhostop = 0.0
    aRho, au, A = SolveSchroedinger(epsilon,L,rhostart,rhostop,nsteps,ndata)
    u = au[-1]  # get last value, which is value of integral, ie value of u at rho=0
    if (verbose):
        print '    u(0)=',u
    return u




def SolveSchroedinger(epsilon,L,rhostart,rhostop,nsteps,ndata=1):
    
    """
    Get the normalized eigenfunction u as a function of position, over the 
    range rho=rhostart to rhostop. Returns arrays aRho and au, and 
    normalization factor A. 
    Note: u=R*r
    ndata determines size of arrays produced.
    """

    if (verbose): print 'SolveSchroedinger(',epsilon,L,rhostart,rhostop,nsteps,ndata,')...'

    # Calculate the start values
    if (rhostart==0.0):
        # Start values for rho=0
        u0 = 0.0
        u0prime = choose(L,(1.0,0.0,0.0))
    else:
        # Start values for large rho
        # These come from the approximation for u at large rho,
        #     u=B*exp(-s*rho)
        # with B a scaling factor that will be accounted for when u is normalized.
        # Note: +/- u0's are equivalent, since only probability density is physically 
        # measurable, but the negative one matches the plot in Griffiths.
        s = sqrt(-epsilon)
        u0 = -exp(-s*rhostart)
        u0prime = -s * u0
        
    if (verbose): 
        print '   u0',u0
        print '   u0prime',u0prime

    # Pack epsilon and L into an argument object that will be passed to 
    # Solve2ndOrder and from there to our function upp. 
    args = EmptyObject()
    args.epsilon = epsilon
    args.L = L

    # Do 2nd order integration from rhostart to rhostop
    #                                 (g, xmin, xmax, y0, yprime0, nsteps, ndata, args)
    aRho, au = Solve2ndOrder(upp, rhostart,rhostop,u0,u0prime, nsteps, ndata, args)

    # Now normalize our u array
    A = GetA(epsilon, aRho, au)
    au = au/sqrt(A)

    return aRho, au, A



def GetA(epsilon, aRho, au):
    
    """
    Calculate the integral of u**2 from rho=0 to infinity.
    This is used in normalizing u (and hence Psi).
    Calculates integral from rho=0 to rho0 numerically using the array au, 
    and rho=rho0 to infinity analytically using an approximation.
    """

    if (verbose): print 'GetA(',epsilon,')...'

    # Make sure the array starts at 0.0 (direction not important)
    rhomin = min(aRho)
    rho0 = max(aRho)
    assert(rhomin == 0.0)
    assert(rho0 > 0.0)

    # Integrate area under curve numerically (rho=0 to rho0)
    ay = au**2  # vectorized
    Aintegral = IntegrateDiscrete(aRho,ay)
    
    # Calculate area from tail analytically (rho=rho0 to infinity).
    # These come from the approximation for u at large rho,
    #     u=B*exp(-s*rho).
    # so Atail = Integral(u**2,rho,rho0,infinity) = B**2/(2s)*exp(-2s rho0)
    # just need to match at rhomax so that u(rho0) is continuous
    # uleft = uright, with uright = B*exp(-s*rho0), so
    # B = uleft / exp(-s*rho0)
    s = sqrt(-epsilon)
    uleft = au[-1]
    B = uleft / exp(-s*rho0)
    Atail = B**2/(2*s)*exp(-2*s*rho0)
    
    # Add both to get total area
    A = Aintegral+Atail

    if (verbose): 
        print '      rho0',rho0
        print '      B',B
        print '      Aintegral',Aintegral
        print '      Atail',Atail
        print '      A',A

    assert(A>0)  # otherwise could end up with complex numbers down the line...
    
    return A



def GetEigenenergy(epsilonMin,epsilonMax,L,nsteps):
    
    """
    Find an eigenenergy between the specified values by obtaining 
    successively better approximations and then extrapolating to a final value.
    """

    if (verbose): print 'GetEigenenergy(',epsilonMin,epsilonMax,L,nsteps,')...'
        
    # Define 5 rho0 values to use for extrapolating to rho=infinity.
    # put in descending order so rho inverse is in ascending order!
    # Dr. Mayes requirement: pick 3 values between 100 and 10. 
    aRho0 = array([100,50,25,16,10],'d')
    aRho0inv = 1.0/aRho0 # get 1/rho0, vectorized!
    npoints = len(aRho0)
    aNegEpsilon = zeros(npoints,'d') # create empty array
    
    # Find energy value estimate for each rho0
    for i in range(npoints): # do i=0,npoints-1 
        rho0 = aRho0[i]
        aNegEpsilon[i] = -GetEnergy(epsilonMin,epsilonMax,L,rho0,nsteps)

    if (verbose):
        print 'aRho0:',aRho0
        print 'aRho0inv:',aRho0inv
        print 'aNegEpsilon:',aNegEpsilon
    
    # Get interpolation values for later plotting
    ax = arange(0.0,0.1,0.0011,'d') 
    ndata = len(ax)
    ay = zeros(ndata,'d')
    if (verbose): print 'ax',ax
    for i in range(ndata):
        x = ax[i]
        ay[i] = Interpolate(aRho0inv, aNegEpsilon, x, npoints)
        
    # Now extrapolate to rho0=infinity to find the best energy estimate.
    # Note: polynomial extrapolation doesn't work very well here!
    rho0inv = 0.0  # for rho0=infinity
    negepsilon = Interpolate(aRho0inv, aNegEpsilon, rho0inv, npoints)
    
    print 'negepsilon:',negepsilon
    print 'aRho0:',aRho0
    print 'aRho0inv:',aRho0inv
    print 'aNegEpsilon:',aNegEpsilon
    print 'ax',ax
    print 'ay',ay

    # plot results
    gplt.plot([rho0inv],[negepsilon],'title "extrapolated" with points pointtype 4 pointsize 1')
    gplt.hold('on')
    gplt.plot(aRho0inv,aNegEpsilon,'title "calculated" with points pointtype 1 pointsize 1')
    gplt.plot(ax,ay,'title "interpolated" with lines')
    gplt.title('-Epsilon vs 1/rho0')
    gplt.xtitle('1/rho0')
    gplt.ytitle('-Epsilon')
    gplt.grid('on')
    gplt.legend('left bottom')
    gplt.output('epsilon vs rho0inv.gif',gifTerminal)
    gplt.hold('off')
    
    epsilon = -negepsilon
    return epsilon
  




def GetRexact(n,L,aRho):
    """
    Analytical solutions to Schroedinger equation.
    See Griffiths qm p 141, rho = r/a. 
    Note: a0^-3/2 disappears in changing variables from R(r) to R(rho).
    """
    if n==1:
        aR = 2*exp(-aRho) # vectorized
    elif n==2:
        if L==0:
            aR = 1.0/sqrt(2.0) * (1-aRho/2.0) * exp(-aRho/2.0) # vectorized
        elif L==1:
            aR = 1.0/sqrt(24.0) * aRho * exp(-aRho/2.0) # vectorized
    elif n==3:
        if L==0:
            aR = 2.0/sqrt(27) * (1.0-2.0*aRho/3.0+2.0*aRho**2/27) * exp(-aRho/3.0)
        elif L==1:
            aR = 8.0/27/sqrt(6) * (1-aRho/6.0) * aRho * exp(-aRho/3.0)
        elif L==2:
            aR = 4.0/81/sqrt(30) * aRho**2 * exp(-aRho/3.0)
    return aR



    

#----------------------------------------------------------------------------
# main program
#----------------------------------------------------------------------------

print 
print 'program hw5 running...'
print 'date: ',now()
print 

#verbose = True    
verbose = False

folder = 'results\\'

# Loop over n (principle quantum number)
#for n in range(2,2):
for n in range(1,4):

    # Specify the range where we want to find an eigenenergy
    # (we know there's one at -1/4, so look in -1/2 to -1/8)
    epsilonMin = -1.0/2  # -0.5
    epsilonMax = -1.0/8  # -0.125
    
    # Find an eigenenergy between the specified range.
    # This will extrapolate from some approximate values to what should
    # be an exact value (though polynomial extrapolation doesn't work 
    # too well here)
    nsteps = 10**3
    L = 0 # doesn't actually matter - energy should be same for all values of L
    epsilon = GetEigenenergy(epsilonMin,epsilonMax,L,nsteps)

    # Just use the exact energy...
#    epsilon = -1.0/n**2  # in Rydbergs
    
    # Get the radial eigenfunctions R associated with this eigenenergy.
    # Note: Psi(r,angles) = R(r) * Y(angles)
    verbose = True
    rhostart = 0.0
    rhostop = 20.0
    nsteps = 10**4
    ndata = 200
    for L in range(n): # do L=0,n-1
    
        state = "%d%d" % (n,L)  # eg 20
        
        # Get eigenfunction u (normalized)
        aRho, au, A = SolveSchroedinger(epsilon,L,rhostart,rhostop,nsteps,ndata)
    
        # Get radial eigenfunction R from R = u/rho
        aR = au / aRho  # vectorized!
        # this will get divide by zero at zero, so need to replace with 
        # correct limiting value.
        # As rho->0, u = rho**(L+1), so for L=0, u=rho, for L=1, u=rho**2
        # and R=u/rho, so for L=0, u=1, for L=1, u=rho=0
        # but these are non-normalized values, so need to divide by sqrt(A)
        aR[0] = choose(L,(1.0,0.0,0.0)) / sqrt(A)
    
        # Get exact solution for comparison
        aRho2 = arange(rhostart, rhostop, 1.0)
        aRexact = GetRexact(n,L,aRho2)

        # Get probability density
        aP = 4*pi*au**2  # vectorized!

        print 'Eigenstate %s:' % state
        print 'aRho:',aRho
        print 'au:',au
        print 'aR:',aR
        print 'aRexact:',aRexact

        uname = 'u' + state
        Rname = 'R' + state
        Pname = 'P' + state
        ## uname = 'u_{%s}' % state
        ## Rname = 'R_{%s}' % state
        ## Pname = 'P_{%s}' % state

        # Plot u, R, and P vs rho
        statename = 'Hydrogen' + state
        title = 'Hydrogen state nl=' + state
        gplt.plot(aRho,au, 'title "%s" with lines linetype 0' % uname)
        gplt.hold('on')
        gplt.plot(aRho,aR, 'title "%s" with lines linetype 0' % Rname)
        gplt.plot(aRho2,aRexact, 'title "%s exact" with points pointtype 1 pointsize 1' % Rname)
        gplt.plot(aRho,aP, 'title "%s" with lines linetype 2' % Pname)
        gplt.title(title)
        gplt.xtitle('rho')
        gplt.ytitle('Dimensionless u, R, P')
        gplt.grid('on')
        gplt.output(folder + statename +'.gif',gifTerminal)
        ## gplt.output(folder + statename +'.ps',psTerminal)
        gplt.hold('off')



# Plot u(0) vs Epsilon
if 1:
    verbose = False
    ndata = 200
    epsilonMin = -1.0
    epsilonMax = 0.0
    h=(epsilonMax-epsilonMin)/ndata
    aEpsilon = arange(epsilonMin, epsilonMax, h)
    au0 = zeros(ndata,'d')
    au1 = zeros(ndata,'d')
    au2 = zeros(ndata,'d')
    rho0 = 20.0
    rhomin = rho0
    rhomax = 0.0
    nsteps = 10**3
    for i in range(ndata):
        epsilon = aEpsilon[i]
        L = 0
        au0[i] = IntegrateSchroedingerToZero(epsilon,L,rho0,nsteps)
        L = 1
        au1[i] = IntegrateSchroedingerToZero(epsilon,L,rho0,nsteps)
        L = 2
        au2[i] = IntegrateSchroedingerToZero(epsilon,L,rho0,nsteps)
        print i,epsilon,au0[i]
    print aEpsilon
    print au0
    print au1
    gplt.plot(aEpsilon, au0, 'title "u0" with lines')
    gplt.hold('on')
    gplt.plot(aEpsilon, au1, 'title "u1" with lines')
    gplt.plot(aEpsilon, au2, 'title "u2" with lines')
    gplt.title('u(0) for different L values vs Epsilon')
    gplt.xtitle('Epsilon')
    gplt.ytitle('u(0)')
    gplt.grid('on')
    gplt.legend('left top')
    gplt.output(folder + 'u vs epsilon.gif',gifTerminal)
    gplt.hold('off')



