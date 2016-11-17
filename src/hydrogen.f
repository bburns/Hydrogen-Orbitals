!----------------------------------------------------------------------------
! Homework 5: Hydrogen Atom Orbitals
! Brian Burns
! PHYS 4350 Computational Physics
! Dr. Mayes
! May 7, 2004
!----------------------------------------------------------------------------
! Find the eigenenergy and eigenfunction for the first excited state of Hydrogen. 
! Some variables: 
!    rho = distance from nucleus in Bohr radii = r/a0
!    epsilon = energy of electron in Rydbergs
!    Psi(r,angles) = R(r) * Y(angles)
!    Use change of variable u(r) = R(r) * r
! We'll solve for u(rho) to get R(rho), then plot R(rho) vs rho.
!----------------------------------------------------------------------------
     
      program hw5
    
      implicit none

      real*8 pi, epsilon0, me, e, hbar, c, eV
      real*8 alpha, a0, E1, Rydberg
    
      real*8 epsilon, rhostart, rhostop
      real*8 epsilonMin, epsilonMax
      integer n, L, nsteps
      character*80 gifTerminal, psTerminal  
      character*80 filename   !,folder
      integer ndata
      parameter(ndata=200) ! number of plot points for numeric solution
      real*8 aRho(0:ndata), au(0:ndata), aR(0:ndata), aP(0:ndata)
      integer ndata2
      parameter(ndata2=20) ! number of plot points for exact solution
      real*8 aRho2(0:ndata2), aRexact(0:ndata2)
      real*8 h
      real*8 A
      integer i
      real*8 GetEigenenergy
    
      ! global variables
      logical verbose
      common /globals/verbose
    
      ! fundamental constants
      pi = 3.1415926535897931d0
      epsilon0 = 8.85419d-12  ! C**2/Jm = (surface charge density / Efield)
      me = 9.10939d-31 ! kg
      e = 1.60218d-19 ! C
      hbar = 1.05457d-34 ! Js
      c = 2.99792d8 ! m/s
      eV = 1.60218d-19 ! J
      
      ! derived constants
      alpha = e**2d0 / (4d0*pi*epsilon0*hbar*c) ! Fine structure const = 1/137 (unitless)
      a0 = hbar / (alpha * me * c)  ! Bohr radius = 0.529e-10m
      E1 = - alpha**2d0 * me * c**2d0 / 2d0 ! Ground state energy of Hydrogen atom = -13.6 eV
      Rydberg = E1  ! -13.6 eV
      ! En = E1/n**2
            
      gifTerminal = 'gif large size 1024,760'
      psTerminal = 'postscript enhanced color'

      ! This will turn on debugging output
      verbose = .false.


!----------------------------------------------------------------------------
!     Main program
!----------------------------------------------------------------------------

      print * 
      print *, 'program hw5 running...'
      print * 
      
      
      nsteps = 10**7
      
      ! Loop over n (principle quantum number)
      do n=1,3
      
          if (n==2) then
              ! Specify the range where we want to find an eigenenergy
              ! (we know there's one at -1/4, so look in -1/2 to -1/8)
              epsilonMin = -0.125d0
              epsilonMax = -0.5d0
              ! Find an eigenenergy between the specified range.
              ! This will extrapolate from some approximate values to what should
              ! be an exact value (though polynomial extrapolation doesn't work 
              ! too well here)
              L = 0 ! doesn't actually matter - energy is same for all values of L
              epsilon = GetEigenenergy(epsilonMin,epsilonMax,L,nsteps)
              print *,'Found energy:'
              print *,n,epsilonMin,epsilonMax,epsilon
          else
              ! Just use the exact energy...
              epsilon = -1.0d0/n**2d0  ! in Rydbergs
          endif

          ! Loop over L (azimuthal quantum number)
          ! Get the radial eigenfunctions R associated with this eigenenergy.
          ! Note: Psi(r,angles) = R(r) * Y(angles)
          rhostart = 0d0
          rhostop = 20d0
          do L=0,n-1
          
              ! Get eigenfunction u (normalized)
              call SolveSchroedinger(epsilon,L,rhostart,rhostop, nsteps, &
                   & ndata, aRho, au, A)
          
              ! Get radial eigenfunction R from R = u/rho
              do i=1,ndata
                  aR(i) = au(i) / aRho(i)
              enddo
              ! this will get divide by zero at zero, so need to replace with 
              ! correct limiting value.
              ! As rho->0, u = rho**(L+1), so for L=0, u=rho, for L=1, u=rho**2
              ! and R=u/rho, so for L=0, u=1, for L=1, u=rho=0
              ! but these are non-normalized values, so need to divide by sqrt(A)
              aR(0) = 0d0
              if (L==0) aR(0) = 1d0/sqrt(A)
          
              ! Get probability density
              do i=0,ndata
                  aP(i) = 4d0*pi*au(i)**2d0
              enddo
      
              ! Get exact solution for comparison
              h = (rhostop - rhostart) / ndata2
              do i=0,ndata2
                  aRho2(i) = rhostart + h * i
                  call GetRexact(ndata2,n,L,aRho2,aRexact)
              enddo
      
              print *, 'Eigenstate ', n,L
              if (verbose) then
                  print *, 'aRho:',aRho
                  print *, 'au:',au
                  print *, 'aR:',aR
                  print *, 'aRexact:',aRexact
              endif 
            
              ! Get filename
              ! g77 doesn't support encode, so have to do this...
100     format(A,I1,I1,a,$)
              write(unit=filename,fmt=100) 'Hydrogen',n,L,'.txt'  ! eg Hydrogen10.txt
              
              ! Save data to file
              open (unit=12,file=filename,status='unknown')
              write (12,*) '# numeric solution (rho, u, R, P)'
              do i=0,ndata
                  write (12,'(4(1PD24.14))') aRho(i),au(i),aR(i),aP(i)
              enddo
              write (12,*) ''
              write (12,*) ''
              write (12,*) '# exact solution (rho, Rexact)'
              do i=0,ndata2
                  write (12,'(2(1PD24.14))') aRho2(i),aRexact(i)
              enddo
              close (12)

          enddo
      enddo
      
      end
      
      

      subroutine GetRexact(ndata,n,L,aRho,aRexact)
          
          ! Get analytical solutions to Schroedinger equation for Hydrogen atom.
          ! See Griffiths qm p141, rho = r/a. 
          ! Note: a0**-3/2 disappears in changing variables from R(r) to R(rho).
          ! ndata = size of aRho and aRexact arrays (0 to ndata)
          ! n = principle quantum number
          ! L = azimuthal quantum number
          ! aRho = array with rho values
          ! aRexact = output array holding exact solution
          
          integer ndata, n, L
          real*8 aRho(0:ndata), aRexact(0:ndata)
          
          ! local variables
          integer i
          
          ! yes this is slow, but lost vectorized arrays that were in python
          do i=0,ndata
              if (n==1) then
                  aRexact(i) = 2d0*exp(-aRho(i))
              elseif (n==2) then
                  if (L==0) then
                      aRexact(i) = 1d0/sqrt(2d0) * (1d0-aRho(i)/2d0) * exp(-aRho(i)/2d0) 
                  else ! L==1
                      aRexact(i) = 1d0/sqrt(24d0) * aRho(i) * exp(-aRho(i)/2d0) 
                  endif
              elseif (n==3) then
                  if (L==0) then
                      aRexact(i) = 2d0/sqrt(27d0) * (1d0-2d0*aRho(i)/3d0+2d0*aRho(i)**2d0/27d0) * exp(-aRho(i)/3d0)
                  elseif (L==1) then
                      aRexact(i) = 8d0/27d0/sqrt(6d0) * (1d0-aRho(i)/6d0) * aRho(i) * exp(-aRho(i)/3d0)
                  elseif (L==2) then
                      aRexact(i) = 4d0/81d0/sqrt(30d0) * aRho(i)**2d0 * exp(-aRho(i)/3d0)
                  endif
              endif
          enddo
          
          end
      
      
      
!--------------------------------------------------------------------------
!     Schroedinger routines
!--------------------------------------------------------------------------
      
      function upp(rho, u, uprime, epsilon, L)   ! ie g(x,y,yprime,args)
      
          ! This is effectively the Shroedinger equation.
          ! upp stands for u''
          ! Calculates u'' at some point rho, u, u'.
          ! From u'' = -(epsilon + 2z/rho -L(L+1)/rho^2)u, with z=1 for Hydrogen.
          ! Note: R(rho) = u(rho)/rho.
          ! args should contain epsilon and L
          
          implicit none
          real*8 upp, rho, u, uprime
          real*8 epsilon
          integer L
          
          ! Handle bad behavior at rho=0, using approximation for small rho,
          !    u = rho^(L+1), so u'' becomes
          !    upp = -epsilon*rho^(L+1) - 2z*rho^L + L(L+1)*rho^(L-1)
          ! which can be simplified for the various cases of L and rho=0:
          if (rho==0d0) then
              if (L==0) then
                  upp = -2d0   ! ie -2z
              else if (L==1) then
                  upp = 2d0
              else
                  upp = 1d0
              endif
          else
              upp = -(epsilon + 2d0/rho - L*(L + 1d0)/rho**2d0)*u
          endif
          
          return
          end
      
      
      function IntegrateSchroedingerToZero(epsilon,L,rho0,nsteps,ndata) 
          
          ! Get value for u at 0 by integrating the Schroedinger equation
          ! from rho0 to 0 (going backwards). 
          ! Note: R(rho) = u(rho)/rho.
          
          implicit none
          real*8 IntegrateSchroedingerToZero
          real*8 epsilon, rho0
          integer L, nsteps, ndata
          logical verbose
          common /globals/verbose
          real*8 aRho(0:ndata),au(0:ndata)
          real*8 A
          
          ! local variables
          real*8 rhostart, rhostop,u
          
          if (verbose) then 
              print *, 'IntegrateSchroedingerToZero(',epsilon,L,rho0,nsteps,ndata,')...'
          endif
              
          ! Must have enough points here for discrete integral to work for normalization,
          ! and nsteps must be evenly divisible by ndata to work properly!
          rhostart = rho0
          rhostop = 0d0
          call SolveSchroedinger(epsilon,L,rhostart,rhostop, &
                       nsteps,ndata,aRho,au,A)
          u = au(ndata)  ! get last value, which is value of integral, ie value of u at rho=0
          if (verbose) then
              print *, '    u(0)=',u
          endif
          IntegrateSchroedingerToZero = u
          return
          end
      
      
      subroutine SolveSchroedinger(epsilon,L,rhostart,rhostop,nsteps,ndata,aRho,au,A) 
          
          ! Get the normalized eigenfunction u as a function of position, over the 
          ! range rho=rhostart to rhostop. returns arrays aRho and au, and 
          ! normalization factor A. 
          ! Note: u=R*r
          ! ndata determines size of arrays produced.
          
          implicit none
          real*8 epsilon, rhostart, rhostop
          integer L, nsteps, ndata
          real*8 aRho(0:ndata),au(0:ndata),A
          logical verbose
          common /globals/verbose
          
          ! local variables
          real*8 u0, u0prime, s
          external upp
          real*8 GetA
          integer i
          
          if (verbose) then 
              print *, 'SolveSchroedinger(',epsilon,L,rhostart,rhostop,nsteps,ndata,')...'
          endif
      
          ! Calculate the start values
          if (rhostart==0d0) then
              ! Start values for rho=0, using approximation for u at small rho,
              !    u = rho^(L+1)
              ! so
              !    u' = (L+1)*rho^L
              u0 = 0d0
              u0prime = 0d0
              if (L==0) u0prime = 1d0
          else
              ! Start values for large rho
              ! These come from the approximation for u at large rho,
              !     u=B*exp(-s*rho)
              ! with B a scaling factor that will be accounted for when u is normalized.
              ! Note: +/- u0's are equivalent, since only probability density is physically 
              ! measurable, but the negative one matches the plot in Griffiths.
              s = sqrt(-epsilon)
              u0 = -exp(-s*rhostart)
              u0prime = -s * u0
          endif 
          
          if (verbose) then 
              print *, '   u0',u0
              print *, '   u0prime',u0prime
          endif 
      
      
          ! Do 2nd order integration from rhostart to rhostop
          !                                 (g, xmin, xmax, y0, yprime0, nsteps, ndata, args)
          call Solve2ndOrder(upp, rhostart,rhostop,u0,u0prime, nsteps, &
                     & ndata, aRho, au, epsilon, L)
      
          ! Now normalize our u array
          A = GetA(epsilon, aRho, au, ndata)
          do i=0,ndata
              au(i) = au(i)/sqrt(A)
          enddo
          if (verbose) print *,'  normalized by A=',A
          end
      
      
      function GetA(epsilon, aRho, au, ndata) 
          
          ! Calculate the integral of u**2 from rho=0 to infinity.
          ! This is used in normalizing u (and hence Psi).
          ! Calculates integral from rho=0 to rho0 numerically using the array au, 
          ! and rho=rho0 to infinity analytically using an approximation.

          implicit none
          real*8 GetA, epsilon
          integer ndata
          real*8 aRho(0:ndata), au(0:ndata)
          logical verbose
          common /globals/verbose
          
          ! local variables
          real*8 rhomin, rhomax, rho0
          real*8 ay(0:ndata)
          real*8 Aintegral, Atail, A, s, uleft, B
          integer i
          real*8 IntegrateDiscrete
          
          if (verbose) then
              print *, 'GetA(',epsilon,')...'
          endif
      
          ! Make sure the array starts at 0.0 (direction not important)
          if (aRho(0)==0d0) then
              rhomin = 0d0
              rhomax = aRho(ndata)
          else
              rhomin = aRho(ndata)
              rhomax = aRho(0)
          endif
      
          ! Integrate area under curve numerically (rho=0 to rho0)
          do i=0,ndata
              ay(i)=au(i)**2d0
          enddo
          Aintegral = IntegrateDiscrete(aRho,ay,ndata)
          
          ! Calculate area from tail analytically (rho=rho0 to infinity).
          ! These come from the approximation for u at large rho,
          !     u=B*exp(-s*rho).
          ! so Atail = Integral(u**2,rho,rho0,infinity) = B**2/(2s)*exp(-2s rho0)
          ! just need to match at rhomax so that u(rho0) is continuous
          ! uleft = uright, with uright = B*exp(-s*rho0), so
          ! B = uleft / exp(-s*rho0)
          s = sqrt(-epsilon)
          uleft = au(ndata)
          B = uleft / exp(-s*rho0)
          Atail = B**2d0/(2d0*s)*exp(-2d0*s*rho0)
          
          ! Add both to get total area
          A = Aintegral+Atail
      
          if (verbose) then 
              print *, '      rho0',rho0
              print *, '      B',B
              print *, '      Aintegral',Aintegral
              print *, '      Atail',Atail
              print *, '      A',A
          endif 
          
          GetA = A
          return
          end


      
      function GetEigenenergy(epsilonMin,epsilonMax,L,nsteps)
          
          ! Find an eigenenergy between the specified values by obtaining 
          ! successively better approximations and then extrapolating to a final value.

          implicit none
          real*8 GetEigenenergy
          real*8 epsilonMin, epsilonMax
          integer L, nsteps
          logical verbose
          common /globals/verbose
          
          ! local variables
          integer i,npoints,np
          parameter (npoints=5,np=103)
          real*8 rho0
          real*8 aRho0(1:npoints), aRho0inv(1:npoints),aNegEpsilon(1:npoints)
          real*8 GetEnergy, Interpolate
          real*8 rho0inv
          real*8 negepsilon
          external Interpolate
          character*80 filename
          integer ndata
          real*8 ax(1:np),ay(1:np)
          real*8 h,x
          
          if (verbose) then 
              print *, 'GetEigenenergy(',epsilonMin,epsilonMax,L,nsteps,')...'
          endif
              
          ! Define 5 rho0 values to use for extrapolating to rho=infinity.
          ! put in descending order so rho inverse is in ascending order.
          ! Dr. Mayes requirement: pick 3 values between 100 and 20 (was 10). 
          data aRho0/100d0,70d0,40d0,30d0,20d0/
          do i=1,npoints
              aRho0inv(i) = 1d0/aRho0(i)
          enddo
          
          ! Find energy value estimate for each rho0
          ! Lots of time gets spent in this loop
          do i=1,npoints
              rho0 = aRho0(i)
              aNegEpsilon(i) = -GetEnergy(epsilonMin,epsilonMax,L,rho0,nsteps)
          enddo
      
          if (verbose) then
              print *, 'aRho0:',aRho0
              print *, 'aRho0inv:',aRho0inv
              print *, 'aNegEpsilon:',aNegEpsilon
          endif
          
          ! Get interpolation values for later plotting (rho0inv = 0 to 1/20)
          h = (1d0/20d0 - 0d0)/np
          ndata = npoints
          do i=1,np
              x = 0.0001d0 + h * i
              ay(i) = Interpolate(ndata,aRho0inv,aNegEpsilon,x,npoints)
              ax(i) = x
          enddo
          if (verbose) print *, 'ax',ax

          ! Now extrapolate to rho0=infinity to find the best energy estimate.
          rho0inv = 0d0  ! for rho0=infinity
          ndata = npoints ! size of arrays aRho0inv and aNegEpsilon
          negepsilon = Interpolate(ndata,aRho0inv,aNegEpsilon,rho0inv,npoints)

          if (verbose) then
              print *, 'negepsilon:',negepsilon
              print *, 'aRho0:',aRho0
              print *, 'aRho0inv:',aRho0inv
              print *, 'aNegEpsilon:',aNegEpsilon
              print *, 'ax',ax
              print *, 'ay',ay
          endif
      
          ! Save data to file
          filename = 'Epsilon.txt'
          open (unit=12,file=filename,status='unknown')
          write (12,*) '# epsilon vs 1/rho0'
          do i=1,npoints
              write (12,'(2(1PD24.14))') aRho0inv(i),aNegEpsilon(i)
          enddo
          write (12,*) ''
          write (12,*) ''
          write (12,*) '# extrapolated to rho0=infinity'
          write (12,'(2(1PD24.14))') rho0inv,negepsilon
          write (12,*) ''
          write (12,*) ''
          write (12,*) '# epsilon vs 1/rho0 (interpolated values)'
          do i=1,np
              write (12,'(2(1PD24.14))') ax(i),ay(i)
          enddo
          close(12)

          GetEigenenergy = -negepsilon
          return
          end


        
      function GetEnergy(epsilonMin,epsilonMax,L,rho0,nsteps) 
          
          ! Find the energy (eigenenergy) for which u(0)=0, by integrating 
          ! our Schroedinger equation from the starting point rho0 to the origin.
          ! If our energy is correct the SE should integrate to zero at the origin,
          ! since that is one of the boundary conditions. The other BC is at the 
          ! starting point rho0, for which we have approximate values for u and u'.
          ! The search will begin between the specified values for the energy. 
          ! nsteps is used in the integrations.
          
          implicit none
          real*8 GetEnergy
          real*8 epsilonMin, epsilonMax, rho0
          integer L, nsteps
          logical verbose
          common /globals/verbose
          
          ! local variables
          integer niterations
          real*8 xchangemax, epsilon
          real*8 FindZero
          external IntegrateSchroedingerToZero
          
          if (verbose) then 
              print *, 'GetEnergy(',epsilonMin,epsilonMax,L,rho0,nsteps,')'
          endif
      
          niterations = 15
          xchangemax = 1d-12
          
          ! Call FindZero to find epsilon
          epsilon = FindZero(IntegrateSchroedingerToZero,epsilonMin,&
                    & epsilonMax,xchangemax,niterations,L,rho0,nsteps)
                
          GetEnergy = epsilon
          return
          end


      
!--------------------------------------------------------------------------
!     Integration routines
!--------------------------------------------------------------------------
      
      subroutine Solve2ndOrder(g,xmin,xmax,y0,yprime0,nsteps,ndata, &
                        & ax,ay,epsilon,L)
                                  
          ! Integrate the function g from xmin to xmax and return arrays for x and y.
          ! Integrate the function g from xmin to xmax and end
          ! Uses nsteps for number of steps, ndata+1 for number of data points in arrays. 
          ! y0 and yprime0 are initial values. 
          ! The function g should return y'' for x,y,y'
          ! The function g should end
          ! args is an arbitrary argument that will be passed to the function g. 
          ! Uses 5th order Runge-Kutta.

          implicit none
          real*8 g, xmin, xmax, y0, yprime0
          real*8 epsilon
          integer L
          integer nsteps, ndata
          real*8 ax(0:ndata),ay(0:ndata)
          logical verbose
          common /globals/verbose
          
          ! local variables
          real*8 h, x, y, yprime
          integer nskip,i,j
          real*8 k0,k1,k2,k3,m0,m1,m2,m3

          if (verbose) then 
              print *, 'Solve2ndOrder(',xmin,xmax,y0,yprime0,nsteps, &
                        & ndata,epsilon,L,')...'
          endif
          
          ! get step size
          h = (xmax - xmin) / nsteps
          
          ! get number of points to skip between writing to arrays
          nskip = nsteps/ndata    
          
          ! initial point
          x = xmin
          y = y0
          yprime = yprime0
          
          if (verbose) then
              print *, 'nsteps=',nsteps
              print *, 'ndata=',ndata
              print *, 'nskip=',nskip
              print *, 'h(step distance)=',h
              print *, '    j,x,y:'
          endif
          
          ! save first point to arrays
          j=0
          ax(j) = x
          ay(j) = y
          if (verbose) print *, '    ',j,x,y
          
          ! loop over number of steps
          do i=1,nsteps
              
              k0 = h * yprime
              m0 = h * g(x, y, yprime, epsilon, L)
              k1 = h * (yprime+m0/2d0)
              m1 = h * g(x+h/2d0, y+h*yprime/2d0, yprime+m0/2d0, epsilon, L)
              k2 = h * (yprime+m1/2d0)
              m2 = h * g(x+h/2d0, y+h*yprime/2d0+h*m0/4d0, yprime+m1/2d0, epsilon, L)
              k3 = h * (yprime+m2)
              m3 = h * g(x+h, y+h*yprime+h*m1/2d0, yprime+m2, epsilon, L)
              
              ! order of y and yprime calculation is important here!!
              x = xmin + h*i
              y = y + h*yprime+h*(m0+m1+m2)/6d0
              yprime = yprime + (m0+2d0*m1+2*m2+m3)/6d0
              
              ! save x,y to arrays every <nskip> data points, including last data point
              if (mod(i,nskip)==0) then
                  j=j+1
                  ax(j)=x
                  ay(j)=y
                  if (verbose) print *, '    ',j,x,y
              endif
          enddo
          
          end
          
      
      
      function IntegrateDiscrete(ax,ay,ndata) 
          
          ! Find the area under the curve specified by the given x and y arrays.
          ! Equal spacing in x array is assumed. 

          implicit none
          real*8 IntegrateDiscrete
          integer ndata
          real*8 ax(0:ndata),ay(0:ndata)
          logical verbose
          common /globals/verbose
          
          ! local variables
          integer i
          real*8 h, area, a
          
          if (verbose) print *,'IntegrateDiscrete(',ndata,')'
          
          !assert(ndata > 1)
          
          h = abs(ax(1)-ax(0)) ! assuming equal spacing between points
          if (verbose) print *,'  h',h

          ! Simpson's 1/3 rule
          area = 0d0
          do i=0, ndata-2, 2  ! Note: must step through every other point
              a = h * (ay(i) + 4*ay(i+1) + ay(i+2)) / 3d0  ! Simpson's 1/3 Rule
              area = area + a
          enddo
          if (verbose) print *,'  area simpson',area

          IntegrateDiscrete = area
          return
          end


      
!--------------------------------------------------------------------------
!     Interpolation routines
!--------------------------------------------------------------------------
      
      function FindZero(f,xmin,xmax,xchangemax,niterations,L,rho0,nsteps) 
          
          ! Find the zero for the given function f that lies between xmin and xmax, 
          ! ie try to find value for x for which f(x)=0. 
          ! xchangemax - max change between iterations
          ! niterations - max number of iterations
          ! args - can be any arguments needed by the function f 

          implicit none
          real*8 FindZero
          real*8 f
          real*8 xmin, xmax, xchangemax
          real*8 rho0
          integer L,nsteps
          integer niterations, ndata
          external f
          logical verbose
          common /globals/verbose
          
          ! local variables
          real*8 x1, x2, x3, f1, f2, f3
          real*8 x, x3old, deltaX, deltaF
          integer i
          real*8 slopeInverse
          
          x1 = xmin
          x2 = xmax
          
          ndata = 200 ! must be large enough for discrete integration of area
          
          ! Get first two values for f
          f1 = f(x1,L,rho0,nsteps,ndata)
          f2 = f(x2,L,rho0,nsteps,ndata)
          if (verbose) then
              print *, 0,x1,f1
              print *, 1,x2,f2
          endif
          
          ! Try to find value for x that will give f=0
          ! ie iterate until f converges to zero
          x3old = 0d0
          do i=0,niterations-1 
              ! find new x and f(x)
              deltaX = x2-x1
              deltaF = f2-f1
              slopeInverse = deltaX / deltaF
              x3 = x1-f1*slopeInverse
              f3 = f(x3,L,rho0,nsteps,ndata)
              if (verbose) print *, i+2,x3,f3
              ! replace one of the endpoints with out new values
              if (f3<0d0 .and. f1<0d0) then
                  x1=x3
                  f1=f3
              else
                  x2=x3
                  f2=f3
              endif
              ! exit when change in x is < xchangemax
              if (abs(x3-x3old) < xchangemax) exit !break
              x3old = x3
          enddo
              
          x = x3
          
          FindZero = x
          return
          end
          


      function Interpolate(ndata,ax,ay,x,npoints)
      
          ! Interpolate or extrapolate from the given data to find the value
          ! for y at the point x.
          ! ndata specifies the size of the given arrays, ax and ay (1 to ndata). 
          ! npoints is the number of points to use in the extrapolation. 
          ! Arbitrary spacing in the ax array is okay, but should be sorted in ascending order!
          ! Uses Lagrangian interpolation with <npoints> data points:
          !   y(x) = y1*l1(x) + y2*l2(x) + ... yn*ln(x)
          !   li(x) = (x-x1)/(xi-x1) * (x-x2)/(xi-x2) * ... (x-xn)/(xi-xn) [excluding ith term]
          ! y(x) will be a polynomial of order <npoints-1>
    
          implicit none
          real*8 interpolate
          integer ndata
          real*8 ax(1:ndata),ay(1:ndata)
          real*8 x
          integer npoints
          logical verbose
          common /globals/verbose
    
          ! local variables
          integer ibest(41)/41*0/ ! 41=2*npointsmax+1
          integer i,j
          integer ni,nj
          real*8 y
          real*8 li
    
          if (verbose) print *,'interpolate/extrapolate to x=',x, &
                    & ', npoints=',npoints
    
          ! Find best data points to use, based on which are closest to 
          ! requested point x. Will find <npoints> best data points and 
          ! store index of points in array ibest(1:npoints)
          call FindBest(ndata,ax,x,npoints,ibest)

          ! Build the polynomial y(x), evaluated at the point x.
          y=0d0
          do i=1,npoints
                li=1d0
                ni=ibest(i) ! index to ith best point
                ! build up li(x) term, evaluated at the point x
                do j=1,npoints
                      if (i .ne. j) then ! exclude j=i term
                            nj=ibest(j) ! index to jth best point
                            li=li*(x-ax(nj))/(ax(ni)-ax(nj))
                      endif
                enddo
                y=y+ay(ni)*li
          enddo
    
          Interpolate=y
    
          return
          end
      

      subroutine FindBest(ndata,ax,x,npoints,ibest)

          ! Find <npoints> data points closest to requested point x 
          ! and store index to points in array ibest(1:npoints).
          ! ndata specifies the size of the given array ax (1 to ndata). 
          
          implicit none
          integer ndata,npoints
          real*8 ax(1:ndata)
          real*8 x
          integer ibest(*)
          logical verbose
          common /globals/verbose
          
          ! local variables
          integer BSearch ! binary search function
          integer i,j,iclosest,imin,imax
          integer ncandidates,npoints2
          real*8 xdelta(41)/41*0d0/
          
          ! find index of point closest to x
          iclosest = BSearch(ndata,ax,x)
          if (verbose) print *,'closest point at x=',ax(iclosest)
    
          ! Get npoints points in each direction, and find distance 
          ! from x for each. 
          ! This will handle cases where point is at start or end of 
          ! data set, or where all closest points lie in one direction.
          imin = iclosest-npoints
          imax = iclosest+npoints
          ! make sure imin and imax are in array range
          if (imin < 1) imin=1
          if (imax > ndata) imax=ndata
          ncandidates = imax-imin+1
          if (verbose) then
                print *,'candidate points:', ncandidates
                print *,'    j,i,xdata(i),xdelta(j):'
          endif
          do i = imin,imax
                j = i-imin+1
                xdelta(j) = abs(ax(i)-x) ! distance from x
                if (verbose) print *,'    ',j,i,ax(i),xdelta(j)
          enddo
    
          ! Sort points by xdelta, in ascending order
          call IndexSort(ncandidates,xdelta,ibest)
          
          ! Exclude closest point if it's actually the point we're searching for
          ! (dr mayes requirement)
          npoints2 = ncandidates
          if (xdelta(ibest(1)) .eq. 0d0) then
                if (verbose) print *,'excluding point with xdelta=0'
                ! reduce number of available candidates by one
                ncandidates = ncandidates - 1
                ! make sure we don't have more points than candidates
                npoints2 = ncandidates 
                ! shift points down by one
                do i=1,ncandidates
                      ibest(i)=ibest(i+1)
                enddo
          endif

          ! make sure we don't have more points than candidates
          if (npoints2 > ncandidates) npoints2 = ncandidates

          ! Adjust ibest array to correct range
          ! note: at this point the first <npoints> is all we need
!          do i=1,npoints
          do i=1,npoints2
                ibest(i)=ibest(i)+imin-1
          end do
    
          if (verbose) then
                print *,'best points (sorted by xdelta):'
                print *,'    i,ibest,xdata,xdelta'
                do i=1,npoints2
                      print *,'    ',i,ibest(i),ax(ibest(i)), &
                                & abs(x-ax(ibest(i)))
                enddo
          endif
    
          return
          end
      

      subroutine IndexSort(ndata,ax,iorder)
      
          ! Sorts array ax(1:ndata) into ascending numerical order,
          ! storing results in index array iorder(1:ndata).
          ! Uses bubble sort, so keep number of points low.
          
          implicit none
          integer ndata
          real*8 ax(*)
          integer iorder(*)
          logical verbose
          common /globals/verbose
          
          ! local variables
          integer i,j
          integer ni,nj
          real*8 xi,xj
          
          ! initialize iorder array
          do i=1,ndata
              iorder(i)=i
          enddo
          
          ! bubble sort
          do i = 1,ndata
                ni = iorder(i)
                xi = ax(ni)
                do j = i+1,ndata
                      nj = iorder(j)
                      if (ni .ne. nj) then
                            xj = ax(nj)
                            if (xi > xj) then
                                  iorder(i) = nj
                                  iorder(j) = ni
                                  ni = nj
                                  xi = xj
                                  ! xj will get reassigned next round so don't worry about it here
                            endif
                      endif
                enddo
          enddo
    
          end
    
    
      function BSearch(ndata,ax,x)
      
          ! Do binary search to find the index of the point closest to 
          ! the given point x.
          ! Array ax(1:ndata) should be in ascending order.
    
          implicit none
          integer BSearch
          integer ndata
          real*8 ax(1:ndata)
          real*8 x
          logical verbose
          common /globals/verbose
          
          ! local variables
          integer i,imin,imax
          
          imin = 1
          imax = ndata
          do
                i = (imin+imax)/2
                if (x < ax(i)) then
                      imax = i
                else
                      imin = i
                endif
                if (imin+1 >= imax) then ! done
                      ! use imin or imax, depending on which is closer to x
                      if (abs(x-ax(imin)) < abs(x-ax(imax))) then
                            i = imin
                      else
                            i = imax
                      endif
                      exit
                endif
          enddo
          
          BSearch=i
          
          return
          end
    
