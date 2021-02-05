c Wes Johnson 10/16/2020
c this file contains functions that are used in external files 
c First we define any functions that will be needed
c 09/30/2020 added gradB function
c 10/16/2020 changed units to b=q=m=vThermal 
      MODULE my_subs
      IMPLICIT NONE
      !change the units for numerical stability: 
      REAL*8, parameter :: B0= 1.0,q0= 1.1,delq=3.0
      REAL*8, parameter :: vThermalSI =  1.0D5             ![m/s]
      REAL*8, parameter :: qSI        =  1.6021766D-19     ![C]
      REAL*8, parameter :: mSI        =  1.6726219D-27     ![kg] 
      REAL*8, parameter :: omegaSI    =  qSI/mSI*B0        ![hz]
      REAL*8, parameter :: rhoSI      =  vThermalSI/omegaSI![m]
      REAL*8, parameter :: k          =  1.380649D-23      ![J*K]
      REAL*8, parameter :: R0         =  8.0/rhoSI         
      REAL*8, parameter :: a          =  2.0/rhoSI
      REAL*8, parameter :: kt         =  0.5*mSI*vThermalSI**2     
      ![s] -> [s]*omegaSI
      ![m] -> [m]/rhoSI

      CONTAINS
      
c     function that finds the cross product
      FUNCTION cross(a, b)
c       INTEGER, DIMENSION(3), INTENT(OUT) :: cross
        REAL*8, DIMENSION(3) :: cross
        REAL*8, DIMENSION(3), INTENT(IN) :: a, b
        cross(1) = a(2) * b(3) - a(3) * b(2)
        cross(2) = a(3) * b(1) - a(1) * b(3)
        cross(3) = a(1) * b(2) - a(2) * b(1)
      END FUNCTION cross

c     function that finds "the zero-beta circular, unshifted mag. f"
c     from Scott Parker's SingleParticle.pdf, the coordiantes are 
c     torodial expressed in cartesian 
      FUNCTION bfeild(r_vec)
        REAL*8, DIMENSION(3) :: bfeild
        REAL*8, DIMENSION(3), INTENT(IN) :: r_vec
        REAL*8 :: r,x,y,z,qr,alpha,beta,rho,rhox
        REAL*8, DIMENSION(3) :: theta_hat,phi_hat
        x = r_vec(1)
        y = r_vec(2)
        z = r_vec(3)
        rho= sqrt(x*x+y*y)
        rhox= rho-R0
        r = sqrt((rhox)*(rhox)+z*z)
        qr = q0 + delq*(r/a)*(r/a)
        alpha= R0*B0*(1./rho)
        beta = r/qr/R0
        phi_hat = [(-1.d0)*y,x,0.0d0]*(1.d0/rho)
        theta_hat = [(1.)*z*x/r/rho,(1.)*z*y/r/rho,(-1.)*rhox/r]
        bfeild = alpha*(phi_hat+beta*theta_hat)
      END FUNCTION bfeild


      !this function returns the magnetic dipole feild given a position
      FUNCTION bMagDipole(r_vec)
        REAL*8, DIMENSION(3) :: bMagDipole
        REAL*8, DIMENSION(3), INTENT(IN) :: r_vec
        REAL*8, DIMENSION(3) :: xHat,yHat,zHat
        REAL*8 :: r,x,y,z,B0
        B0 = 300.0! = 3 mu
        x = r_vec(1)
        y = r_vec(2)
        z = r_vec(3)
        xHat = [1.0,0.0,0.0]
        yHat = [0.0,1.0,0.0]
        zHat = [0.0,0.0,1.0]
        r = NORM2(r_vec)
        bMagDipole = B0*(3*x*z/r**5*xHat + 3*y*z/r**5*yHat + 
     &2*(z*z-x*x-y*y)/r**5*zHat)
      END FUNCTION bMagDipole

c     Function that uses the boris method to update velosity given 
c     the electromagnetic feild, particl velocity, particle charge
c     to mass ratio and the time step.
      FUNCTION updateVelocityBoris(vel,E,B,qm,dt)
        !update velosity using the boris method: birdsall 
        !and langdon pg 62 "physics via computer simulation" 
        REAL*8, DIMENSION(3) :: updateVelocityBoris
        REAL*8, DIMENSION(3), INTENT(IN) :: vel,E,B
        REAL*8, INTENT(IN) :: qm,dt
        REAL*8 :: t_mag2
        REAL*8, DIMENSION(3) :: t,s,velMinus,velMinusCrossT,velPrime
        REAL*8, DIMENSION(3) :: velPrimeCrossS,velPlus
        t = qm * B * 0.5 * dt 
        t_mag2 = t(1)*t(1)+t(2)*t(2)+t(3)*t(3)!t squared 
        !calculate s vector 
        s = 2.0*t/(1.0+t_mag2)
        !calculate v minus 
        velMinus=vel+qm*E*0.5*dt 
        !calculate v prime 
        velMinusCrossT=cross(velMinus,t)
        velPrime=velMinus+velMinusCrossT
        !calculate v plus 
        velPrimeCrossS=cross(velPrime,s)
        velPlus=velMinus+velPrimeCrossS
        !calculate v n+1/2a
        updateVelocityBoris=velPlus+qm*E*0.5*dt
      END FUNCTION updateVelocityBoris
      
      !The following function integrates the guiding center equation of
      !motion using forward euler (leap frog if set back 1/2 dt)
      FUNCTION updateGuidingCenter(vPara,bNorm,dt)
        REAL*8, DIMENSION(3) :: updateGuidingCenter
        REAL*8, DIMENSION(3), INTENT(IN) :: vPara,bNorm
        REAL*8, INTENT(IN) :: dt
        updateGuidingCenter = vPara + (NORM2(vPara))*bNorm*dt
      END FUNCTION updateGuidingCenter
      
      !The following function finds the parallel and perpendicular
      !components of the velosity given the vel and B. the output
      !is an array that contains the perp comp first and para
      !comp second and bNrom as thrid 
      FUNCTION getPerpAndParallel(vel,B)
        REAL*8, DIMENSION(3,3) :: getPerpAndParallel
        REAL*8, DIMENSION(3), INTENT(IN) :: vel,B
        REAL*8, DIMENSION(3) :: bNorm,vPerp,vPara
        bNorm = B/NORM2(B)
        vPara = DOT_PRODUCT(vel,bNorm)*bNorm
        vPerp = vel - vPara
        getPerpAndParallel(1,:) = vPerp
        getPerpAndParallel(2,:) = vPara
        getPerpAndParallel(3,:) = bNorm
      END FUNCTION getPerpAndParallel


      !this function calculates the gradient of the torodial B feild 
      !given the position of the particle
      FUNCTION gradB(r_vec)
        REAL*8, DIMENSION(3) :: gradB
        REAL*8, DIMENSION(3), INTENT(IN) :: r_vec
        REAL*8 :: r,x,y,z,qr,rho,rhox
        REAL*8 :: gBx,gBy,gBz
        REAL*8 :: t1,t2,t3,t4,t5,t6,t7 !define repeated terms
        x = r_vec(1)
        y = r_vec(2)
        z = r_vec(3)
        rho= sqrt(x*x+y*y)
        rhox= rho-R0
        r = sqrt((rhox)*(rhox)+z*z)
        qr = q0 + delq*(r/a)*(r/a)

        t1 = (B0*R0)/(2*rho)
        t2 = 2/((R0**2)*(qr**2))
        t3 = (4*delq*qr**2)/((a**2)*(R0**2)*(qr**3))
        t4 = sqrt((r**2)/((R0**2)*(qr**2))+1)
        t5 = (B0*R0)/(rho**3)
        t6 = t1*(t2-t3)/t4
        t7 = t6*(1-R0/rho)-t5*t4

        gBx = t7*x
        gBy = t7*y
        gBz = t6*z
        
        gradB = [gBx,gBy,gBz]
      END FUNCTION gradB


      END MODULE my_subs
