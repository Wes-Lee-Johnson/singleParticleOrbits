c Wes Johnson
c 11/30/2020
c this fortran module contains functions to calculate the 
c linear eigenmodes in a circular tokamak


cccc5ccc10cccc5ccc20cccc5ccc30cccc5ccc40cccc5ccc50cccc5ccc60cccc5ccc70cc
      MODULE lem_funcs
      USE my_subs 
      IMPLICIT NONE
      REAL*8, parameter :: omega = 1D5/omegaSI,r0Cent=a/2.0,
     &                     amplitudeE = 0.01*1D4 *qSI/kt*rhoSI
 
      CONTAINS 

      FUNCTION fj_inside(r,j,m0)
         REAL*8 :: qPrime,n0,delR
         REAL*8 :: fj_inside 
         INTEGER*8,INTENT(IN) :: j,m0
         REAL*8, INTENT(IN) :: r 
         qPrime    = qrDiv(r0Cent)
         n0        = DBLE(NINT(m0/qr(r0Cent))) 
         delR      = 1.0/(n0*qPrime) 
         fj_inside = (r/delR-r0Cent/delR-DBLE(j))
      END FUNCTION fj_inside

      FUNCTION fj(r,j,m0)
         REAL*8 :: fj
         REAL*8, INTENT(IN) :: r
         INTEGER*8,INTENT(IN) :: j,m0 
         fj = exp((-1.0)/(2.0*9.0)*(fj_inside(r,j,m0))**2)
      END FUNCTION fj

      FUNCTION r(rho,z)
         REAL*8, INTENT(IN) :: rho,z
         REAL*8 :: r
         r = sqrt((rho-R0)**2+(z)**2)
      END FUNCTION r

      FUNCTION theta(rho,z)
         REAL*8, INTENT(IN) :: rho,z
         REAL*8 :: theta
         theta = atan2(z,rho-R0)
      END FUNCTION theta

      FUNCTION bigR(x,y)
         REAL*8, INTENT(IN) :: x,y
         REAL*8 :: bigR
         bigR = sqrt(x*x + y*y)
      END FUNCTION bigR

      FUNCTION polodial(x,y)
         REAL*8, INTENT(IN) :: x,y
         REAL*8 :: polodial
         polodial = atan2(y,x)
      END FUNCTION polodial 

      FUNCTION pol2Rect(Ezeta,ER,x,y)
         REAL*8, INTENT(IN) :: Ezeta,ER,x,y
         REAL*8, DIMENSION(2) :: pol2Rect
         pol2Rect(1) = (-y*Ezeta+x*ER)/bigR(x,y)
         pol2Rect(2) = (x*Ezeta+y*ER)/bigR(x,y)
      END FUNCTION pol2Rect

      FUNCTION phi_LEM(r,zeta,theta,time,N_theta,m0)
         REAL*8, INTENT(IN) :: r,zeta,theta,time
         INTEGER*8, INTENT(IN) :: N_theta,m0
         INTEGER*8 :: j
         REAL*8 :: n0,phi_LEM_sum,phi_LEM
         n0 =  DBLE(NINT(DBLE(m0)/qr(r0Cent)))
         do j = -N_theta,N_theta,1
           phi_LEM_sum = phi_LEM_sum+fj(r,j,m0)*
     &              cos(DBLE(m0+j)*theta-omega*time-n0*zeta)
         enddo 
         phi_LEM = An0(r)*phi_LEM_sum 
      END FUNCTION phi_LEM 

      FUNCTION An0(r)
         REAL*8, INTENT(IN) :: r
         REAL*8 :: An0 
         An0 = exp((-1.0)/(2.0*4.0)*(r/a-1.0)**2)
      END FUNCTION An0

      FUNCTION qr(r)
         REAL*8, INTENT(IN) :: r 
         REAL*8 :: qr 
         qr = q0 + delq * (r/a)**2 
      END FUNCTION qr

      FUNCTION qrDiv(r)
         REAL*8, INTENT(IN) :: r
         REAL*8 :: qrDiv
         qrDiv = 2*delq*r/a**2 
      END FUNCTION qrDiv

      FUNCTION calcBmag(rho,z)
         REAL*8, INTENT(IN) :: rho,z 
         REAL*8 :: rs,qrs,grs,calcBmag
         rs = r(rho,z)
         qrs = qr(rs)
         grs = rs/(R0*qrs) 
         calcBmag=(B0*R0/rho)*sqrt(1*grs**2)
      END FUNCTION calcBmag

      FUNCTION bFeildRZ(rho,z)
         REAL*8, INTENT(IN) :: rho,z
         REAL*8 :: alpha,rs,qrs,grs,bRho,bZ
         REAL*8, DIMENSION(2) :: bFeildRZ
         alpha =  R0*B0/rho
         rs    = r(rho,z)
         qrs   = qr(rs)
         grs   = rs/(R0*qrs)
         bRho  = alpha*grs*(-z/rs)
         bZ    = alpha*grs*((rho-R0)/rs)
         bFeildRZ(1) = bRho ; bFeildRZ(2) = bZ
      END FUNCTION bFeildRZ

c####5###10####5###20####5###30####5###40####5###50####5###60####5###70
c#                       Electric Feild 
c####5###10####5###20####5###30####5###40####5###50####5###60####5###70
      FUNCTION eFeildRect(pos,time)
         REAL*8, INTENT(IN) :: time
         REAL*8, DIMENSION(3), INTENT(IN) :: pos 
         INTEGER*8 :: m0,N_theta 
         REAL*8, DIMENSION(3) :: eFeildRect 
         REAL*8, DIMENSION(2) :: temp
         REAL*8 :: rho,zeta 
         REAL*8 :: x,y,z 
         x = pos(1);y = pos(2);z = pos(3)
         rho     = bigR(x,y)
         zeta    = polodial(x,y)
         m0      = 15
         N_theta = 2 
         eFeildRect    = eFeild(rho,zeta,z,time,m0,N_theta) 
         !coordinate transformation:
         temp          = pol2Rect(eFeildRect(2),eFeildRect(1),x,y)
         eFeildRect(1) = temp(1) !x comp
         eFeildRect(2) = temp(2) !y comp 
         eFeildRect    = amplitudeE*eFeildRect
      END FUNCTION eFeildRect

      FUNCTION eFeild(rho,zeta,z,time,m0,N_theta)
         REAL*8, INTENT(IN) :: rho,zeta,z,time
         INTEGER*8, INTENT(IN) :: m0,N_theta
         REAL*8 :: rs,thetas,phiLem,E_rho,E_z,E_zeta
         REAL*8, DIMENSION(3) :: eFeild
         rs    = r(rho,z)
         thetas= theta(rho,z)
         phiLem= phi_LEM(rs,zeta,thetas,time,N_theta,m0)
         E_rho = An0DivRho(rs)*rDivRho(rho,rs)*phiLem+
     &      An0(rs)*Erho(rho,zeta,z,time,m0,N_theta)
         E_z   = An0DivZ(rs)*rDivZ(rho,rs)*phiLem+
     &      An0(rs)*Ez(rho,zeta,z,time,m0,N_theta)
         E_zeta= An0(rs)*Ezeta(rho,zeta,z,time,m0,N_theta)
         eFeild(1) = (-1D0)*E_rho !has wrong direction
         eFeild(2) = (-1D0)*E_zeta!needs to be rotated 
         eFeild(3) = (-1D0)*E_z
      END FUNCTION eFeild

      FUNCTION Erho(rho,zeta,z,time,m0,N_theta)
         REAL*8, INTENT(IN) :: rho,zeta,z,time
         INTEGER*8, INTENT(IN) :: m0,N_theta
         INTEGER*8 :: j 
         REAL*8 :: rs,thetas,n0,Erho
         n0 =  DBLE(NINT(DBLE(m0)/qr(r0Cent)))
         rs     = r(rho,z)
         thetas = theta(rho,z)
         Erho = 0.0
         do j = -N_theta,N_theta,1
           Erho = Erho + fjDivRho(rs,j,m0)*rDivRho(rho,rs)*
     &            cos(DBLE(m0+j)*thetas-n0*zeta-omega*time)+
     &            fj(rs,j,m0)*
     &            sin(DBLE(m0+j)*thetas-n0*zeta-omega*time)*
     &            DBLE(m0+j)*thetaDivRho(rho,z)
         enddo 
      END FUNCTION Erho
    
      FUNCTION rDivRho(rho,rs)
         REAL*8, INTENT(IN) :: rho,rs
         REAL*8 :: rDivRho
         rDivRho = (rho-R0)/rs
      END FUNCTION rDivRho

      FUNCTION thetaDivRho(rho,z)
         REAL*8, INTENT(IN) :: rho,z
         REAL*8 :: rhox,x,numer,denom
         REAL*8 :: thetaDivRho 
         rhox  = rho-R0
         x     = z/rhox
         numer = x/rhox
         denom = 1+x**2
         thetaDivRho = numer/denom
      END FUNCTION thetaDivRho

      FUNCTION fjDivRho(r,j,m0)
         REAL*8, INTENT(IN) :: r
         INTEGER*8, INTENT(IN) :: j,m0
         REAL*8 :: qPrime,n0,delR,gaussian
         REAL*8 :: fjDivRho
         qPrime = qrDiv(r0Cent)
         n0 =  DBLE(NINT(DBLE(m0)/qr(r0Cent)))
         qPrime = qrDiv(r0Cent)
         delR   = 1/(n0*qPrime)
         gaussian = 
     &        exp((-1.0)*(fj_inside(r,j,m0))**2)
         fjDivRho = gaussian*
     &              ((-2.0)*(fj_inside(r,j,m0))/delR)
      END FUNCTION fjDivRho

      FUNCTION Ez(rho,zeta,z,time,m0,N_theta)
         REAL*8, INTENT(IN) :: rho,zeta,z,time
         INTEGER*8, INTENT(IN) :: m0,N_theta
         INTEGER*8 :: j 
         REAL*8 :: rs,thetas,n0,Ez
         rs     = r(rho,z)
         thetas = theta(rho,z)
         n0 =  DBLE(NINT(DBLE(m0)/qr(r0Cent)))
         Ez = 0.0
         do j = -N_theta,N_theta,1
           Ez = Ez + fjDivZ(rs,j,m0)*rDivZ(z,rs)*
     &         cos((DBLE(m0+j))*thetas-n0*zeta-omega*time)+
     &         (-1D0)*fj(rs,j,m0)*
     &         sin(DBLE(m0+j)*thetas-n0*zeta-omega*time)*
     &         DBLE(m0+j)*thetaDivZ(rho,z)
         enddo 
      END FUNCTION Ez

      FUNCTION rDivZ(z,r)
         REAL*8, INTENT(IN) :: z,r
         REAL*8 :: rDivZ 
         rDivZ = z/r
      END FUNCTION rDivZ

      FUNCTION thetaDivZ(rho,z)
         REAL*8, INTENT(IN) :: rho,z
         REAL*8 :: rhox,x,numer,denom
         REAL*8 :: thetaDivZ
         rhox  = rho-R0
         x     = z/rhox
         numer = 1D0/rhox
         denom = 1D0+x**2
         thetaDivZ = numer/denom
      END FUNCTION thetaDivZ
    
      FUNCTION fjDivZ(rs,j,m0) 
         REAL*8, INTENT(IN) :: rs
         INTEGER*8, INTENT(IN) :: j,m0
         REAL*8 :: fjDivZ 
         fjDivZ = fjDivRho(rs,j,m0)
      END FUNCTION fjDivZ

      FUNCTION An0DivZ(r)
         REAL*8, INTENT(IN) :: r
         REAL*8 :: gauss,An0DivZ
         gauss  = exp((-1.0)/(2.0*4.0)*(r/a-1.0)**2) 
         An0DivZ= gauss*(-1.0)/(2.0*4.0)*(r/a-1.0)*2.0
      END FUNCTION An0DivZ

      FUNCTION An0DivRho(r)
         REAL*8, INTENT(IN) :: r
         REAL*8 :: An0DivRho
         An0DivRho = An0DivZ(r)
      END FUNCTION An0DivRho

      FUNCTION Ezeta(rho,zeta,z,time,m0,N_theta)
         REAL*8, INTENT(IN) :: rho,zeta,z,time
         INTEGER*8, INTENT(IN) :: m0,N_theta
         INTEGER*8 :: j
         REAL*8 :: rs,thetas,n0,Ezeta
         Ezeta = 0.0
         n0 =  DBLE(NINT(DBLE(m0)/qr(r0Cent)))
         rs     = r(rho,z)
         thetas = theta(rho,z)
         do j = -N_theta,N_theta,1
           Ezeta = Ezeta + 
     &            fj(rs,j,m0)*
     &            sin(DBLE(m0+j)*thetas-n0*zeta-omega*time)*n0
         enddo 
         Ezeta = Ezeta/rho
      END FUNCTION Ezeta

c####5###10####5###20####5###30####5###40####5###50####5###60####5###70
c#                 Derivative of Electric Feild 
c####5###10####5###20####5###30####5###40####5###50####5###60####5###70


      FUNCTION divTeFeildRect(pos,time)
         REAL*8, INTENT(IN) :: time
         REAL*8, DIMENSION(3), INTENT(IN) :: pos 
         INTEGER*8 :: m0,N_theta 
         REAL*8, DIMENSION(3) :: divTeFeildRect 
         REAL*8, DIMENSION(2) :: temp
         REAL*8 :: rho,zeta 
         REAL*8 :: x,y,z 
         m0      = 15 
         N_theta = 2 
         x = pos(1);y = pos(2);z = pos(3)
         rho     = bigR(x,y)
         zeta    = polodial(x,y)
         divTeFeildRect    = divTeFeild(rho,zeta,z,time,m0,N_theta) 
         !coordinate transformation:
         temp    = pol2Rect(divTeFeildRect(2),divTeFeildRect(1),x,y)
         divTeFeildRect(1) = temp(1) !x comp
         divTeFeildRect(2) = temp(2) !y comp 
         divTeFeildRect    = amplitudeE*divTeFeildRect
      END FUNCTION divTeFeildRect
      
      FUNCTION divTeFeild(rho,zeta,z,time,m0,N_theta) 
         REAL*8, INTENT(IN) :: rho,zeta,z,time
         INTEGER*8, INTENT(IN) :: m0,N_theta
         REAL*8 :: rs,thetas,divT_E_rho,divT_E_z,divT_E_zeta
         REAL*8, DIMENSION(3) :: divTeFeild
         rs    = r(rho,z)
         thetas= theta(rho,z)
         divT_E_rho = divTErho(rho,zeta,z,time,m0,N_theta)
         divT_E_z   = divTEz(rho,zeta,z,time,m0,N_theta)
         divT_E_zeta= divTEzeta(rho,zeta,z,time,m0,N_theta)
         divTeFeild(1) = (-1D0)*divT_E_rho !has wrong direction
         divTeFeild(2) = (-1D0)*divT_E_zeta!needs to be rotated 
         divTeFeild(3) = (-1D0)*divT_E_z
      END FUNCTION divTeFeild

      FUNCTION divTErho(rho,zeta,z,time,m0,N_theta)
        REAL*8, INTENT(IN) :: rho,zeta,z,time
        INTEGER*8, INTENT(IN) :: m0,N_theta
        INTEGER*8 :: j,n0 
        REAL*8 :: divTErho,rs,thetas,thetaDivRho_,rDivRho_
        REAL*8 :: sin_,cos_,An0_,fj_,fjDivRho_,An0DivRho_,inside
        !need to get rDivRho multiplied 
        rs          = r(rho,z)
        thetas      = theta(rho,z)
        n0 =  DBLE(NINT(DBLE(m0)/qr(r0Cent)))
        An0_        = An0(rs)
        An0DivRho_  = An0DivRho(rs)
        thetaDivRho_= thetaDivRho(rho,z)
        rDivRho_    = rDivRho(rho,rs)
        divTErho = 0d0 
        do j=-N_theta,N_theta
          inside    = DBLE(m0+j)*thetas-n0*zeta-omega*time
          sin_      = dsin(inside)
          cos_      = dcos(inside)
          fj_       = fj(rs,j,m0)
          fjDivRho_ = fjDivRho(rs,j,m0)
          divTErho = divTErho
     &    +          An0DivRho_*rDivRho_*fj_*omega*sin_
     &    +          An0_*fjDivRho_*rDivRho_*omega*sin_
     &    +          An0_*fj_*omega*cos_*thetaDivRho_
        enddo  
c       if(divTErho.gt.1D-3)then 
c               write(*,*)divTErho,'divTErho'
c       endif
      END FUNCTION divTErho

      FUNCTION divTEz(rho,zeta,z,time,m0,N_theta)
        REAL*8, INTENT(IN) :: rho,zeta,z,time
        INTEGER*8, INTENT(IN) :: m0,N_theta
        INTEGER*8 :: j,n0 
        REAL*8 :: divTEz,rs,thetas,thetaDivZ_,rDivZ_
        REAL*8 :: sin_,cos_,An0_,fj_,fjDivZ_,An0DivZ_,inside
        rs          = r(rho,z)
        thetas      = theta(rho,z)
        n0 =  DBLE(NINT(DBLE(m0)/qr(r0Cent)))
        An0_        = An0(rs)
        An0DivZ_    = An0DivZ(rs)
        thetaDivZ_  = thetaDivZ(rho,z)
        rDivZ_      = rDivZ(rho,rs)
        divTEz = 0d0 
        do j=-N_theta,N_theta
          inside    = DBLE(m0+j)*thetas-n0*zeta-omega*time
          sin_      = dsin(inside)
          cos_      = dcos(inside)
          fj_       = fj(rs,j,m0)
          fjDivZ_   = fjDivZ(rs,j,m0)
          divTEz    = divTEz
     &    +          An0DivZ_*rDivZ_*fj_*omega*sin_
     &    +          An0_*fjDivZ_*rDivZ_*omega*sin_
     &    +          An0_*fj_*omega*cos_*thetaDivZ_
        enddo 
c       if(divTEz.gt.1D-3)then 
c               write(*,*)divTEz,'divTEz'
c       endif
      END FUNCTION divTEz

      FUNCTION divTEzeta(rho,zeta,z,time,m0,N_theta)
        REAL*8, INTENT(IN) :: rho,zeta,z,time
        INTEGER*8, INTENT(IN) :: m0,N_theta
        INTEGER*8 :: n0 
        REAL*8 :: divTEzeta,rs,thetas
        REAL*8 :: An0_,phiLem
        rs          = r(rho,z)
        thetas      = theta(rho,z)
        n0 =  DBLE(NINT(DBLE(m0)/qr(r0Cent)))
        An0_        = An0(rs)
        phiLem = phi_LEM(rs,zeta,thetas,time,N_theta,m0)
        divTEzeta = An0_*phiLem*n0*omega/rho 
c       if (divTEzeta.gt.1D-3) then 
c               write(*,*)divTEzeta,'divTEzeta'
c       endif
      END FUNCTION divTEzeta

      END MODULE lem_funcs 
