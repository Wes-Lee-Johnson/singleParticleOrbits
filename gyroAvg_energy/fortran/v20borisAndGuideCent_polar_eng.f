c Wes Johnson 01/12/2020 
c this program uses the boris integration method to find the motion of a
c charge in a simple torodial feild, zero order ect. 
c this program also integrates the guiding center equations of motion 
c the integration no longer saves every step by swapping old and new  
c calculates the energy of the guiding center particles
c the program uses units q=m=b=1 vThermal=1
c the prgram uses external module to calculate the linear 
c eigenmodes of a cicular tokamak. The partciles are pushed
c in resulting electric feild. 
c this program contains functions to gyro average for the guiding
c center. 01/15/21 

      MODULE internal_funcs
      USE my_subs
      USE lem_funcs
      IMPLICIT NONE
      
      REAL*8, PARAMETER :: PI=4.D0*DATAN(1.D0)
      INTEGER, PARAMETER :: num = 4
       
      CONTAINS 

      FUNCTION vParaDot(mu,q,m,pos)
        REAL*8, DIMENSION(3), INTENT(IN) :: pos
        REAL*8, INTENT(IN) :: mu,m,q
        REAL*8, DIMENSION(3) :: bNorm,gradientB,B
        REAL*8 :: bMag,vParaDot
        INTEGER :: i
        B        = bfeild(pos)
        bMag     = NORM2(B)
        bNorm    = B / bMag
        gradientB= gradB(pos)
        vParaDot = ((-1d0)*mu/m)*dot_product(gradientB,bNorm)
      END FUNCTION vParaDot

      FUNCTION xDot(vParaMag,mu,q,m,pos,time)
        REAL*8, DIMENSION(3), INTENT(IN) :: pos
        REAL*8, INTENT(IN) :: vParaMag,mu,q,m,time
        REAL*8, DIMENSION(3) :: bNorm,gradientB,BcrossGradB,B,E
        REAL*8, DIMENSION(3) :: EcrossB_drift,onRing,polarDrift
        REAL*8, DIMENSION(3) :: xDot,ddtE
        REAL*8 :: bMag,omega,vPerpMag,coeff,gyroRad
        REAL*8, DIMENSION(3,num) :: ring 
        INTEGER :: i 
        B        = bfeild(pos)
        E        = eFeildRect(pos,time)
        bMag     = NORM2(B)
        bNorm    = B / bMag
        omega    = gyroFreq(q,m,bMag)
        vPerpMag = vPerpMagnitude(mu,q,m,bMag)
        gyroRad  = vPerpMag / omega
        ring     = getRing(gyroRad,bNorm,pos) 
        gradientB = gradB(pos)
c       !gyro avg 
        E         = [0d0,0d0,0d0]
        do i = 1,num 
          E        = E        + eFeildRect(ring(:,i),time)
        enddo 
        E        = E        / dble(num) 
        ddtE       = divTeFeildRect(pos,time)
        polarDrift = (m/q/bMag**2)*(ddtE - dot_product(ddtE,bNorm))!perp comp
        BcrossGradB = cross(B,gradientB)/(bMag**2.0)
        coeff = ((vParaMag)**2.0+0.5*(vPerpMag)**2.0)/omega
        EcrossB_drift = cross(E,B)/(bMag**2.0)
        if(norm2(polarDrift) > 1d0)then 
          write(*,*)norm2(polarDrift)/norm2(EcrossB_drift)
        endif 
        xDot = coeff*BcrossGradB+vParaMag*bNorm+EcrossB_drift+polarDrift
      END FUNCTION xDot

      FUNCTION vPerpMagnitude(mu,q,m,bMag)
        REAL*8, INTENT(IN) :: mu,q,m,bMag
        REAL*8 :: vPerpMagnitude
        vPerpMagnitude = sqrt(2.0*mu*bMag/m)
      END FUNCTION vPerpMagnitude

      FUNCTION gyroFreq(q,m,bMag)
        REAL*8, INTENT(IN) :: q,m,bMag
        REAL*8 :: gyroFreq 
        gyroFreq = q / m * bMag
      END FUNCTION gyroFreq
      
cccc5ccc10cccc5ccc20cccc5ccc30cccc5ccc40cccc5ccc50cccc5ccc60cccc5ccc70cc
      FUNCTION getRing(gyroRad,bNorm,pos)
        REAL*8, DIMENSION(3), INTENT(IN) :: pos,bNorm
        REAL*8, INTENT(IN) :: gyroRad
        INTEGER*8 :: i
        REAL*8, DIMENSION(3,num) :: getRing 
        REAL*8, DIMENSION(3) :: vPrime1,vPrime2,v1,v2,rand1,rand2
        REAL*8 :: arg,dot 
        call random_number(rand1)
        call random_number(rand2) 
        !vectors perp to bNorm
        vPrime1 = cross(bNorm,rand1)
        vPrime2 = cross(bNorm,rand2)
        !normalize 
        v1      = vPrime1 / norm2( vPrime1 ) 
        !orthogonalize 
        vPrime2 = vPrime2 - dot_product( vPrime2,v1 ) * v1
        !normalize 
        v2      = vPrime2 / norm2( vPrime2 )
        if (dot_product(v1,v2)>1D-5.or.(dot_product(v1,v1)-1d0)>1D-5
     &      .or.(dot_product(v2,v2)-1d0)>1D-5)then
          write(*,*)dot_product(v1,v2),dot_product(v1,v1),
     &            dot_product(v2,v2)
        endif 
        !calculate ring positions
        do i = 1,num 
          arg = 2d0*PI/dble(num)*dble(i) 
          getRing(:,i) = gyroRad*(dcos(arg)*v1+dsin(arg)*v2)+pos
        enddo 
      END FUNCTION getRing

      END MODULE internal_funcs


c MAIN PROGRAM
      PROGRAM cleanTorodialBfeild
      USE internal_funcs
      USE my_subs
      USE lem_funcs
      IMPLICIT NONE
c VARIABLE DECLARATIONS:
      REAL*8 :: q = 1.0D0,m=1.0D0,qm
      REAL*8 :: bMag,mu,vPerpMag,eng0,engBor0
      REAL*8 :: start,finish
      REAL*8 :: dt,bor_dt,time
      REAL*8 :: vParaDot1st,vParaDot2nd,vParaMagPred
      REAL*8, dimension(3) :: xDot1st,xDot2nd,posPred
      INTEGER :: numSteps = 10**4,numStepsBor = 10**6
      INTEGER :: saveStepBor,saveStep,step,skip=10**0,skipBor=10**2
      REAL*8, dimension(3,1000000)  :: pos,vel
      REAL*8, dimension(3,1000000)  :: posBor,velBor
      REAL*8, dimension(1000000) :: times,vParaMag
      REAL*8, dimension(1000000) :: timesBor,eng,engBor
      REAL*8, dimension(3) :: E,B,vPerp,vPara,bNorm,vParaNorm,vPerpNorm
      REAL*8, dimension(3) :: posBorOld,posBorNew,velBorOld,velBorNew
      REAL*8, dimension(3) :: posOld,posNew
      REAL*8 :: vParaMagOld,vParaMagNew

      call cpu_time(start)
c INTIALIZE:
      write(*,*)"Initialize"
      !set up time of the simulation
      bor_dt = 10D-10*omegaSI
      dt     = 10D-8*omegaSI

      qm = q/m
      posBor(:,1)=[9.5,0.0,0.0]/rhoSI
      velBor(:,1)=[1.0,1.0,1.0]
c     call SRAND(289)
c     velBor(:,1)=[rand(),rand(),rand()]
c     velBor(:,1) = velBor(:,1)/norm2(velBor(:,1))*5.0
      B=bfeild(posBor(:,1))
      bMag = norm2(B)
      bNorm = B / bMag

      vParaMag(1) = 2.0D0
      vPerpMag    = 8.0D0

      vPara = ABS(DOT_PRODUCT(velBor(:,1),bNorm))*bNorm
      vPerp = velBor(:,1)-vPara

      vParaNorm = vPara/norm2(vPara)
      vPerpNorm = vPerp/norm2(vPerp)

      velBor(:,1) = vParaNorm*vParaMag(1) + vPerpNorm*vPerpMag 
      pos(:,1) = posBor(:,1)+cross(vPerpNorm*vPerpMag,bNorm)*m/(bMag*q)
      mu = 0.5*m*(vPerpMag**2.0)/bMag

      !set the particle back 0.5 * dt 
      velBor(:,1)=updateVelocityBoris(velBor(:,1),E,B,qm,(-0.5)*bor_dt)

      !set up temp variables for stepping
      posBorNew    = posBor(:,1)
      velBorNew    = velBor(:,1)
      saveStepBor  = 0 

      posNew       = pos(:,1)
      vParaMagNew  = vParaMag(1)
      saveStep     = 0

      !set the initial electric feild: 
      E = eFeildRect(posBorNew,0D0)

      !calculate the intital energy of the particles: 
      eng0     = 0.5*m*vParaMag(1)**2 + mu*bMag 
      engBor0  = 0.5*m*norm2(velBor(:,1))**2
      eng(1)   = eng0
      engBor(1)= engBor0

c INTEGRATE THE TRAJECTORY 
      write(*,*)"Calculating Guiding Center Trajectory"
      !Loop over the timesteps 
      do step = 1,numSteps -1
        posOld = PosNew 
        vParaMagOld = vParaMagNew 
        
        time = float(step)*dt 
        !update the guiding center using predictor corrector method
        !prediction phase 
        vParaDot1st = vParaDot(mu,q,m,posOld)
        vParaMagPred = vParaMagOld+vParaDot1st*dt
        xDot1st = xDot(vParaMagOld,mu,q,m,posOld,time)
        posPred=posOld+xDot1st*dt 
        !correction phase 
        vParaDot2nd = vParaDot(mu,q,m,posPred)
        vParaMagNew=vParaMagOld+
     &                        0.5*dt*(vParaDot2nd+vParaDot1st)
        xDot2nd=xDot(vParaMagPred,mu,q,m,posPred,time+dt)
        posNew=posOld+0.5*dt*(xDot1st+xDot2nd)
        if (mod(step,skip).eq.0) then
          saveStep = saveStep +1
          !make time for plotting
          times(saveStep)=float(step)*dt
          !save the updated pos and vel 
          pos(:,saveStep)=posNew
          vParaMag(saveStep)=vParaMagNew
          bMag = norm2(bfeild(posNew))
          eng(saveStep) = 0.5*m*vParaMagNew**2 + mu*bMag
        end if 
      end do 

      write(*,*)"Calculating Lorenz Ion Trajectory" 
      do step = 1,numStepsBor -1 
        posBorold = posBorNew
        velBorOld = velBorNew
        !update the velocity boris 
        B = bfeild(posBorOld)
        time = float(step)*bor_dt
        E = eFeildRect(posBorold,time)
        velBorNew = 
     &            updateVelocityBoris(velBorOld,E,B,qm,bor_dt)
        !update position
        posBorNew = posBorOld + velBorNew*bor_dt


        if (mod(step,skipBor).eq.0) then
          saveStepBor = saveStepBor +1
          !make time for plotting
          timesBor(saveStepBor)=float(step)*bor_dt
          !save the updated pos and vel 
          !subtract the gyro radius vector to get to center of orbit
          bMag = norm2(B)
          bNorm= B / bMag
          vPerp= velBorNew - dot_product(bNorm,velBorNew)*bNorm
          posBor(:,saveStepBor)=
     &           posBorNew + cross(vPerp,bNorm)*m/(bMag*q)
          velBor(:,saveStepBor)=velBorNew
          engBor(saveStepBor) = 0.5*m*(norm2(velBorNew))**2
        end if 
      end do 
        
      write(*,*)"lorenz ion:",timesBor(saveStepBor),size(timesBor),
     &          "guiding center:",times(saveStep),size(times)

c SAVE DATA TO FILE:
      write(*,*)"Writting Guiding Center to File"
      OPEN(unit=1,file="../data/guidingCenter.txt",
     &status="UNKNOWN")
      do step = 1,saveStep
        if (mod(step,1).eq.0) then
          WRITE(1,*)times(step)/omegaSI,vel(:,step)*vThermalSI,
     &    pos(:,step)*rhoSI,eng(step)*(vThermalSI**2)*mSI
        end if 
      end do 
      CLOSE(1)
c SAVE DATA TO FILE:
      write(*,*)"Writting Lorenz Ion to File"
      OPEN(unit=1,file="../data/lorenzIon_subRad.txt",
     &status="UNKNOWN")
      do step = 1,saveStepBor
        if (mod(step,1).eq.0) then
          WRITE(1,*)timesBor(step)/omegaSI,velBor(:,step)*vThermalSI,
     &    posBor(:,step)*rhoSI,engBor(step)*(vThermalSI**2)*mSI
        end if 
      end do 
      CLOSE(1)
c TERMINATION:
      call cpu_time(finish)
      write(*,*)"RUNTIME:",finish-start,"[S]"
      STOP
      END PROGRAM cleanTorodialBfeild

