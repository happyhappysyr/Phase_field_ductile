c     Author: Mojtaba Abdolkhani 
c     Email: mojtababdolkhani@gmail.com
      SUBROUTINE UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL,
     1 DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP,
     2 PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATV, PROPS, NPROPS,
     3 COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER,
     4 KSPT, KSTEP, KINC)
c
      INCLUDE 'ABA_PARAM.INC'
c
      CHARACTER*20 CMNAME
c
      DIMENSION STRESS(NTENS), STATEV(NSTATV), DDSDDE(NTENS, NTENS),
     1 STRAN(NTENS), DSTRAN(NTENS),PROPS(NPROPS), DROT(3, 3),TIME(2),
     2 EPLAS(NTENS), FLOW(NTENS)
c
      PARAMETER(zero=0.D0, one=1.D0, two=2.D0, three=3.D0)
c      
      parameter (six = 6.d0, newton = 20, dsqrt23 = dsqrt(two/three),
     1 toler = 1.0d-8, zeta = 1.d0)
c
      DIMENSION one_hydro(ntens,ntens), dstran_trace(ntens),
     1 one_dev(ntens,ntens), STRESS_MINUS(NTENS), 
     2 DDSDDE_MINUS(NTENS,NTENS)
c      
      REAL*8 EMOD,ENU,EBULK3,EBULK,EG,EG2,EG3,ELAM  
C ======================================================================
c
c     The energy decomposition is triggered via material PROPS(9) in input.
c     Two options are available:   ELASTOPLASTIC_ISO(41) or ELASTOPLASTIC_SD(40)  
c     
c             ELASTOPLASTIC_ISO - no deformation energy decomposition
c             ELASTOPLASTIC_SD  - spherical-deviatoric energy decomposition
C ======================================================================     
      phi=temp+dtemp 
c      
c     
C ======================================================================       
      !Material properties from input file
      EMOD=PROPS(1)
      ENU=PROPS(2)
	  xl=props(3) ! Phase field length scale
      Gc=props(4) ! Toughness
      xK = 0.0001d0
c      
      !Degradation function
      degfnc = zero
      degfnc = (one - xK)*(one - phi)**2 + xK
	  dg=-2.d0*(1.d0-phi)
      ddg=2.d0
c      
      !Computing lame constants
      EBULK3=EMOD/(one-two*ENU)
      EBULK = EBULK3/three
      EG2=EMOD/(one+ENU)
      EG=EG2/two
      EG3=three*EG
      ELAM=(EBULK3-EG2)/three
c      
C ======================================================================  
      !Building material matrix
      DO 20 K1=1,NTENS
        DO 11 K2=1,NTENS
           DDSDDE(K2,K1)=0.0
 11     CONTINUE
 20   CONTINUE
C
      DO 40 K1=1,NDI
        DO 30 K2=1,NDI
           DDSDDE(K2,K1)=ELAM
 30     CONTINUE
        DDSDDE(K1,K1)=EG2+ELAM
 40   CONTINUE
      DO 50 K1=NDI+1,NTENS
        DDSDDE(K1,K1)=EG
 50   CONTINUE
C ======================================================================          
      !RECOVER PLASTIC STRAINS AND ROTATE FORWARD
      !RECOVER EQUIVALENT PLASTIC STRAIN
c            
      eqplas = statev(1)
      call rotsig(statev(1+1:ntens), drot, eplas, 2, ndi, nshr)
C ======================================================================
c      
      STRAN = STRAN + DSTRAN
      STRESS = zero
c
      !Stress predictor
      DO K1=1, NTENS
         DO K2=1, NTENS
             STRESS(K2)=STRESS(K2)+DDSDDE(K2, K1)*(STRAN(K1)-EPLAS(K1))
         END DO
      END DO
c      
      if (type2.eq.40) then
          ! DECOMPOSITION - VOLUMETRIC DEVIATORIC
          ONE_hydro = zero
          A_trace = zero
          dstran_trace = zero
c          
          do i = 1, ndi
              do j = 1, ndi
                  ONE_hydro(i,j) = one/three
              end do
          end do
c          
          ONE_Dev = zero
          do i = 1, ntens
              do j = 1, ntens
                  if (i .eq. j) ONE_dev (i,j) = ONE
              end do
          end do
c      
          ONE_dev = ONE_dev - ONE_hydro      
          dstran_trace = matmul(ONE_hydro, STRAN)
c         
          if (dstran_trace(1) .le. zero) then
              A_trace = one 
c      
              DDSDDE_MINUS = A_trace*matmul(ONE_hydro, DDSDDE)
              DDSDDE = DDSDDE - DDSDDE_MINUS
c      
              STRESS_MINUS = A_trace*matmul(ONE_hydro, STRESS)
              STRESS = STRESS - STRESS_MINUS
          else
              STRESS_MINUS = zero
              DDSDDE_MINUS = zero
          end if
c
      else if (type2.eq.41) then     
          STRESS_MINUS = zero
          DDSDDE_MINUS = zero
          A_trace = zero
      else
          write(6,*)""
          write(6,*)"WARNING: CHECK MATERIAL NAME!"
      end if
C ======================================================================
c      
      !CALCULATE EQUIVALENT VON MISES STRESS
c      
      smises = (stress(1) - stress(2))**2 + (stress(2) - stress(3))**2
     1       + (stress(3) - stress(1))**2
      
      do k1 = ndi + 1, ntens
          smises = smises + six*stress(k1)**2
      enddo
c      
      smises = dsqrt(smises/two)*dsqrt23
C ======================================================================
      !GET YIELD STRESS FROM THE SPECIFIED HARDENING CURVE
c      
      SYIEL0 = PROPS(5)
      Qinf = PROPS(6)
      b = PROPS(7)
	  jeltype1 = PROPS(8) !AT-2 (30) OR TH (31)
	  type2 = PROPS(9) ! SD(40) OR ISO(41)
c      
      SYIELD = SYIEL0 + Qinf*(ONE - EXP(-b*eqplas))
C ======================================================================
c
      !DETERMINE ACTIVE YIELDING
c      
      if (smises .gt. (one + toler)*syield*dsqrt23) then
c          
          !ACTIVE YIELDING
          !SEPARATE HYDROSTATIC STRESS FROM DEVIATORIC
          !CALCULATE FLOW DIRECTION
c          
          shydro = (stress(1) + stress(2) + stress(3))/three
          onesy = one/smises
          do k1 = 1, ndi
              flow(k1) = onesy*(stress(k1) - shydro)
          enddo
c      
          do k1 = 1 + ndi, ntens
              flow(k1) = stress(k1)*onesy
          enddo      
C ======================================================================
c      
          !SOLVE FOR EQUIVALENT VON MISES STRESS AND EQUIVALENT PLASTIC STRAIN INCREMENT USING NEWTON ITERATION
c      
          syiel_0 = syield
          dgama = zero
          hard = zero
c      
          do 130 knewton = 1, newton
c          
              rhs = smises - two*EG*dgama - dsqrt23*syield   
              dgama = dgama + rhs/(two*EG + (two/three)*hard)
c          
              SYIELD = SYIEL0 +Qinf*(ONE-EXP(-b*(eqplas+dsqrt23*dgama)))
              HARD = Qinf*b*EXP(-b*(eqplas+dsqrt23*dgama))
c          
              if (abs(rhs) .lt. toler*syiel0*dsqrt23) go to 10
c      
 130      CONTINUE
          WRITE(6,2) NEWTON
 2        FORMAT(//,30X,'***WARNING - PLASTICITY ALGORITHM DID NOT ',
     1        'CONVERGE AFTER ',I3,' ITERATIONS')
c
          PNEWDT = 0.65d0
c      
10        continue    
C ======================================================================
c
          !UPDATE STRESS, ELASTIC AND PLASTIC STRAINS AND EQUIVALENT PLASTIC STRAIN
c
          do k1 = 1, ndi
              stress(K1) = stress(k1) - 2.d0*EG*dgama*FLOW(K1)
              EPLAS(K1)=EPLAS(K1)+FLOW(K1)*dgama
          enddo
c      
          do k1 = ndi + 1, ntens  
              stress(K1) = stress(k1) - 2.d0*EG*dgama*FLOW(K1)
              EPLAS(K1)=EPLAS(K1)+FLOW(K1)*dgama*TWO
          enddo
c
          eqplas = eqplas + dsqrt23*dgama
c      
          !CALCULATE PLASTIC DISSIPATION
          spd = spd + dsqrt23*dgama*(SYIEL0 + Qinf) + 
     1    (Qinf/b)*(exp(-b*(eqplas + dsqrt23*dgama)) - exp(-b*eqplas))  
C ======================================================================
c      
          !FORMULATE THE JACOBIAN (MATERIAL TANGENT)
c      
          DDSDDE = zero
c      
          theta = one - 2.d0*EG*dgama/(smises)
          theta_bar = (one/(one + hard/EG3)) - (one-theta)
c      
          do k1 = 1, ndi
              do k2 = 1, ndi
                  DDSDDE(k1,k2) = EBULK - 2.d0*EG*theta/three - 
     1        2.d0*EG*theta_bar*flow(k1)*flow(k2)
              end do
              DDSDDE(k1,k1) = DDSDDE(k1,k1) + 2.d0*EG*theta
          end do
c            
          do k1 = ndi + 1, ntens
                  DDSDDE(k1,k1) = EG*theta
     1        - EG*theta_bar*flow(k1)*flow(k1)
          end do
      end if     
C ======================================================================         
      DDSDDE = degfnc*DDSDDE
      STRESS = degfnc*STRESS
c         
      DDSDDE = DDSDDE + DDSDDE_MINUS
      STRESS = STRESS + STRESS_MINUS
C ======================================================================        
c          
      !STORE STRAINS IN STATE VARIABLE ARRAY
c          
      do k1 = 1, ntens
          statev(1 + k1) = eplas(k1)    
      end do
      statev(1) = eqplas  
c      
      SSE = zero      !elastic deformation energy array
c      
      SSE = 0.5d0*ELAM*(one-A_trace)*((stran(1)-eplas(1)) + 
     1 (stran(2)-eplas(2)) + (stran(3)-eplas(3)))**2 + 
     2 EG*((stran(1)-eplas(1))**2 + (stran(2)-eplas(2))**2 + 
     3 (stran(3)-eplas(3))**2) +
     4 0.5d0*EG*((stran(4)-eplas(4))**2)
c	
c   ************************************************************************************************************   
c      
	  !SSE - Elastic energy
	  psiE_old = statev(2+ntens)                   
	  psiE_new = SSE 
c      
	  !SPD - Plastic energy
	  psiE_SPD = SPD
c                  
	  psiE_SSE = zero
	  psiE = zero
c                 
c	  ! AT-2 or TH model
c	  if (jeltype1.eq.30)then      !AT-2 model
c		  psiE_crit = zero
c	  elseif (jeltype1.eq.31)then  !TH model
		  psiE_crit = 0.2652d0*(Gc/xL)
c	  else
c		  write(6,*)"***ERROR - check element types"
c	  end if
c                 
	  ! Elastic energy history field
	  if (psiE_new.gt.psiE_old) then
		  psiE_SSE = psiE_new
	  else
		  psiE_SSE = psiE_old
	  end if
c                  
	  if ((psiE_SSE + psiE_SPD).gt.psiE_crit) then
		  psiE = psiE_SSE + psiE_SPD - psiE_crit
		  psiE = zeta*psiE
	  else
		  psiE = zero
	  end if
c                  
	  statev(2+ntens) = psiE_SSE
c ***************************************************************************************				  
      w=phi**2
      dw=2.d0*phi
      ddw=2.d0
      cw=0.5d0
c     based on 30 and 31 paper:
c     rpl=-(dg*psiE*2.d0*cw/(xl*Gc)+dw/(2.d0*xl**2))
c     drpldt=-(ddg*psiE*2.d0*cw/(xl*Gc)+ddw/(2.d0*xl**2))
      rpl=((psiE/psiE_crit)*((1.d0-phi)/(xl**2)))-(phi/(xl**2))
      drpldt=-((1.d0+(psiE/psiE_crit))/(xl**2))
c 
      RETURN
      END