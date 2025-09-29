      subroutine vumat(nblock, ndir, nshr, nstatev, nfieldv, nprops,
     1  lanneal, stepTime, totalTime, dt, cmname, coordMp, charLength,
     2  props, density, strainInc, relSpinInc,
     3  tempOld, stretchOld, defgradOld, fieldOld,
     4  stressOld, stateOld, enerInternOld, enerInelasOld,
     5  tempNew, stretchNew, defgradNew, fieldNew,
     6  stressNew, stateNew, enerInternNew, enerInelasNew)
C***********************************************************************
C     
      include 'vaba_param.inc'
C
      dimension props(nprops), density(nblock), coordMp(nblock,*),
     1  charLength(nblock), strainInc(nblock,ndir+nshr),
     2  relSpinInc(nblock,nshr), tempOld(nblock),
     3  stretchOld(nblock,ndir+nshr),
     4  defgradOld(nblock,ndir+nshr+nshr),
     5  fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatev), enerInternOld(nblock),
     7  enerInelasOld(nblock), tempNew(nblock),
     8  stretchNew(nblock,ndir+nshr),
     8  defgradNew(nblock,ndir+nshr+nshr),
     9  fieldNew(nblock,nfieldv),
     1  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     2  enerInternNew(nblock), enerInelasNew(nblock)
C
      character*80 cmname

      real*8 F(3,3), U(3,3), UInv(3,3), R(3,3), FOld(3,3)

      real*8 F_h0ApOld(3,3), F_h0ApNew(3,3)
      real*8 F_h0ApNewInv(3,3), F_h0Ae(3,3)


      real*8 sigma_h0AOld(3,3), sigma_h0ANew(3,3)
      real*8 sigma_h0BOld(3,3), sigma_h0BNew(3,3)

      real*8 sigma_h1AOld(3,3), sigma_h1ANew(3,3)

      real*8 sigma_sAOld(3,3), sigma_sANew(3,3)

      real*8 sigmaVol_h0A(3,3), sigmaVol_h0B(3,3)
      real*8 sigmaVol_h1A(3,3)
      real*8 sigmaVol_sA(3,3)

      real*8 scorotNew(3,3), sigmaNew(3,3)
      real*8 EYE(3,3)


      real*8 mu_h0A, s0_h0A, G_h0A, gammaDot0_h0A
      real*8 s_h0AOld, s_h0ANew, gammaDot_h0A
      real*8 epspl_h0AOld, epspl_h0ANew

      real*8 mu_h0B, Im_h0B

      real*8 mu_h1A, Im_h1A

      real*8 mu_sA, Im_sA

      real*8 k_sA, k_h0A, k_h0B, k_h1A

      real*8 enerDev_h0A, enerDev_h0B
      real*8 enerVol_h0A, enerVol_h0B

      real*8 enerDev_h1A 
      real*8 enerVol_h1A

      real*8 enerDev_sA
      real*8 enerVol_sA

      real*8 vf_param1, vf_param2, rho_h0, rho_h1

      real*8 vf_hInit, xi_h0Old, xi_h0New, vf_h0, vf_h1, vf_s

      real*8 J, J_dens

      real*8 kb,theta
      real*8 zero

C     variables for reading vf data from text file
      integer nLines
      real*8, dimension (:,:), allocatable :: vfData
      real*8  Xpt, Ypt, Zpt

c     initialize zero parameter
      parameter(zero=0.d0)
C     identity matrix
      call ONEM(EYE)

C     material properties
      mu_sA = props(1)
      Im_sA = props(2)

      mu_h0A = props(3)
      s0_h0A = props(4)
      G_h0A = props(5)
      gammaDot0_h0A = props(6)
      
      mu_h0B = props(7)

      mu_h1A = props(8)
      Im_h1A = props(9)

      k_sA = props(10)
      k_h0A = props(11)
      k_h0B = props(12)

      k_h1A = props(13)

      vf_param1 = props(14)
      vf_param2 = props(15)
      rho_h0 = props(16)
      rho_h1 = props(17)

      kb = props(18)
      theta = props(19)

C     read vf data from text file
      if((totalTime.eq.zero).and.(stepTime.eq.zero)) then
        call numLinesVolfracData('volume fraction data file',nLines)
        allocate (vfData(nLines,5))
        call readVolfracData('volume fraction data file',nLines,vfData)
      endif

      do 100 km = 1,nblock

C       stores new deformation gradient in F
        F(1,1) = defgradNew(km,1)
        F(2,2) = defgradNew(km,2)
        F(3,3) = defgradNew(km,3)        
        F(1,2) = defgradNew(km,4)
        F(2,3) = defgradNew(km,5)
        F(3,1) = defgradNew(km,6)
        F(2,1) = defgradNew(km,7)
        F(3,2) = defgradNew(km,8)
        F(1,3) = defgradNew(km,9)
        
C       stores old deformation graident in FOld
        FOld(1,1) = defgradOld(km,1)
        FOld(2,2) = defgradOld(km,2)
        FOld(3,3) = defgradOld(km,3)
        FOld(1,2) = defgradOld(km,4)
        FOld(2,3) = defgradOld(km,5)
        FOld(3,1) = defgradOld(km,6)
        FOld(2,1) = defgradOld(km,7)
        FOld(3,2) = defgradOld(km,8)
        FOld(1,3) = defgradOld(km,9)

C       stores new stretch tensor in U
        U(1,1) = stretchNew(km,1)
        U(2,2) = stretchNew(km,2)
        U(3,3) = stretchNew(km,3)
        U(1,2) = stretchNew(km,4)
        U(2,1) = stretchNew(km,4)
        U(2,3) = stretchNew(km,5)
        U(3,2) = stretchNew(km,5)
        U(3,1) = stretchNew(km,6)
        U(1,3) = stretchNew(km,6)       

C       R=F*UInv
        call M3INV(U,UInv)
        R = matmul(F,UInv)

        if((totalTime.eq.zero).and.(stepTime.eq.zero)) then

          Xpt = coordMp(km,1)
          Ypt = coordMp(km,2)
          Zpt = coordMp(km,3)

          call initVolfrac(Xpt,Ypt,Zpt,nLines,vfData,vf_hInit)

C         print some of the data during initialization 
          if (km.eq.nblock) then
            write(6,*) km
            write(6,*) coordMp(km,1)
            write(6,*) coordMp(km,2)
            write(6,*) coordMp(km,3)
            write(6,*) vf_hInit
          endif

C         initialize F_h0ApOld
          stateOld(km,1:3) = 1.0D0
          stateOld(km,4:6) = 0.0D0

C         initialize sigma_h0AOld
          stateOld(km,7:12) = 0.0D0
C         initialize sigma_h0BOld
          stateOld(km,13:18) = 0.0D0
C         initialize sigma_h1AOld
          stateOld(km,19:24) = 0.0D0     
C         initialie sigma_sAOld
          stateOld(km,25:30) = 0.0D0

C         initialize sigmaVol_h0A
          stateOld(km,31:36) = 0.0D0
C         initialize sigmaVol_h0B
          stateOld(km,37:42) = 0.0D0
C         initialize sigmaVol_h1A
          stateOld(km,43:48) = 0.0D0
C         initialize sigmaVol_sA
          stateOld(km,49:54) = 0.0D0


C         initialize enerDev_h0A
          stateOld(km,55) = 0.0D0
C         initialize enerDev_h0B
          stateOld(km,56) = 0.0D0
C         initialize enerDev_h1A
          stateOld(km,57) = 0.0D0
C         initialize enerDev_sA
          stateOld(km,58) = 0.0D0


C         initialize enerVol_h0A
          stateOld(km,59) = 0.0D0
C         initialize enerVol_h0B
          stateOld(km,60) = 0.0D0
C         initialize enerVol_h1A
          stateOld(km,61) = 0.0D0
C         initialize enerVol_sA
          stateOld(km,62) = 0.0D0


C         initialize s_h0AOld
          stateOld(km,63) = s0_h0A
C         initialize gammaDot_h0A
          stateOld(km,64) = 0.0D0
C         initialize epspl_h0AOld
          stateOld(km,65) = 0.0D0
C         initialize vf_hInit
          stateOld(km,66) = vf_hInit
C         initialize xi_h0Old
          stateOld(km,67) = 1.0D0
C         intialize xi_h0*vf_hInit
          stateOld(km,68) = vf_hInit
C         initialize J
          stateOld(km,69) = 1.0D0
C         initialize J_dens
          stateOld(km,70) = 1.0D0

C         intial elastic response
          
C         deviatoric stress for sA
          call gentDeviatoric(F,Im_sA,mu_sA,sigma_sANew,enerDev_sA)

C         deviatoric stress for h0A
          call elDeviatoric(F,mu_h0A,sigma_h0ANew,enerDev_h0A)
C         deviatoric stress for h0B
          call nhDeviatoric(F,mu_h0B,sigma_h0BNew,enerDev_h0B)

C         volumetric stress for sA
          call volumetric2(F,k_sA,sigmaVol_sA,enerVol_sA)

C         volumetric stress for h0A
          call volumetric2(F,k_h0A,sigmaVol_h0A,enerVol_h0A)
C         volumetric stress for h0B
          call volumetric2(F,k_h0B,sigmaVol_h0B,enerVol_h0B)

          vf_h0 = vf_hInit
          vf_s = 1.0D0 - vf_hInit

          sigmaNew = 
     1    vf_h0*(sigma_h0ANew+sigma_h0BNew+sigmaVol_h0A+sigmaVol_h0B)+
     2    vf_s*(sigma_sANew+sigmaVol_sA)

        else

          F_h0ApOld(1,1) = stateOld(km,1)
          F_h0ApOld(2,2) = stateOld(km,2)
          F_h0ApOld(3,3) = stateOld(km,3)        
          F_h0ApOld(1,2) = stateOld(km,4)
          F_h0ApOld(2,1) = stateOld(km,4)
          F_h0ApOld(2,3) = stateOld(km,5)
          F_h0ApOld(3,2) = stateOld(km,5)        
          F_h0ApOld(3,1) = stateOld(km,6)
          F_h0ApOld(1,3) = stateOld(km,6)

          sigma_h0AOld(1,1) = stateOld(km,7)
          sigma_h0AOld(2,2) = stateOld(km,8)
          sigma_h0AOld(3,3) = stateOld(km,9)        
          sigma_h0AOld(1,2) = stateOld(km,10)
          sigma_h0AOld(2,1) = stateOld(km,10)
          sigma_h0AOld(2,3) = stateOld(km,11)
          sigma_h0AOld(3,2) = stateOld(km,11)        
          sigma_h0AOld(3,1) = stateOld(km,12)
          sigma_h0AOld(1,3) = stateOld(km,12)

          s_h0AOld = stateOld(km,63)
          epspl_h0AOld = stateOld(km,65)
          vf_hInit = stateOld(km,66)
          xi_h0Old = stateOld(km,67)

          enerDev_h0A = stateOld(km,55)
          enerDev_h0B = stateOld(km,56)
          enerDev_h1A = stateOld(km,57)


          call PLINC4(F,F_h0ApOld,vf_hInit*xi_h0Old,sigma_h0AOld,
     1      dt,s_h0AOld,s_h0ANew,
     2      F_h0ApNew,gammaDot_h0A,G_h0A,kb,theta,gammaDot0_h0A,
     3      epspl_h0AOld,epspl_h0ANew)
          call M3INV(F_h0ApNew,F_h0ApNewInv)
          F_h0Ae = matmul(F,F_h0ApNewInv)


          call gentDeviatoric(F,Im_sA,mu_sA,sigma_sANew,enerDev_sA)
          call elDeviatoric(F_h0Ae,mu_h0A,sigma_h0ANew,enerDev_h0A)
          call nhDeviatoric(F,mu_h0B,sigma_h0BNew,enerDev_h0B)
          call gentDeviatoric(F,Im_h1A,mu_h1A,sigma_h1ANew,enerDev_h1A)

          call volumetric2(F,k_sA,sigmaVol_sA,enerVol_sA)
          call volumetric_hard(F,k_h0A,rho_h0,rho_h1,xi_h0Old,         
     1                         sigmaVol_h0A,enerVol_h0A)
          call volumetric_hard(F,k_h0B,rho_h0,rho_h1,xi_h0Old,         
     1                         sigmaVol_h0B,enerVol_h0B)
          call volumetric_hard(F,k_h1A,rho_h0,rho_h1,xi_h0Old,         
     1                         sigmaVol_h1A,enerVol_h1A)


          call volfrac15(F,gammaDot_h0A,xi_h0Old,xi_h0New,dt,
     1                   rho_h0,rho_h1,vf_param1,vf_param2,J,
     2                   enerDev_h0A+enerDev_h0B,enerDev_h1A,J_dens)

          vf_h0 = vf_hInit*xi_h0New
          vf_h1 = vf_hInit*(1.0D0 - xi_h0New)
          vf_s =  1.0D0 - vf_hInit

          sigmaNew = vf_h0*(sigma_h0ANew+sigma_h0BNew+sigmaVol_h0A+
     1             sigmaVol_h0B)+
     2             vf_h1*(sigma_h1ANew+sigmaVol_h1A)+
     3             vf_s*(sigma_sANew+sigmaVol_sA)

C          enerInternNew(km) = enerDev_sA + enerDev_sB + enerDev_h0A +
C     1                        enerDev_h0B + enerDev_h1B +
C     2                        enerVol_sA + enerVol_sB +
C     3                        enerVol_h0A + enerVol_h0B + 
C     4                        enerVol_h1B          

          call enerIntern(F,FOld,sigmaNew,density(km),dt,
     1       enerInternOld(km),enerInternNew(km))


        endif

C       pass state variables to Abaqus *********************************

        stateNew(km,1) = F_h0ApNew(1,1)
        stateNew(km,2) = F_h0ApNew(2,2)
        stateNew(km,3) = F_h0ApNew(3,3)
        stateNew(km,4) = F_h0ApNew(1,2)
        stateNew(km,5) = F_h0ApNew(2,3)
        stateNew(km,6) = F_h0ApNew(3,1)


        stateNew(km,7) = sigma_h0ANew(1,1)
        stateNew(km,8) = sigma_h0ANew(2,2)
        stateNew(km,9) = sigma_h0ANew(3,3)
        stateNew(km,10) = sigma_h0ANew(1,2)
        stateNew(km,11) = sigma_h0ANew(2,3)
        stateNew(km,12) = sigma_h0ANew(3,1)

        stateNew(km,13) = sigma_h0BNew(1,1)
        stateNew(km,14) = sigma_h0BNew(2,2)
        stateNew(km,15) = sigma_h0BNew(3,3)
        stateNew(km,16) = sigma_h0BNew(1,2)
        stateNew(km,17) = sigma_h0BNew(2,3)
        stateNew(km,18) = sigma_h0BNew(3,1)

        stateNew(km,19) = sigma_h1ANew(1,1)
        stateNew(km,20) = sigma_h1ANew(2,2)
        stateNew(km,21) = sigma_h1ANew(3,3)
        stateNew(km,22) = sigma_h1ANew(1,2)
        stateNew(km,23) = sigma_h1ANew(2,3)
        stateNew(km,24) = sigma_h1ANew(3,1)

        stateNew(km,25) = sigma_sANew(1,1)
        stateNew(km,26) = sigma_sANew(2,2)
        stateNew(km,27) = sigma_sANew(3,3)
        stateNew(km,28) = sigma_sANew(1,2)
        stateNew(km,29) = sigma_sANew(2,3)
        stateNew(km,30) = sigma_sANew(3,1)


        stateNew(km,31) = sigmaVol_h0A(1,1)
        stateNew(km,32) = sigmaVol_h0A(2,2)
        stateNew(km,33) = sigmaVol_h0A(3,3)
        stateNew(km,34) = sigmaVol_h0A(1,2)
        stateNew(km,35) = sigmaVol_h0A(2,3)
        stateNew(km,36) = sigmaVol_h0A(3,1)

        stateNew(km,37) = sigmaVol_h0B(1,1)
        stateNew(km,38) = sigmaVol_h0B(2,2)
        stateNew(km,39) = sigmaVol_h0B(3,3)
        stateNew(km,40) = sigmaVol_h0B(1,2)
        stateNew(km,41) = sigmaVol_h0B(2,3)
        stateNew(km,42) = sigmaVol_h0B(3,1)

        stateNew(km,43) = sigmaVol_h1A(1,1)
        stateNew(km,44) = sigmaVol_h1A(2,2)
        stateNew(km,45) = sigmaVol_h1A(3,3)
        stateNew(km,46) = sigmaVol_h1A(1,2)
        stateNew(km,47) = sigmaVol_h1A(2,3)
        stateNew(km,48) = sigmaVol_h1A(3,1)

        stateNew(km,49) = sigmaVol_sA(1,1)
        stateNew(km,50) = sigmaVol_sA(2,2)
        stateNew(km,51) = sigmaVol_sA(3,3)
        stateNew(km,52) = sigmaVol_sA(1,2)
        stateNew(km,53) = sigmaVol_sA(2,3)
        stateNew(km,54) = sigmaVol_sA(3,1)


        stateNew(km,55) = enerDev_h0A
        stateNew(km,56) = enerDev_h0B
        stateNew(km,57) = enerDev_h1A
        stateNew(km,58) = enerDev_sA

        stateNew(km,59) = enerVol_h0A
        stateNew(km,60) = enerVol_h0B
        stateNew(km,61) = enerVol_h1A
        stateNew(km,62) = enerVol_sA

        stateNew(km,63) = s_h0ANew
        stateNew(km,64) = gammaDot_h0A
        stateNew(km,65) = epspl_h0ANew
        stateNew(km,66) = vf_hInit
        stateNew(km,67) = xi_h0New
        stateNew(km,68) = xi_h0New*vf_hInit

        stateNew(km,69) = J
        stateNew(km,70) = J_dens

C       pass sigmaNew to Abaqus **************************************** 

C       store components of sigmaNew in stressNew
        stressNew(km,1) = sigmaNew(1,1)
        stressNew(km,2) = sigmaNew(2,2)
        stressNew(km,3) = sigmaNew(3,3)

        stressNew(km,4) = sigmaNew(1,2)
        stressNew(km,5) = sigmaNew(2,3)
        stressNew(km,6) = sigmaNew(3,1)

  100 continue

      return
      end

      include '20190910_UR.for'

      include 'numLinesVolfracData.for'
      include 'readVolfracData.for'
      include 'initVolfrac.for'

      include 'elDeviatoric.for'
      include 'gentDeviatoric.for'
      include 'nhDeviatoric.for'
      include 'volumetric2.for'
      include 'volumetric_hard.for'

      include 'PLINC4.for'

      include 'volfrac15.for'

      include 'enerIntern.for'
      include 'enerInelas.for'