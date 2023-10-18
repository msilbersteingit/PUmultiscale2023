C     Requires UR.for utility subroutine
      
      subroutine PLINC4(F,FpOld,vf_h,sigmaOld,dt,sOld,sNew,FpNew,
     1       gammaDot,G,kb,theta,gammaDot0,epsplOld,epsplNew)

      implicit none

      real*8 G, kb, theta, gammaDot0 
      real*8 vf_h, tau, dt, gammaDot, sigmaOldDevMag, zero
      real*8 sOld, sNew, sdot, sScaled
      real*8 F(3,3)
      real*8 Fe(3,3), FeInv(3,3), FpOldInv(3,3)
      real*8 FpOld(3,3), sigmaOld(3,3), FpNew(3,3), FpDot(3,3), Dp(3,3)
      real*8 sigmaOldDev(3,3)
      real*8 epspl_dot, Dp_dd_Dp
      real*8 epsplOld, epsplNew

      parameter(zero=0.d0)

C       G = 2.0D-19
C       kb = 1.38064852D-23
C       theta = 293
C       gammaDot0 = 1.0D8

 

C     obtain deviatoric old Cauchy stress (sigmaOldDev)
      call DEVM(sigmaOld,sigmaOldDev)

C     obtain magnitude of sigmaOldDev
      call MAGM(sigmaOldDev,sigmaOldDevMag)      
 
C     calculate equivalent shear stress from old cauchy stress
C     tau = sqrt(0.5 * sigmaOldDev DD sigmaOldDev) 
      call EQUIVS(sigmaOldDev,tau)
       
C      gammaDot = SINH(tau/tau0) / A

      
C      sScaled = 0.1*sOld + 3.0*vf_h*(sOld - 0.1*sOld)

      gammaDot = gammaDot0*EXP( -(G/(kb*theta)) * (1 - (tau/sOld)))


C     increase with scalar plastic strain 
C      sdot = gammaDot*3.0D6
      sdot = 0.0D0

      sNew = sOld + sdot*dt

      if(sigmaOldDevMag.eq.zero) then
        call ZEROM(Dp)
      else
         Dp = gammaDot * sigmaOldDev / sigmaOldDevMag
C        call ZEROM(Dp)
      endif


      call M3INV(FpOld,FpOldInv)
      Fe = matmul(F,FpOldInv)
      call M3INV(Fe,FeInv)

      FpDot = matmul(FeInv,matmul(Dp,F))

      FpNew = FpOld + FpDot*dt


      call DOTPM(Dp,Dp,Dp_dd_Dp)
      epspl_dot = ((2.0D0/3.0D0) * Dp_dd_Dp)**(0.5D0)

      epsplNew = epsplOld + epspl_dot*dt


      end subroutine PLINC4







