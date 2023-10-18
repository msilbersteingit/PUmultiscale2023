C     requires UR.for utility subroutine

      subroutine enerInelas(Fe,FpOld,FpNew,sigma,dens,dt,oldEn,newEn)

      implicit none

      real*8 dens, dt, oldEn, newEn
      real*8 Fe(3,3), FpOld(3,3), FpNew(3,3), sigma(3,3)
      real*8 FpDot(3,3), FpInv(3,3), Lp(3,3), Dp(3,3), FeInv(3,3)
      real*8 sigmaDoubleDotDp 


      FpDot = (FpNew - FpOld)/dt

      call M3INV(FpNew,FpInv)
      call M3INV(Fe,FeInv)

      Lp = matmul(Fe,matmul(FpDot,matmul(FpInv,FeInv)))
      Dp = 0.5*(Lp + transpose(Lp))

      
      call DOTPM(sigma,Dp,sigmaDoubleDotDp)

      newEn = oldEn + (sigmaDoubleDotDp/dens)*dt

      end