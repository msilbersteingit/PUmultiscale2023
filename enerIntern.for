C     requires UR.for utility subroutine

      subroutine enerIntern(F,FOld,sigma,dens,dt,oldEnInt,newEnInt)

      implicit none

      real*8 dens, dt, oldEnInt, newEnInt, sigmaDoubleDotD
      real*8 F(3,3), FOld(3,3), sigma(3,3)
      real*8 FDot(3,3), FInv(3,3), L(3,3), D(3,3)

      FDot = (F - FOld)/dt

      call M3INV(F,FInv)

      L = matmul(FDot,FInv)
      D = 0.5*(L + transpose(L))

      call DOTPM(sigma,D,sigmaDoubleDotD)

      newEnInt = oldEnInt + (sigmaDoubleDotD/dens)*dt

      end