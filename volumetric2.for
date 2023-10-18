C     requires UR.for utility subroutine
      subroutine volumetric2(F,k,sigma,enerVolumetric)

      implicit none

      real*8 k
      real*8 F(3,3), sigma(3,3), J, EYE(3,3)
      real*8 enerVolumetric

      call ONEM(EYE)
      call MDET(F,J)

      sigma = k*(J - 1.0D0)*EYE

      enerVolumetric = (k/2.0D0)* ((J - 1.0D0)**2.0D0)

      end subroutine volumetric2
