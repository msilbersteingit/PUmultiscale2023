C     requires UR.for utility subroutine
      subroutine volumetric_hard(F,k,rho_h0,rho_h1,xi,
     1                       sigma,enerVolumetric)
      implicit none

      real*8 k_h0, k_h1, k
      real*8 F(3,3), sigma(3,3), EYE(3,3)
      real*8 J, J_el, J_dens, xi, rho_h0, rho_h1
      real*8 enerVolumetric

      call ONEM(EYE)
      call MDET(F,J)

      J_dens = rho_h0/ ( (rho_h1-rho_h0)*(1.0D0-xi) + rho_h0 ) 

      J_el = J/J_dens

      sigma = k*(J_el - 1.0D0)*EYE

      enerVolumetric = (k/2.0D0)* ((J_el - 1.0D0)**2.0D0)

      end subroutine volumetric_hard
