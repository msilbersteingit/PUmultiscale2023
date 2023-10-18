C     Requires UR.for utility subroutine

      subroutine elDeviatoric(F,amu,sigma_bar,enerElDeviatoric)

      implicit none

      real*8 amu, J
      real*8 F(3,3), F_bar(3,3), B(3,3), B_bar(3,3), sigma_bar(3,3)
      real*8 B_barInv(3,3), e_bar(3,3)
      real*8 EYE(3,3)
      real*8 enerElDeviatoric, sigma_barDDe_bar

      call ONEM(EYE)
      call MDET(F,J)

      F_bar = J**(-1.0D0/3.0D0) * F
      B_bar = matmul(F_bar,transpose(F_bar))
      
      call M3INV(B_bar,B_barInv)
      e_bar = 0.5D0*(EYE - B_barInv)

      sigma_bar = 2.0D0*amu*e_bar

      call DOTPM(sigma_bar,e_bar,sigma_barDDe_bar)
      
      enerElDeviatoric = 0.5D0*sigma_barDDe_bar


      end subroutine elDeviatoric

