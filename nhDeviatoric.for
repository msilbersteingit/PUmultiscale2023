C     requires UR.for utility subroutine

      subroutine nhDeviatoric(F,amu,sigma,enerNhDeviatoric)

      implicit none

      real*8 amu, J, I1, I1_bar
      real*8 F(3,3), F_bar(3,3), B(3,3), B_bar(3,3), sigma(3,3)
      real*8 EYE(3,3)
      real*8 enerNhDeviatoric

      call ONEM(EYE)
      call MDET(F,J)

      B = matmul(F,transpose(F))

      F_bar = J**(-1.0D0/3.0D0) * F
      B_bar = matmul(F_bar,transpose(F_bar))

      call TRACEM(B_bar,I1_bar)


C     amu = 2C

      sigma = J**(-5.0D0/3.0D0) * amu * B
     1 - J**(-1.0D0) * (1.0D0/3.0D0) * I1_bar
     2 * amu * EYE

      enerNhDeviatoric = (amu/2.0D0)*(I1_bar - 3.0D0)


      end