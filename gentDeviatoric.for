C     requires UR.for utility subroutine

      subroutine gentDeviatoric(F,Im,amu,sigma,enerGentDeviatoric)

      implicit none

      real*8 Im, Jm, amu, J, I1, I1_bar, enerGentDeviatoric
      real*8 F(3,3), F_bar(3,3), B(3,3), B_bar(3,3), sigma(3,3)
      real*8 EYE(3,3)

      Jm = Im - 3.0D0

      call ONEM(EYE)
      call MDET(F,J)

      B = matmul(F,transpose(F))

      F_bar = J**(-1.0D0/3.0D0) * F
      B_bar = matmul(F_bar,transpose(F_bar))

      call TRACEM(B_bar,I1_bar)

      sigma = J**(-5.0D0/3.0D0) * (amu*Jm)/(Jm - I1_bar + 3.0D0) * B
     1 - J**(-1.0D0) * (1.0D0/3.0D0) * I1_bar
     2 * (amu*Jm)/(Jm - I1_bar + 3.0D0) * EYE

      enerGentDeviatoric = 
     1(-amu*Jm/2.0D0)*LOG(1.0D0-((I1_bar-3.0D0)/Jm))

      end