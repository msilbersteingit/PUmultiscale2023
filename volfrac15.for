C     Requires UR.for utility subroutine

      subroutine volfrac15(F,gamaDot,xi_h0Old,
     1                     xi_h0New,dt,rho_h0,rho_h1,param1,
     2                     param2,J,enerDev_h0,enerDev_h1,J_dens)

      implicit none

      real*8 F(3,3)
      real*8 gamaDot
      real*8 xi_h0Old, xi_h0New, dt
      real*8 rho_h0, rho_h1
      real*8 param1, param2
      real*8 xi_h0Dot
      real*8 J, J_dens
      real*8 enerDev_h0, enerDev_h1

      xi_h0Dot = param1*(gamaDot)**param2

      xi_h0New = xi_h0Old + xi_h0Dot*dt

      call MDET(F,J)
      J_dens = rho_h0/ ( (rho_h1-rho_h0)*(1.0D0-xi_h0New) +  rho_h0 ) 

      if (J_dens .GT. J) then
            xi_h0New = xi_h0Old
      else if (enerDev_h1 .GT. enerDev_h0) then
            xi_h0New = xi_h0Old
      end if


C     keep this incase ALE rescales vf
      if (xi_h0Old .GT. 1.0D0) then
            xi_h0New = 1.0D0
      else if (xi_h0Old .LE. 0.0D0) then
            xi_h0New = 0.0D0
      end if

      if (xi_h0New .LE. 0.0D0) then
            xi_h0New = 0.0D0
      else if (xi_h0New .GT. 1.0D0) then
            xi_h0New = 1.0D0
      end if



      end subroutine volfrac15