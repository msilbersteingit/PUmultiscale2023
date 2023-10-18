      subroutine initVolfrac(x,y,z,nLines,vfData,vf)
        implicit none
        integer i, nLines
        real*8 vfData(nLines,5), x, y, z, vf
        real*8 distMin, dist

C        do i = 1,nLines
C            if (abs(x-vfData(i,2)) .LE. 0.0000001 .and.  
C     1         abs(y-vfData(i,3)) .LE. 0.0000001 .and.
C     2         abs(z-vfData(i,4)) .LE. 0.0000001) then
C               vf = vfData(i,5)
C               exit
C            end if
C        end do

        distMin = ((x-vfData(1,2))**2 + (y-vfData(1,3))**2
     1              + (z-vfData(1,4))**2)**(1.0D0/2.0D0)
        
        vf = vfData(1,5)


        do i = 2,nLines
          dist = ((x-vfData(i,2))**2 + (y-vfData(i,3))**2
     1              + (z-vfData(i,4))**2)**(1.0D0/2.0D0)

          if (dist .LT. distMin) then
            distMin = dist
            vf = vfData(i,5)
          end if
        end do

      end subroutine initVolfrac