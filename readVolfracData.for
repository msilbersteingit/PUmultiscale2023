      subroutine readVolfracData(fileName,nLines,data)
        implicit none
        integer nLines, i
        real*8 data(nLines,5)
        character(*) :: fileName

        open (105, file = fileName)
        do i=1,nLines
        read (105,*) data(i,1), data(i,2), data(i,3),
     1               data(i,4), data(i,5)
        enddo
        close (105)

      end subroutine readVolfracData