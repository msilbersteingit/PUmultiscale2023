      subroutine numLinesVolfracData(fileName,nLines)
        implicit none
        integer nLines, io
        character(*) :: fileName

        nLines = 0 
        open (105, file = fileName)
        DO 
          READ (105,*, END=10) 
           nLines = nLines + 1 
        END DO
   10   CLOSE (105) 






      end subroutine numLinesVolfracData

