      subroutine openfi
      use param
      use mpih
      implicit none

      IF(myid.eq.0)then

        open(95, file='./data/vmax.out',status='unknown',access='sequential',position='append')
        open(96, file='./cfl.out',status='unknown',access='sequential',position='append')
        open(97, file='./data/nusse_walls.out',status='unknown',access='sequential',position='append')
        open(98, file='./data/Savg.out',status='unknown',access='sequential',position='append')
        open(99, file='./data/co2avg.out',status='unknown',access='sequential',position='append')
        open(100, file='./data/Tavg.out',status='unknown',access='sequential',position='append')
      ! reset the time history
      if(ireset.eq.1 .or. nread.eq.0)then
          rewind(95)
          rewind(96)
          rewind(97)
          rewind(98)
          rewind(99)
          rewind(100)
      endif

      ENDIF

      return
      end   
      
!==============================================

      subroutine closefi
      use mpih
      implicit none
      
      if(myid.eq.0)then

      !EP  in vmaxv.F
      close(95)
      close(96)
      close(97)
      close(98)
      close(99)
      close(100)

      endif
     
      return      
      end
