      subroutine volavg
      use param
      use local_arrays, only: dsal, co2,isbody1,isbody2,dens
      use mpih
      implicit none
      integer :: j,i,k,kstart,kend,kstartr,kendr
      real :: sumbody1, sumbody2, sumsal, sumco2, sumdens
      real :: my_sumbody1, my_sumbody2, my_sumsal, my_sumco2, my_sumdens
      real :: vol, volr
      real :: Savg, co2avg, Tavg

      vol = dble(n1m*n2m*n3m)      
      volr = dble(n1mr*n2mr*n3mr) 

      my_sumsal=sum(dsal(:,:,:))
      my_sumbody1=sum(isbody1(:,:,:))
      my_sumbody2=sum(isbody2(:,:,:))
      my_sumco2=sum(co2(:,:,:))
      my_sumdens=sum(dens(:,:,:))
      call MPI_REDUCE(my_sumsal,sumsal,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(my_sumco2,sumco2,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(my_sumdens,sumdens,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(my_sumbody1,sumbody1,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(my_sumbody2,sumbody2,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      sumsal = sumsal   - 0.0*sumbody1 - 0.0*sumbody2
      Savg = sumsal/(volr-sumbody1-sumbody2)
      sumco2 = sumco2   - 0.0*sumbody1 - 0.0*sumbody2
      co2avg = sumco2/(vol-sumbody1-sumbody2)
      sumdens = sumdens - 1.0*sumbody1 - 1.0*sumbody2
      Tavg = sumdens/(vol-sumbody1-sumbody2)

      if(myid.eq.0)then
        write(98,547) time,Savg
      endif
547   format(1x,f10.4,1x,ES20.8)

      if(myid.eq.0)then
        write(99,548) time,co2avg
      endif
548   format(1x,f10.4,1x,ES20.8)

      if(myid.eq.0)then
        write(100,549) time,Tavg
      endif
549   format(1x,f10.4,1x,ES20.8)

      return         
      end               
