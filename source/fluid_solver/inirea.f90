      subroutine inirea
      use mpih
      use mpi_param, only: kstart,kend,kstartr,kendr
      use local_arrays, only: dens,q2,q3,q1,dsal,isbody1,isbody2,co2
      use param
      IMPLICIT NONE
      real :: dummyr 
      character(50) :: filnam
      integer :: i,j
      integer :: n1o,n2o,n3o,istr3o,str3o,mref1o,mref2o,mref3o
      integer, parameter :: ghosts = 4
      
      
      ! Reading old grid and time information by rank0
      if (myid .eq. 0) then
        filnam = 'continua_grid.h5'
        call HdfStart
        call HdfSerialReadIntScalar('n1',filnam,n1o)
        call HdfSerialReadIntScalar('n2',filnam,n2o)
        call HdfSerialReadIntScalar('n3',filnam,n3o)
        call HdfSerialReadRealScalar('time',filnam,time)
        call HdfSerialReadIntScalar('istr3',filnam,istr3o)
        call HdfSerialReadIntScalar('str3',filnam,str3o)
        call HdfSerialReadIntScalar('mref1',filnam,mref1o)
        call HdfSerialReadIntScalar('mref2',filnam,mref2o)
        call HdfSerialReadIntScalar('mref3',filnam,mref3o)
        call HdfClose
        write(*,'(5x,a)') 'old grid info read'
      endif

      call MPI_BCAST(n1o,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(n2o,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(n3o,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(istr3o,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(str3o,1,MDP,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(time,1,MDP,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(mref1o,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(mref2o,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(mref3o,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 
      
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      !  One to one HDF read
      call mpi_read_continua(n1,n2,n3,kstart,kend,1,q1)
      call mpi_read_continua(n1,n2,n3,kstart,kend,2,q2)
      call mpi_read_continua(n1,n2,n3,kstart,kend,3,q3)
      call mpi_read_continua(n1,n2,n3,kstart,kend,4,dens)
      call mpi_read_continua(n1r,n2r,n3r,kstartr,kendr,5,dsal)
      call mpi_read_continua(n1r,n2r,n3r,kstartr,kendr,8,co2)
      call mpi_read_continua(n1,n2,n3,kstart,kend,9,isbody1)
      call mpi_read_continua(n1,n2,n3,kstart,kend,10,isbody2)

      if (ireset.eq.1) then                                             
        time=0.d0
      endif                                                             

      return                                                            
      end                                                                                                                                  
