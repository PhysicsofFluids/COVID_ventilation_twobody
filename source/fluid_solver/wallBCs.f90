!===========================================================
! memory allocations
!***********************************************************
      subroutine wallBCs_mem_alloc
      use mpih
      use param
      use mpi_param
      use walls_vars
      implicit none

      allocate( q1ybn(1:n1,kstart-lvlhalo:kend+lvlhalo), q2ybn(1:n1,kstart-lvlhalo:kend+lvlhalo), q3ybn(1:n1,kstart-lvlhalo:kend+lvlhalo) )
      allocate( densybn(1:n1,kstart-lvlhalo:kend+lvlhalo),dsalybn(1:n1r,kstartr-lvlhalo:kendr+lvlhalo),co2ybn(1:n1,kstart-lvlhalo:kend+lvlhalo) )
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      q1ybn(:,:)=0.0d0;q2ybn(:,:)=0.0d0;q3ybn(:,:)=0.0d0
      densybn(:,:)=0.0d0;dsalybn(:,:)=0.0d0;co2ybn(:,:)=0.0d0

      allocate( q1ybs(1:n1,kstart-lvlhalo:kend+lvlhalo), q2ybs(1:n1,kstart-lvlhalo:kend+lvlhalo), q3ybs(1:n1,kstart-lvlhalo:kend+lvlhalo) )
      allocate( densybs(1:n1,kstart-lvlhalo:kend+lvlhalo), dsalybs(1:n1r,kstartr-lvlhalo:kendr+lvlhalo), co2ybs(1:n1,kstart-lvlhalo:kend+lvlhalo) )
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      q1ybs(:,:)=0.0d0;q2ybs(:,:)=0.0d0;q3ybs(:,:)=0.0d0
      densybs(:,:)=0.0d0;dsalybs(:,:)=0.0d0;co2ybs(:,:)=0.0d0

      return
      end

