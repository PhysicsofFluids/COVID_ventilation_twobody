      subroutine tschem
      use param
      use local_arrays
      use mgrd_arrays
      use mpih
      use mpi_param, only: kstart,kend,kstartr,kendr
      implicit none
      integer :: ns,inp,ntr
      integer :: j,k,i,mstep

      dts = dt/dble(lmax)
!====== Begin RK scheme ======!
      do ns=1,nsst
        al=alm(ns)
        ga=gam(ns)
        ro=rom(ns)

        !-- Set var to zero
        dph(:,:,:)=0.0d0
        dq(:,:,:)=0.0d0
        qcap(:,:,:)=0.0d0

        !-- Outflow BCs update through wave eq.
        call calcOutflowBCs
 
        !-- update Outflow BCs
        call updateOutflowBCs
        call updateInflowBCs

        !-- Explicit terms
        call hdnlq1
        call hdnlq2
        call hdnlq3
        call hdnlsa
        call hdnlte
        call hdnlco2

        !-- Implicit terms
        call invtrq1
        call invtrq2
        call invtrq3
        call update_upper_ghost(n1,n2,q3)
        call invtrsa
        call update_both_ghosts(n1r,n2r,dsal,kstartr,kendr)
        call invtrte
        call update_both_ghosts(n1,n2,dens,kstart,kend)
        call invtrco2
        call update_both_ghosts(n1,n2,co2,kstart,kend)

        !-- velbc
        call velbc

        !-- add body obj
        call addbodyobj

        !-- Swap before pressure solver !-- Important
        !call update_both_ghosts(n1,n2,q1,kstart,kend) !not needed for divg
        !call update_both_ghosts(n1,n2,q2,kstart,kend) !not needed for divg
        !call update_both_ghosts(n1,n2,q3,kstart,kend) !not needed for divg
        !call update_both_ghosts(n1r,n2r,dsal,kstartr,kendr) !not needed for multiresol
        !call update_both_ghosts(n1,n2,dens,kstart,kend)
        !call update_both_ghosts(n1,n2,co2,kstart,kend)

        !-- Velocity pressure correction
        dph(:,:,:)=0.0d0
        call update_upper_ghost(n1,n2,q3)
        call divg
        call phcalc
        call update_both_ghosts(n1m,n2m+2,dph,kstart,kend)
        call updvp                 ! SOLENOIDAL VEL FIELD
        call update_both_ghosts(n1,n2,q1,kstart,kend)
        call update_both_ghosts(n1,n2,q2,kstart,kend)
        call update_both_ghosts(n1,n2,q3,kstart,kend)
        call prcalc                         !! PRESSURE FIELD
        call update_both_ghosts(n1,n2,pr,kstart,kend)

        !==== Multiresol info exchange ====!
        call mgrd_velitp
        call mgrd_dsalc
        call update_both_ghosts(n1,n2,dsalc,kstart,kend)
      enddo
!====== End RK scheme ======!

      return                                                            
      end                                                               
