!===========================================================
! Declaration of global variables
!***********************************************************      
      module param
        implicit none
        !===========================================================
        !      grid size
        integer,parameter :: mref1=1
        integer,parameter :: mref2=1
        integer,parameter :: mref3=1
        integer :: m1,m2,m3

        integer :: m2m,m3m,m2mh,m1m
        integer :: m1r,m2r,m3r
        integer :: m1mr,m2mr,m3mr,m2mhr

       !==============================================================
       !      inital condition
        integer,parameter :: tag_ini = 1
        integer,parameter :: Nitf = 1

       !==============================================================
       !      bou.in
        integer   :: n1, n2, n3, nsst, nwrit, nread
        integer   :: ntst, ireset
        real      :: tpin, trestart, tmax
        real      :: alx3, str3, rext1, rext2
        integer   :: istr3, lmax
        real      :: Ra,Sc,Lambdaco2,Lambdavap,InOutvel,InOutzlen,Inzlevel,Outzlevel,kernel_width_space,kernel_width_time
        integer   :: idtv
        real      :: dtmax, resid, cflmax, dt,dt_o, cfllim, dts
        real      :: tframe, tslice
        integer   :: mov_zcut_k

        integer   :: ubctop, ubcbot

       !==============================================================
       !      part.in
        integer   :: imlsfor,imlssca,pread
        real      :: sclf1,sclf2,xpos1,xpos2, objinitheight
        character*50  objfx1,objfx2

       !==============================================================
       !     simulation parameters
        real      :: time
        real      :: dsaltop, dsalbot, denstop, densbot, co2top,co2bot
        real :: dx2,dx3,dx1
        real :: dx2q,dx3q,dx1q

        real :: dx2r,dx3r,dx1r
        real :: dx2qr,dx3qr,dx1qr
    
        real, allocatable, dimension(:) :: xc,xm
        real, allocatable, dimension(:) :: yc,ym
        real, allocatable, dimension(:) :: zc,zm,g3rc,g3rm
    
        real, allocatable, dimension(:) :: xcr,xmr
        real, allocatable, dimension(:) :: ycr,ymr
        real, allocatable, dimension(:) :: zcr,zmr,g3rcr,g3rmr

       !==============================================================
       !     quantities for derivatives
        real, allocatable, dimension(:) :: udx3c,udx3m
        real, allocatable, dimension(:) :: udx3cr,udx3mr

       !==============================================================
       !     grid indices
        integer, allocatable, dimension(:) :: jmv,jpv
        integer, allocatable, dimension(:) :: imv,ipv
        integer, allocatable, dimension(:) :: jmhv
        integer, allocatable, dimension(:) :: kmc,kpc,kmv,kpv,kup,kum

        integer, allocatable, dimension(:) :: imvr,ipvr
        integer, allocatable, dimension(:) :: jmvr,jpvr
        integer, allocatable, dimension(:) :: kmvr,kpvr
  
       !==============================================================
       !     metric coefficients
        real, allocatable, dimension(:) :: ap3j,ac3j,am3j
        real, allocatable, dimension(:) :: ap3ck,ac3ck,am3ck
        real, allocatable, dimension(:) :: ap3sk,ac3sk,am3sk
        real, allocatable, dimension(:) :: ap3ssk,ac3ssk,am3ssk   
        real, allocatable, dimension(:) :: ap3sskr,ac3sskr,am3sskr

       !==============================================================
       !     variables for FFTW and Poisson solver
        real, dimension(13) :: ifx1
        real, allocatable, dimension(:) :: trigx1
        real, allocatable, dimension(:) :: ak2,ap
        real, allocatable, dimension(:) :: ak1,ao
        real, allocatable, dimension(:) :: amphk,acphk,apphk
        integer*8 :: fwd_plan,bck_plan
        
       !==============================================================
       !     other variables
        integer  :: n2m, n3m, n1m
        integer :: n1r,n2r,n3r,n1mr,n2mr,n3mr
        real :: byct, bycs, bycco2, nu,Dw
        real :: pi
        real :: al,ga,ro
        real :: beta
        integer :: ntime
        integer, parameter:: ndv=3
        real, dimension(1:ndv) :: vmax
        real, dimension(1:3) :: gam,rom,alm
        real :: usref1, usref2, usref3              
      end module param
      

!===========================================================
! 3D arrays, dynamically allocated by each process
!***********************************************************
      module local_arrays
      use param
        implicit none
        real,allocatable,dimension(:,:,:) :: q1,q2,q3,dens,dsal,co2
        real,allocatable,dimension(:,:,:) :: hro,hsal,rhs,rhsr,hco2
        real,allocatable,dimension(:,:,:) :: ru1,ru2,ru3,ruro,rusal,ruco2
        real,allocatable,dimension(:,:,:) :: pr,qcap,dph,dq
        real,allocatable,dimension(:,:,:) :: isbody1,isbody2 !CS Heated body
        real,allocatable,dimension(:,:,:) :: exhale_signal !CS Breath signal
      end module local_arrays
    
      module mpih
        implicit none
        include 'mpif.h'
        integer :: myid, numtasks, numthreads, ierr
        integer, parameter :: master=0
        integer, parameter :: lvlhalo=1
        integer :: MDP = MPI_DOUBLE_PRECISION
        integer :: STATUS(MPI_STATUS_SIZE,4)
        integer :: req(1:4)
        integer(kind=MPI_OFFSET_KIND) :: disp, offset
      end module mpih
      
!===========================================================
! mpi param
!***********************************************************
      module mpi_param
        implicit none
        integer :: istart,iend, jstart,jend, kstart,kend
        integer :: jstartr,jendr, kstartr,kendr
        integer :: jstartp,jendp
        integer :: dj,dk,mydata,mydatam
        integer :: djp,djr,dkr
        integer, allocatable, dimension(:) :: offsetj,offsetk
        integer, allocatable, dimension(:) :: offsetjr,offsetkr
        integer, allocatable, dimension(:) :: offsetjp
        integer, allocatable, dimension(:) :: countj,countk
        integer, allocatable, dimension(:) :: countjr,countkr
        integer, allocatable, dimension(:) :: countjp
        integer, allocatable, dimension(:) :: countf
        integer(8), allocatable, dimension(:) :: offsetf 
      end module mpi_param

!===========================================================
! multiresolution
!***********************************************************
      module mgrd_arrays
        use param
        implicit none
        integer, parameter :: mrefa=mref1*mref2*mref3
        integer, allocatable, dimension(:) :: irangs,jrangs,krangs
        real, allocatable, dimension(:,:) :: cxq1, cxq2, cxq3, cxrs
        real, allocatable, dimension(:,:) :: cyq1, cyq2, cyq3, cyrs
        real, allocatable, dimension(:,:) :: czq1, czq2, czq3, czrs
        real, allocatable,dimension(:,:,:) :: q1lr,q2lr,q3lr,dsalc
      end module mgrd_arrays

!===========================================================
! walls BCs
!***********************************************************
      module walls_vars
        use param
        implicit none
        real, allocatable, dimension(:,:) :: q1ybn,q2ybn,q3ybn
        real, allocatable, dimension(:,:) :: densybn,dsalybn,co2ybn
        real, allocatable, dimension(:,:) :: q1ybs,q2ybs,q3ybs
        real, allocatable, dimension(:,:) :: densybs,dsalybs,co2ybs
      end module walls_vars

!===========================================================
! outflow BCs
!***********************************************************
      module outflow_vars
        use param
        implicit none
        real, allocatable, dimension(:,:) :: q1be,dq1be
        real, allocatable, dimension(:,:) :: q2be,dq2be
        real, allocatable, dimension(:,:) :: q3be,dq3be
        real, allocatable, dimension(:,:) :: dq1beo
        real, allocatable, dimension(:,:) :: dq2beo
        real, allocatable, dimension(:,:) :: dq3beo
        real, allocatable, dimension(:,:) :: dsalbe,ddsalbe
        real, allocatable, dimension(:,:) :: ddsalbeo
        real, allocatable, dimension(:,:) :: densbe,ddensbe
        real, allocatable, dimension(:,:) :: ddensbeo
        real, allocatable, dimension(:,:) :: co2be,dco2be
        real, allocatable, dimension(:,:) :: dco2beo
      end module outflow_vars

!===========================================================
! inflow BCs
!***********************************************************
      module inflow_vars
        use param
        implicit none
        real, allocatable, dimension(:,:) :: q1bw,dq1bw
        real, allocatable, dimension(:,:) :: q2bw,dq2bw
        real, allocatable, dimension(:,:) :: q3bw,dq3bw
        real, allocatable, dimension(:,:) :: dsalbw,ddsalbw
        real, allocatable, dimension(:,:) :: densbw,ddensbw
        real, allocatable, dimension(:,:) :: co2bw,dco2bw
      end module inflow_vars

