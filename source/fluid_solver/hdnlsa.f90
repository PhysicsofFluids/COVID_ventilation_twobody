      subroutine hdnlsa
      use param
      use local_arrays, only: hsal,dsal
      use mgrd_arrays, only: q1lr,q2lr,q3lr
      use mpi_param, only: kstartr,kendr
      use walls_vars
      implicit none
      integer :: jc,kc,ic,js
      integer :: kpp,kmm,kp,jp,jmm,jpp,ip,imm,ipp
      real    :: h31,h32,h33,udx2r,udx1r
      real    :: dmpcoef, ydmp, ys
      real    :: fp1d2,fm1d2,phip1d2,phim1d2

      udx1r=dx1r*0.5d0
      udx2r=dx2r*0.5d0
      do kc=kstartr,kendr
        kmm=kc-1
        kpp=kc+1
        kp=kc+1
        do jc=1,n2mr
          jp=jc+1
          jmm=jc-1
          jpp=jc+1
          do ic=1,n1mr
            ip=ic+1
            ipp=ic+1
            imm=ic-1
            if(ic.eq.1) then
               h31=( q1lr(ip,jc,kc)*(dsal(ipp,jc,kc)+dsal(ic,jc,kc)))*udx1r
            elseif(ic.eq.n1mr) then
               h31=(-q1lr(ic,jc,kc)*(dsal(ic,jc,kc)+dsal(imm,jc,kc)))*udx1r
            else
               h31=( q1lr(ip,jc,kc)*(dsal(ipp,jc,kc)+dsal(ic,jc,kc)) &
        &           -q1lr(ic,jc,kc)*(dsal(ic,jc,kc)+dsal(imm,jc,kc)))*udx1r
            endif

            if(jc.eq.1) then
               h32=( q2lr(ic,jp,kc)*(dsal(ic,jpp,kc)+dsal(ic,jc,kc)) &
        &           -q2ybs(ic,kc)*(dsal(ic,jc,kc)+dsalybs(ic,kc)))*udx2r
            elseif(jc.eq.n2mr) then
               h32=( q2ybn(ic,kc)*(dsal(ic,jc,kc)+dsalybn(ic,kc)) &
        &           -q2lr(ic,jc,kc)*(dsal(ic,jc,kc)+dsal(ic,jmm,kc)))*udx2r
            else
               h32=( q2lr(ic,jp,kc)*(dsal(ic,jpp,kc)+dsal(ic,jc,kc)) &
        &           -q2lr(ic,jc,kc)*(dsal(ic,jc,kc)+dsal(ic,jmm,kc)))*udx2r
            endif

            if(kc.eq.1) then
               h33=( (dsal(ic,jc,kpp)*g3rmr(kc)+dsal(ic,jc,kc)*g3rmr(kpp))&
        &              /(g3rmr(kc)+g3rmr(kpp))*q3lr(ic,jc,kp))*udx3mr(kc)
            elseif(kc.eq.n3m) then
               h33=(-(dsal(ic,jc,kc)*g3rmr(kmm)+dsal(ic,jc,kmm)*g3rmr(kc))&
        &              /(g3rmr(kc)+g3rmr(kmm))*q3lr(ic,jc,kc)&
        &          )*udx3mr(kc)
            else
               h33=( (dsal(ic,jc,kpp)*g3rmr(kc)+dsal(ic,jc,kc)*g3rmr(kpp))&
        &              /(g3rmr(kc)+g3rmr(kpp))*q3lr(ic,jc,kp)&
        &           -(dsal(ic,jc,kc)*g3rmr(kmm)+dsal(ic,jc,kmm)*g3rmr(kc))&
        &              /(g3rmr(kc)+g3rmr(kmm))*q3lr(ic,jc,kc)&
        &          )*udx3mr(kc)
            endif

            hsal(ic,jc,kc)= -(h31+h32+h33)
     
          enddo
        enddo
      enddo

      return
      end

