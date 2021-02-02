      subroutine soltek
      use param
      use local_arrays, only : dens,rhs
      use mpi_param
      use mpih
      implicit none
      integer :: jc,kc,ic,info,ipkv(m3mr)
      real :: betadx,fkl(m3mr,m1mr)
      real :: amkT(m3mr-1),ackT(m3mr),apkT(m3mr-1)
      real, allocatable, dimension(:,:,:) :: rhst

      allocate(rhst(1:n3m,1:n1m,jstart:jend))

      call PackZ_UnpackR(rhs,rhst)

      betadx=0.5d0*al*dts*Dw
     
      do kc=1,n3m
        ackT(kc)=(1.d0-ac3ssk(kc)*betadx)
      enddo

      amkT(1:n3m-1)=-betadx*am3ssk(2:n3m)
      apkT(1:n3m-1)=-betadx*ap3ssk(1:n3m-1)

      call ddttrfb(n3m,amkT,ackT,apkT,info)

      do jc=jstart,jend
        do ic=1,n1m
          do kc=1,n3m
            fkl(kc,ic)=rhst(kc,ic,jc)
          end do
        enddo

       call ddttrsb('N',n3m,n1m,amkT,ackT,apkT,fkl,n3m,info)
          
        do ic=1,n1m
          do kc=1,n3m
            rhst(kc,ic,jc)= fkl(kc,ic)
          end do
        enddo
      end do

      call PackR_UnpackZ(rhst,rhs)

      do kc=kstart,kend
        do jc=1,n2m
          do ic=1,n1m
            dens(ic,jc,kc) = dens(ic,jc,kc) + rhs(ic,jc,kc)
          enddo
        enddo
      enddo

      if(allocated(rhst)) deallocate(rhst)

      return
      end
