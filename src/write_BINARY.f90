!******************************************************************************

      subroutine write_3D_Binary(cio,ns)

!     Subroutine that writes out 3-D arrays in binary format

      use precis_param

      use mesh,          only : &
         nxmax,nymax,nzmax,ts1,dz_vec_pd

      use solution,      only : &
         concen_pd

      implicit none

      character(len=13) ,intent(in) :: cio
      integer           ,intent(in) :: ns

      integer :: i,j,k,n
      real(kind=op), dimension(:,:,:),allocatable :: ashcon
      real(kind=op), dimension(:,:)  ,allocatable :: depocon

      ! Write out data in raw binary form
      allocate(ashcon(nxmax,nymax,nzmax))
      allocate(depocon(nxmax,nymax))
      ashcon = 0.0_op
      do n=1,ns
        ashcon(1:nxmax,1:nymax,1:nzmax) =  &
         ashcon(1:nxmax,1:nymax,1:nzmax) + &
         real(concen_pd(1:nxmax,1:nymax,1:nzmax,n,ts1),kind=op)
      enddo
        ! depocon
        depocon = 0.0_op
      do n=1,ns
        depocon(1:nxmax,1:nymax) =  &
         depocon(1:nxmax,1:nymax) + &
         real(concen_pd(1:nxmax,1:nymax,0,n,ts1),kind=op) * &
         real(dz_vec_pd(0),kind=op)*1.0e-6_op
      enddo
      if(op.eq.4)then
        open(unit=20,file='3d_tephra_fall_'//cio//'.raw', &
          status='replace', &
          access='direct',recl=4*nxmax*nymax*nzmax)
      else
        open(unit=20,file='3d_tephra_fall_'//cio//'.raw', &
          status='replace', &
          access='direct',recl=8*nxmax*nymax*nzmax)
      endif
      write(20,rec=1)(((ashcon(i,j,k),i=1,nxmax),j=1,nymax),k=1,nzmax)
      close(20)
      deallocate(ashcon)
      deallocate(depocon)

      end subroutine write_3D_Binary


