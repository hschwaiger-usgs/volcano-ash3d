!******************************************************************************

      subroutine write_3D_Binary(cio,ns)

!     Subroutine that writes out 3-D arrays in binary format

      use precis_param

      use mesh,          only : &
         nxmax,nymax,nzmax,ts1

      use solution,      only : &
         concen_pd

      use Output_Vars,   only : &
         DepositThickness

      implicit none

      character(len=13) ,intent(in) :: cio
      integer           ,intent(in) :: ns

      integer :: i,j,k,n
      real(kind=op), dimension(:,:,:),allocatable :: ashcon

      ! Write out data in raw binary form
      allocate(ashcon(nxmax,nymax,nzmax))
      ashcon = 0.0_op
      do n=1,ns
        ashcon(1:nxmax,1:nymax,1:nzmax) =  &
         ashcon(1:nxmax,1:nymax,1:nzmax) + &
         real(concen_pd(1:nxmax,1:nymax,1:nzmax,n,ts1),kind=op)
      enddo

      ! 3-D total tephra concentration
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

      ! 2-D Deposit thickness
      if(op.eq.4)then
        open(unit=21,file='2d_tephra_depo_'//cio//'.raw', &
          status='replace', &
          access='direct',recl=4*nxmax*nymax)
      else
        open(unit=20,file='2d_tephra_depo_'//cio//'.raw', &
          status='replace', &
          access='direct',recl=8*nxmax*nymax)
      endif
      write(21,rec=1)((DepositThickness(i,j),i=1,nxmax),j=1,nymax)
      close(21)

      end subroutine write_3D_Binary


