!******************************************************************************

      subroutine write_3D_Binary(cio,nx,ny,nz,ashcon_tot)

!     Subroutine that writes out 3-D arrays in binary format

      use precis_param

      implicit none

      character(len=13) ,intent(in) :: cio
      integer           ,intent(in) :: nx
      integer           ,intent(in) :: ny
      integer           ,intent(in) :: nz
      real(kind=op)     ,intent(in) :: ashcon_tot(nx,ny,nz)

      integer :: i,j,k

      ! Write out data in raw binary form

      ! 3-D total tephra concentration
      if(op.eq.4)then
        open(unit=20,file='3d_tephra_fall_'//cio//'.raw', &
          status='replace', &
          access='direct',recl=4*nx*ny*nz)
      else
        open(unit=20,file='3d_tephra_fall_'//cio//'.raw', &
          status='replace', &
          access='direct',recl=8*nx*ny*nz)
      endif
      write(20,rec=1)(((ashcon_tot(i,j,k),i=1,nx),j=1,ny),k=1,nz)
      close(20)

      end subroutine write_3D_Binary

!******************************************************************************

      subroutine write_2D_Binary(nx,ny,OutVar,VarMask,Fill_Value,filename_root)

!     Subroutine that writes out 2-D arrays in binary format

      use precis_param

      use io_units

      use io_data,       only : &
         isFinal_TS,iout3d,WriteTimes

      implicit none

      integer          ,intent(in) :: nx
      integer          ,intent(in) :: ny
      real(kind=ip)    ,intent(in) :: OutVar(nx,ny)
      logical          ,intent(in) :: VarMask(nx,ny)
      character(len=6) ,intent(in) :: Fill_Value
      character(len=20),intent(in) :: filename_root

      real(kind=op)  :: OVar(nx,ny)
      real(kind=op)  :: FValue
      integer :: fid
      integer :: i,j
      character (len=9)  :: cio
      character(len=50)  :: filename_out

      read(Fill_Value,*)FValue

      fid = 30

      if(isFinal_TS)then
        cio='____final'
      else
        if (WriteTimes(iout3d).lt.10.0_ip) then
          write(cio,1) WriteTimes(iout3d)
1         format('00',f4.2,'hrs')
        elseif (WriteTimes(iout3d).lt.100.0_ip) then
          write(cio,2) WriteTimes(iout3d)
2         format('0',f5.2,'hrs')
        else
          write(cio,3) WriteTimes(iout3d)
3         format(f6.2,'hrs')
        endif
      endif
      if(INDEX(filename_root,'ArrivalTime').gt.0)then
          ! For the special cases of DepositArrivalTime.dat and
          ! CloudArrivalTime.dat
        write(filename_out,*)trim(adjustl(filename_root)),'.raw'
      else
        write(filename_out,*)trim(adjustl(filename_root)),cio,'.raw'
      endif

      ! Apply threshold mask to output variable
      do i=1,nx
        do j=1,ny
          if(VarMask(i,j))then
            OVar(i,j) = real(OutVar(i,j),kind=op)
          else
            OVar(i,j) = FValue
          endif
        enddo
      enddo

      if(op.eq.4)then
        open(unit=21,file=trim(adjustl(filename_out)), &
          status='replace', &
          access='direct',recl=4*nx*ny)
      else
       open(unit=20,file=trim(adjustl(filename_out)), &
          status='replace', &
          access='direct',recl=8*nx*ny)
      endif
      write(21,rec=1)((OVar(i,j),i=1,nx),j=1,ny)
      close(21)

      end subroutine write_2D_Binary
