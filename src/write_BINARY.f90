!##############################################################################
!
!  Ash3d_Binary_IO module
!
!  This module manages all output to kml files
!
!      subroutine write_2D_Binary
!      subroutine read_2D_Binary(filename)
!      subroutine write_3D_Binary
!      subroutine read_3D_Binary(filename)
!
!##############################################################################

      module Ash3d_Binary_IO

      use precis_param

      use io_units

      implicit none

        ! Set everything to private by default
      private

        ! Publicly available subroutines/functions
      public deallocate_Binary,  &
             write_2D_Binary,    &
             read_2D_Binary,     &
             write_3D_Binary,     &
             read_3D_Binary

        ! Publicly available variables
        ! These arrays are only used when reading an output file of unknown size
      real(kind=ip), dimension(:,:)  ,allocatable,public :: B_XY
      real(kind=ip), dimension(:,:,:),allocatable,public :: B_XYZ

      contains
      !------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  deallocate_Binary
!
!  Called from: Ash3d_PostProc.F90
!  Arguments:
!    none
!
!  This subroutine deallocates Binary variables
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine deallocate_Binary

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine deallocate_Binary"
      endif;enddo

      if(allocated(B_XY)) deallocate(B_XY)

      end subroutine deallocate_Binary

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  write_2D_Binary
!
!  Called from: output_results and Ash3d_PostProc.F90
!  Arguments:
!    nx            = x length of output array OutVar
!    ny            = y length of output array OutVar
!    OutVar        = 2-d array to be written to binary file
!    VarMask       = 2-d logical that toggles data v.s. Fill_Value
!    Fill_Value    = number used for No-data
!    filename_root = root name of file (20 characters)
!
!  Subroutine that writes out 2-D arrays in binary format
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine write_2D_Binary(nx,ny,OutVar,VarMask,Fill_Value,filename_root)

      use io_data,       only : &
         isFinal_TS,iout3d,WriteTimes

      integer          ,intent(in) :: nx
      integer          ,intent(in) :: ny
      real(kind=ip)    ,intent(in) :: OutVar(nx,ny)
      logical          ,intent(in) :: VarMask(nx,ny)
      character(len=6) ,intent(in) :: Fill_Value
      character(len=20),intent(in) :: filename_root

      real(kind=op)  :: OVar(nx,ny)
      real(kind=op)  :: FValue
      integer :: i,j
      character(len= 9)  :: cio
      character(len=50)  :: filename_out
      integer            :: iostatus
      character(len=120) :: iomessage
      character(len=80)  :: linebuffer080

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine write_2D_Binary"
      endif;enddo

      read(Fill_Value,*,iostat=iostatus,iomsg=iomessage)FValue
      linebuffer080 = Fill_Value
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer080,iomessage)

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
      if(index(filename_root,'ArrivalTime').gt.0)then
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
        open(unit=fid_bin2dout,file=trim(adjustl(filename_out)), &
          status='replace', action='write', &
          access='direct',recl=4*nx*ny)
      else
       open(unit=fid_bin2dout,file=trim(adjustl(filename_out)), &
          status='replace', action='write', &
          access='direct',recl=8*nx*ny)
      endif
      write(fid_bin2dout,rec=1)((OVar(i,j),i=1,nx),j=1,ny)
      close(fid_bin2dout)

      end subroutine write_2D_Binary

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  read_2D_Binary
!
!  Called from: Ash3d_PostProc.F90
!  Arguments:
!    nx       = x length of output array OutVar
!    ny       = y length of output array OutVar
!    OutVar   = 2-d array to be written to binary file
!    filename = name of file (80 characters)
!
!  Subroutine that writes out 2-D arrays in binary format
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine read_2D_Binary(nx,ny,filename)

      integer          ,intent(in)  :: nx
      integer          ,intent(in)  :: ny
      character(len=80),intent(in)  :: filename

      real(kind=op)      :: OVar(nx,ny)
      integer            :: i,j
      integer            :: iostatus
      character(len=120) :: iomessage
      character(len=80)  :: linebuffer080

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine read_2D_Binary"
      endif;enddo

      if(op.eq.4)then
        open(unit=fid_bin2dout,file=trim(adjustl(filename)), &
          status='old', action='read', &
          access='direct',recl=4*nx*ny)
      else
       open(unit=fid_bin2dout,file=trim(adjustl(filename)), &
          status='old', action='read', &
          access='direct',recl=8*nx*ny)
      endif
      read(fid_bin2dout,rec=1,iostat=iostatus,iomsg=iomessage)((OVar(i,j),i=1,nx),j=1,ny)
      linebuffer080 = "Binary read"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
      close(fid_bin2dout)

      if(.not.allocated(B_XY)) allocate(B_XY(1:nx,1:ny))
      B_XY(1:nx,1:ny) = real(OVar(1:nx,1:ny),kind=ip)

      end subroutine read_2D_Binary

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  write_3D_Binary
!
!  Called from: output_results and Ash3d_PostProc.f90
!  Arguments:
!    cio        = time string to be inserted into filename; either '________final'
!                 or yyyymmddhh.h
!    nx         = x length of output array OutVar
!    ny         = y length of output array OutVar
!    nz         = z length of output array OutVar
!    ashcon_tot = 3d array containing the sum of all tephra concentration bins (1:n_gs_max)
!
!  This subroutine writes out 3-D arrays in binary format
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine write_3D_Binary(cio,nx,ny,nz,ashcon_tot)

      character(len=13) ,intent(in) :: cio
      integer           ,intent(in) :: nx
      integer           ,intent(in) :: ny
      integer           ,intent(in) :: nz
      real(kind=op)     ,intent(in) :: ashcon_tot(nx,ny,nz)

      integer :: i,j,k

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine write_3D_Binary"
      endif;enddo

      ! Write out data in raw binary form

      ! 3-D total tephra concentration
      if(op.eq.4)then
        open(unit=fid_bin3dout,file='3d_tephra_fall_'//cio//'.raw', &
          status='replace', action='write', &
          access='direct',recl=4*nx*ny*nz)
      else
        open(unit=fid_bin3dout,file='3d_tephra_fall_'//cio//'.raw', &
          status='replace', action='write', &
          access='direct',recl=8*nx*ny*nz)
      endif
      write(fid_bin3dout,rec=1)(((ashcon_tot(i,j,k),i=1,nx),j=1,ny),k=1,nz)
      close(fid_bin3dout)

      end subroutine write_3D_Binary

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  read_3D_Binary
!
!  Called from: Ash3d_PostProc.F90
!  Arguments:
!    nx       = x length of output array OutVar
!    ny       = y length of output array OutVar
!    OutVar   = 2-d array to be written to binary file
!    filename = name of file (80 characters)
!
!  Subroutine that writes out 2-D arrays in binary format
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine read_3D_Binary(nx,ny,nz,filename)

      integer          ,intent(in)  :: nx
      integer          ,intent(in)  :: ny
      integer          ,intent(in)  :: nz
      character(len=80),intent(in)  :: filename

      real(kind=op)      :: OVar3d(nx,ny,nz)
      integer            :: i,j,k
      integer            :: iostatus
      character(len=120) :: iomessage
      character(len=80)  :: linebuffer080

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine read_3D_Binary"
      endif;enddo

      if(op.eq.4)then
        open(unit=fid_bin3dout,file=trim(adjustl(filename)), &
          status='old', action='read', &
          access='direct',recl=4*nx*ny*nz)
      else
       open(unit=fid_bin3dout,file=trim(adjustl(filename)), &
          status='old', action='read', &
          access='direct',recl=8*nx*ny*nz)
      endif
      read(fid_bin3dout,rec=1,iostat=iostatus,iomsg=iomessage) &
           (((OVar3d(i,j,k),i=1,nx),j=1,ny),k=1,nz)
      linebuffer080 = "Binary read"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer080,iomessage)

      close(fid_bin3dout)

      if(.not.allocated(B_XYZ)) allocate(B_XYZ(1:nx,1:ny,1:nz))
      B_XYZ(1:nx,1:ny,1:nz) = real(OVar3d(1:nx,1:ny,1:nz),kind=ip)

      end subroutine read_3D_Binary

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module Ash3d_Binary_IO

!##############################################################################
