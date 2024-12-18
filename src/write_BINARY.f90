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
             write_3D_Binary,    &
             read_3D_Binary,     &
             BigEnd_2int,        &
             BigEnd_4int,        &
             LitEnd_2int,        &
             LitEnd_4int,        &
             BigEnd_4real,       &
             BigEnd_8real,       &
             LitEnd_4real,       &
             LitEnd_8real

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
      character(len= 50) :: linebuffer050 
      character(len= 80) :: linebuffer080

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine write_2D_Binary"
      endif;enddo

      read(Fill_Value,*,iostat=iostatus,iomsg=iomessage)FValue
      linebuffer080 = Fill_Value
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)

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
      character(len= 50) :: linebuffer050 
      character(len= 80) :: linebuffer080

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
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
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
      character(len= 50) :: linebuffer050 
      character(len= 80) :: linebuffer080

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
      linebuffer050 = "Reading file from binary file"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)

      close(fid_bin3dout)

      if(.not.allocated(B_XYZ)) allocate(B_XYZ(1:nx,1:ny,1:nz))
      B_XYZ(1:nx,1:ny,1:nz) = real(OVar3d(1:nx,1:ny,1:nz),kind=ip)

      end subroutine read_3D_Binary

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  BigEnd_2int
!
!  Called from: Currently not called
!  Arguments:
!    isLit = logical variable that is true if the input r is Little-endian
!    r     = 2-byte integer to be converted
!
!  This function returns the input variable 'r' in 2-byte integer in Big-endian
!  format.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      function BigEnd_2int(isLit,r)

      implicit none

      integer(kind=2) :: BigEnd_2int
      logical         :: isLit
      integer(kind=2) :: r

      integer(kind=2) :: s  = 0
      integer(kind=2) :: t  = 0
      integer         :: bl = 8  ! bit length of the move

      ! Older compilers accepted this equivalence approach, but some newer ones
      ! complained.  The transfer->mvbits->transfer approach seems more accepted,
      ! though it required fortran 95 or newer.

      !integer(kind=2) :: ii(2), jj(2)

      !equivalence (s,ii)
      !equivalence (t,jj)
      !if(isLit)then
      !  ! We need r to be big-endian, but the system is little-endian
      !  ! swap the bytes
      !  s = r
      !  jj(1) = ii(2)
      !  jj(2) = ii(1)
      !  BigEnd_2int = t
      !else
      !  BigEnd_2int = r
      !endif

      if(isLit)then
        ! We need r to be big-endian, but the system is little-endian
        ! swap the bytes
        !  map r onto the 2-byte dummy integer
        s = transfer(r,s)
        ! move 8-bit packets to the reverse positions
        call mvbits(s,0*bl,bl,t,1*bl)
        call mvbits(s,1*bl,bl,t,0*bl)
        ! map t onto the desired output
        BigEnd_2int = transfer (t,BigEnd_2int)
      else
        BigEnd_2int = r
      endif

      end function BigEnd_2int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  BigEnd_4int
!
!  Called from: write_ShapeFile_Polyline
!  Arguments:
!    isLit = logical variable that is true if the input r is Little-endian
!    r     = 4-byte integer to be converted
!
!  This function returns the input variable 'r' in 4-byte integer in Big-endian
!  format.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      function BigEnd_4int(isLit,r)

      implicit none

      integer(kind=4) :: BigEnd_4int
      logical         :: isLit
      integer(kind=4) :: r

      integer(kind=4) :: s  = 0
      integer(kind=4) :: t  = 0
      integer         :: bl = 8  ! bit length of the move

      ! Older compilers accepted this equivalence approach, but some newer ones
      ! complained.  The transfer->mvbits->transfer approach seems more accepted,
      ! though it required fortran 95 or newer.

!      integer(kind=1) :: ii(4), jj(4)
!      equivalence (s,ii)
!      equivalence (t,jj)

!      if(isLit)then
!        ! We need r to be big-endian, but the system is little-endian
!        ! swap the bytes
!        s = r
!        jj(1) = ii(4)
!        jj(2) = ii(3)
!        jj(3) = ii(2)
!        jj(4) = ii(1)
!        BigEnd_4int = t
!      else
!        BigEnd_4int = r
!      endif

      if(isLit)then
        ! We need r to be big-endian, but the system is little-endian
        ! swap the bytes
        !  map r onto the 2-byte dummy integer
        s = transfer(r,s)
        ! move 8-bit packets to the reverse positions
        call mvbits(s,0*bl,bl,t,3*bl)
        call mvbits(s,1*bl,bl,t,2*bl)
        call mvbits(s,2*bl,bl,t,1*bl)
        call mvbits(s,3*bl,bl,t,0*bl)
        ! map t onto the desired output
        BigEnd_4int = transfer (t,BigEnd_4int)
      else
        BigEnd_4int = r
      endif

      end function BigEnd_4int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  LitEnd_2int
!
!  Called from: write_ShapeFile_Polyline
!  Arguments:
!    isLit = logical variable that is true if the input r is Little-endian
!    r     = 2-byte integer to be converted
!
!  This function returns the input variable 'r' in 2-byte integer in Little-endian
!  format.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      function LitEnd_2int(isLit,r)

      implicit none

      integer(kind=2) :: LitEnd_2int
      logical         :: isLit
      integer(kind=2) :: r

      integer(kind=2) :: s  = 0
      integer(kind=2) :: t  = 0
      integer         :: bl = 8  ! bit length of the move

      ! Older compilers accepted this equivalence approach, but some newer ones
      ! complained.  The transfer->mvbits->transfer approach seems more accepted,
      ! though it required fortran 95 or newer.

!      integer(kind=1) :: ii(2), jj(2)

!      equivalence (s,ii)
!      equivalence (t,jj)
!      if(isLit)then
!        LitEnd_2int = r
!      else
!        ! We need r to be little-endian, but the system is big-endian
!        ! swap the bytes
!        s = r
!        jj(1) = ii(2)
!        jj(2) = ii(1)
!        LitEnd_2int = t
!      endif

      if(isLit)then
        LitEnd_2int = r
      else
        ! We need r to be little-endian, but the system is big-endian
        ! swap the bytes
        !  map r onto the 2-byte dummy integer
        s = transfer(r,s)
        ! move 8-bit packets to the reverse positions
        call mvbits(s,0*bl,bl,t,1*bl)
        call mvbits(s,1*bl,bl,t,0*bl)
        ! map t onto the desired output
        LitEnd_2int = transfer (t,LitEnd_2int)
      endif

      end function LitEnd_2int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  LitEnd_4int
!
!  Called from: write_ShapeFile_Polyline
!  Arguments:
!    isLit = logical variable that is true if the input r is Little-endian
!    r     = 4-byte integer to be converted
!
!  This function returns the input variable 'r' in 4-byte integer in Little-endian
!  format.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      function LitEnd_4int(isLit,r)

      implicit none

      integer(kind=4) :: LitEnd_4int
      logical         :: isLit
      integer(kind=4) :: r

      integer(kind=4) :: s  = 0
      integer(kind=4) :: t  = 0
      integer         :: bl = 8  ! bit length of the move

      ! Older compilers accepted this equivalence approach, but some newer ones
      ! complained.  The transfer->mvbits->transfer approach seems more accepted,
      ! though it required fortran 95 or newer.

!      integer(kind=1) :: ii(4), jj(4)

!      equivalence (s,ii)
!      equivalence (t,jj)
!      if(isLit)then
!        LitEnd_4int = r
!      else
!        ! We need r to be little-endian, but the system is big-endian
!        ! swap the bytes
!        s = r
!        jj(1) = ii(4)
!        jj(2) = ii(3)
!        jj(3) = ii(2)
!        jj(4) = ii(1)
!        LitEnd_4int = t
!      endif

      if(isLit)then
        LitEnd_4int = r
      else
        ! We need r to be little-endian, but the system is big-endian
        ! swap the bytes
        s = transfer(r,s)
        ! move 8-bit packets to the reverse positions
        call mvbits(s,0*bl,bl,t,3*bl)
        call mvbits(s,1*bl,bl,t,2*bl)
        call mvbits(s,2*bl,bl,t,1*bl)
        call mvbits(s,3*bl,bl,t,0*bl)
        ! map t onto the desired output
        LitEnd_4int = transfer (t,LitEnd_4int)
      endif

      end function LitEnd_4int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  BigEnd_4real
!
!  Called from: Currently not called
!  Arguments:
!    isLit = logical variable that is true if the input r is Little-endian
!    r     = 4-byte real to be converted
!
!  This function returns the input variable 'r' in 4-byte real in Big-endian
!  format.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      function BigEnd_4real(isLit,r)

      implicit none

      real(kind=4)    :: BigEnd_4real
      logical         :: isLit
      real(kind=4)    :: r

      integer(kind=4) :: s  = 0
      integer(kind=4) :: t  = 0
      integer         :: bl = 8  ! bit length of the move

      ! Older compilers accepted this equivalence approach, but some newer ones
      ! complained.  The transfer->mvbits->transfer approach seems more accepted,
      ! though it required fortran 95 or newer.

!      integer(kind=1) :: ii(4), jj(4)
!      real(kind=4)    :: s, t
!
!      equivalence (s,ii)
!      equivalence (t,jj)
!
!      if(isLit)then
!        ! We need r to be big-endian, but the system is little-endian
!        ! swap the bytes
!        s = r
!        jj(1) = ii(4)
!        jj(2) = ii(3)
!        jj(3) = ii(2)
!        jj(4) = ii(1)
!        BigEnd_4real = t
!      else
!        BigEnd_4real = r
!      endif


      if(isLit)then
        ! We need r to be Big-endian, but the system is Little-endian
        ! swap the bytes
        s = transfer(r,s)
        ! move 8-bit packets to the reverse positions
        call mvbits(s,0*bl,bl,t,3*bl)
        call mvbits(s,1*bl,bl,t,2*bl)
        call mvbits(s,2*bl,bl,t,1*bl)
        call mvbits(s,3*bl,bl,t,0*bl)
        ! map t onto the desired output
        BigEnd_4real = transfer (t,BigEnd_4real)
      else
        BigEnd_4real = r
      endif

      end function BigEnd_4real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  BigEnd_8real
!
!  Called from: Currently not called
!  Arguments:
!    isLit = logical variable that is true if the input r is Little-endian
!    r     = 8-byte real to be converted
!
!  This function returns the input variable 'r' in 8-byte real in Big-endian
!  format.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      function BigEnd_8real(isLit,r)

      implicit none

      real(kind=8)    :: BigEnd_8real
      logical         :: isLit
      real(kind=8)    :: r

      integer(kind=8) :: s  = 0
      integer(kind=8) :: t  = 0
      integer         :: bl = 8  ! bit length of the move

      ! Older compilers accepted this equivalence approach, but some newer ones
      ! complained.  The transfer->mvbits->transfer approach seems more accepted,
      ! though it required fortran 95 or newer.

!      integer(kind=1) :: ii(8), jj(8)
!      real(kind=8)    :: s, t
!
!      equivalence (s,ii)
!      equivalence (t,jj)
!
!      if(isLit)then
!        ! We need r to be big-endian, but the system is little-endian
!        ! swap the bytes
!        s = r
!        jj(1) = ii(8)
!        jj(2) = ii(7)
!        jj(3) = ii(6)
!        jj(4) = ii(5)
!        jj(5) = ii(4)
!        jj(6) = ii(3)
!        jj(7) = ii(2)
!        jj(8) = ii(1)
!        BigEnd_8real = t
!      else
!        BigEnd_8real = r
!      endif


      if(isLit)then
        ! We need r to be Big-endian, but the system is Little-endian
        ! swap the bytes
        s = transfer(r,s)
        ! move 8-bit packets to the reverse positions
        call mvbits(s,0*bl,bl,t,7*bl)
        call mvbits(s,1*bl,bl,t,6*bl)
        call mvbits(s,2*bl,bl,t,5*bl)
        call mvbits(s,3*bl,bl,t,4*bl)
        call mvbits(s,4*bl,bl,t,3*bl)
        call mvbits(s,5*bl,bl,t,2*bl)
        call mvbits(s,6*bl,bl,t,1*bl)
        call mvbits(s,7*bl,bl,t,0*bl)
        ! map t onto the desired output
        BigEnd_8real = transfer (t,BigEnd_8real)
      else
        BigEnd_8real = r
      endif

      end function BigEnd_8real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  LitEnd_4real
!
!  Called from: Currently not called
!  Arguments:
!    isLit = logical variable that is true if the input r is Little-endian
!    r     = 4-byte real to be converted
!
!  This function returns the input variable 'r' in 4-byte real in Little-endian
!  format.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      function LitEnd_4real(isLit,r)

      implicit none

      real(kind=4)    :: LitEnd_4real
      logical         :: isLit
      real(kind=4)    :: r

      integer(kind=4) :: s  = 0
      integer(kind=4) :: t  = 0
      integer         :: bl = 8  ! bit length of the move

      ! Older compilers accepted this equivalence approach, but some newer ones
      ! complained.  The transfer->mvbits->transfer approach seems more accepted,
      ! though it required fortran 95 or newer.

!      integer(kind=1) :: ii(4), jj(4)
!
!      equivalence (s,ii)
!      equivalence (t,jj)
!      if(isLit)then
!        LitEnd_4real = r
!      else
!        ! We need r to be little-endian, but the system is big-endian
!        ! swap the bytes
!        s = r
!        jj(1) = ii(4)
!        jj(2) = ii(3)
!        jj(3) = ii(2)
!        jj(4) = ii(1)
!        LitEnd_4real = t
!      endif

      if(isLit)then
        LitEnd_4real = r
      else
        ! We need r to be little-endian, but the system is big-endian
        ! swap the bytes
        s = transfer(r,s)
        ! move 8-bit packets to the reverse positions
        call mvbits(s,0*bl,bl,t,3*bl)
        call mvbits(s,1*bl,bl,t,2*bl)
        call mvbits(s,2*bl,bl,t,1*bl)
        call mvbits(s,3*bl,bl,t,0*bl)
        ! map t onto the desired output
        LitEnd_4real = transfer (t,LitEnd_4real)
      endif

      end function LitEnd_4real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  LitEnd_8real
!
!  Called from: write_ShapeFile_Polyline
!  Arguments:
!    isLit = logical variable that is true if the input r is Little-endian
!    r     = 8-byte real to be converted
!
!  This function returns the input variable 'r' in 8-byte real in Little-endian
!  format.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      function LitEnd_8real(isLit,r)

      implicit none

      real(kind=8)    :: LitEnd_8real
      logical         :: isLit
      real(kind=8)    :: r

      integer(kind=8) :: s  = 0
      integer(kind=8) :: t  = 0
      integer         :: bl = 8  ! bit length of the move

      ! Older compilers accepted this equivalence approach, but some newer ones
      ! complained.  The transfer->mvbits->transfer approach seems more accepted,
      ! though it required fortran 95 or newer.

!      integer(kind=1) :: ii(8), jj(8)
!
!      equivalence (s,ii)
!      equivalence (t,jj)
!      if(isLit)then
!        LitEnd_8real = r
!      else
!        ! We need r to be little-endian, but the system is big-endian
!        ! swap the bytes
!        s = r
!        jj(1) = ii(8)
!        jj(2) = ii(7)
!        jj(3) = ii(6)
!        jj(4) = ii(5)
!        jj(5) = ii(4)
!        jj(6) = ii(3)
!        jj(7) = ii(2)
!        jj(8) = ii(1)
!        LitEnd_8real = t
!      endif

      if(isLit)then
        LitEnd_8real = r
      else
        ! We need r to be little-endian, but the system is big-endian
        ! swap the bytes
        s = transfer(r,s)
        ! move 8-bit packets to the reverse positions
        call mvbits(s,0*bl,bl,t,7*bl)
        call mvbits(s,1*bl,bl,t,6*bl)
        call mvbits(s,2*bl,bl,t,5*bl)
        call mvbits(s,3*bl,bl,t,4*bl)
        call mvbits(s,4*bl,bl,t,3*bl)
        call mvbits(s,5*bl,bl,t,2*bl)
        call mvbits(s,6*bl,bl,t,1*bl)
        call mvbits(s,7*bl,bl,t,0*bl)
        ! map t onto the desired output
        LitEnd_8real = transfer (t,LitEnd_8real)
      endif

      end function LitEnd_8real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module Ash3d_Binary_IO

!##############################################################################
