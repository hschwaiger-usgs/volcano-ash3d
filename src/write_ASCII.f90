!##############################################################################
!
!      Ash3d_ASCII_IO module
!
!  This module manages all output to ASCII files
!
!      subroutine vprofileopener
!      subroutine vprofilewriter(itime)
!      subroutine vprofilecloser
!      subroutine write_2D_ASCII(nx,ny,OutVar,VarMask,Fill_Value,filename_root)
!      subroutine write_2D_ASCII_flt(nx,ny,IsLL,x1,y1,IsCC,dx,dy,Fill_Value,OutVar,filename)
!      subroutine write_2D_ASCII_flt_regular(nx,ny,IsLL,x1,y1,IsCC,dx,dy,Fill_Value,OutVar,filename)
!      subroutine write_2D_ASCII_int(nx,ny,IsLL,x1,y1,IsCC,dx,dy,Fill_Value,OutVar,filename)
!      subroutine write_2D_ASCII_csv(nx,ny,IsLL,x1,y1,IsCC,dx,dy,Fill_Value,OutVar,filename)
!      subroutine read_2D_ASCII(filename)
!      subroutine write_3D_ASCII(cio)
!      subroutine read_3D_ASCII(filename)
!      subroutine Write_PointData_Airports_ASCII
!
!##############################################################################

      module Ash3d_ASCII_IO

      use precis_param

      use io_units

      implicit none
      !implicit none (type, external)

        ! Set everything to private by default
      private

        ! Publicly available subroutines/functions
      public deallocate_ASCII,  &
             vprofileopener,    &
             vprofilewriter,    &
             vprofilecloser,    &
             write_2D_ASCII,    &
             write_2D_ASCII_flt,&
             write_2D_ASCII_flt_regular,&
             !write_2D_ASCII_int,&
             write_2D_ASCII_csv,&
             read_2D_ASCII,     &
             write_3D_ASCII,    &
             read_3D_ASCII,     &
             Write_PointData_Airports_ASCII

        ! Publicly available variables
        ! These arrays are only used when reading an output file of unknown size
      integer,       dimension(:,:)  ,allocatable,public :: A_XY_int
      real(kind=ip), dimension(:,:)  ,allocatable,public :: A_XY
      real(kind=ip), dimension(:,:,:),allocatable,public :: A_XYZ
      integer      ,public :: A_nx
      integer      ,public :: A_ny
      integer      ,public :: A_nz
      real(kind=ip),public :: A_xll
      real(kind=ip),public :: A_yll
      real(kind=ip),public :: A_zll
      real(kind=ip),public :: A_dx
      real(kind=ip),public :: A_dy
      real(kind=ip),public :: A_dz
      real(kind=ip),public :: A_Fill
      integer      ,public :: A_Fill_int
      logical      ,public :: A_IsInt    = .false.

      contains
      !------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  deallocate_ASCII
!
!  Called from: Ash3d_PostProc.F90
!  Arguments:
!    none
!
!  This subroutine deallocates ASCII variables
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine deallocate_ASCII

      do io=1,2;if(VB(io) <= verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine deallocate_ASCII"
      endif;enddo

      if(allocated(A_XY))     deallocate(A_XY)
      if(allocated(A_XY_int)) deallocate(A_XY_int)
      if(allocated(A_XYZ))    deallocate(A_XYZ)

      do io=1,2;if(VB(io) <= verbosity_debug1)then
        write(outlog(io),*)"     Exited Subroutine deallocate_ASCII"
      endif;enddo

      return

      end subroutine deallocate_ASCII

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  vprofileopener
!
!  Called from: output_results and Ash3d_PostProc.F90
!  Arguments:
!    none
!
!  This subroutine opens files for vertical profiles
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine vprofileopener

      use io_data,       only : &
         MAXPROFILES,nvprofiles,Site_vprofile,x_vprofile, y_vprofile

      use mesh,          only : &
         nzmax,z_cc_pd

      integer :: i,j

      integer ::  ionumber            ! number of output file
      character(len=14)  :: cio

      do io=1,2;if(VB(io) <= verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine vprofileopener"
      endif;enddo

      do i=1,nvprofiles
        ionumber = fid_vprofbase + i-1
        if (i < 10) then
          write(cio,1) i
1         format('vprofile0',i1,'.txt')
        elseif (i < MAXPROFILES) then
          write(cio,2) i
2         format('vprofile',i2,'.txt')
        else
          do io=1,2;if(VB(io) <= verbosity_error)then
            write(errlog(io),*)"ERROR: Too many vertical profiles."
            write(errlog(io),*)"nvprofiles must be < ",MAXPROFILES
              write(errlog(io),*)"       Please increase MAXNPROFILES and recompile."
              write(errlog(io),*)" Ash3d_VariableModules.f90:io_data:MAXPROFILES"
          endif;enddo
          stop 1
        endif
        open(unit=ionumber,file=cio,status='replace',action='write')
        write(ionumber,*)'Vertical profile data for location: ',trim(adjustl(Site_vprofile(i)))
        write(ionumber,3)x_vprofile(i), y_vprofile(i)
3       format( &
       'x:',f10.3,/, &
       'y:',f10.3,/, &
       '                          Output is ash concentration in mg/m3',/, &
       '                          elevation (km) ---->',/, &
       'date-time         hrs')
        write(ionumber,4) (z_cc_pd(j), j=1,nzmax)
4       format(50f16.3)
      enddo

      do io=1,2;if(VB(io) <= verbosity_debug1)then
        write(outlog(io),*)"     Exited Subroutine vprofileopener"
      endif;enddo

      return

      end subroutine vprofileopener

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  vprofilewriter
!
!  Called from: Ash3d.F90 and Ash3d_PostProc.F90
!  Arguments:
!    itime = time step index
!
!  This subroutine writes data on vertical profiles
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine vprofilewriter(itime)

      use global_param,  only : &
         KG_2_MG,KM3_2_M3

      use io_data,       only : &
         nvprofiles

      use mesh,          only : &
         nzmax

      use time_data,     only : &
         SimStartHour,time,BaseYear,useLeap,OutputOffset

      use Output_Vars,    only : &
         CLOUDCON_THRESH,pr_ash

      integer, intent(in) :: itime

      integer :: i,k
      integer   :: ionumber
      character(len=13)  :: cio

      INTERFACE
        character (len=13) function HS_yyyymmddhh_since(HoursSince,byear,useLeaps)
          implicit none
          !implicit none (type, external)
          integer        ,parameter   :: dp        = 8  ! double precision
          real(kind=dp)  ,intent(in)  :: HoursSince
          integer        ,intent(in)  :: byear
          logical        ,intent(in)  :: useLeaps
        end function HS_yyyymmddhh_since
      END INTERFACE

      do io=1,2;if(VB(io) <= verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine vprofilewriter"
      endif;enddo

      do i=1,nvprofiles
        ! don't write if there's no ash
        if(maxval(pr_ash(1:nzmax,itime,i)) < CLOUDCON_THRESH*KG_2_MG/KM3_2_M3) cycle
        ionumber = fid_vprofbase + i-1
        cio = HS_yyyymmddhh_since(SimStartHour+time+OutputOffset,&
                                  BaseYear,useLeap)
        write(ionumber,1) cio, time, (pr_ash(k,itime,i), k=1,nzmax)  ! write tot. concen in mg/m3

1       format(a13,',',f10.3,',',50(e15.3,','))
      enddo

      do io=1,2;if(VB(io) <= verbosity_debug1)then
        write(outlog(io),*)"     Exited Subroutine vprofilewriter"
      endif;enddo

      return

      end subroutine vprofilewriter

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  vprofilecloser
!
!  Called from: output_results and Ash3d_PostProc.F90
!  Arguments:
!    none
!
!  This subroutine closes vertical profile files
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine vprofilecloser

      use io_data,       only : &
         nvprofiles

      integer  ::  i, ionumber

      do io=1,2;if(VB(io) <= verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine vprofilecloser"
      endif;enddo

      do i=1,nvprofiles
        ionumber = fid_vprofbase + i-1
        close(ionumber)
      enddo

      do io=1,2;if(VB(io) <= verbosity_debug1)then
        write(outlog(io),*)"     Exited Subroutine vprofilecloser"
      endif;enddo

      return

      end subroutine vprofilecloser

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  write_2D_ASCII
!
!  Called from: output_results and Ash3d_PostProc.F90
!  Arguments:
!    nx            = x length of output array OutVar
!    ny            = y length of output array OutVar
!    OutVar        = 2-d array to be written to ASCII file
!    VarMask       = 2-d logical that toggles data v.s. Fill_Value
!    Fill_Value    = number used for No-data
!    filename_root = root name of file (20 characters)
!
!  Subroutine that writes out 2-D arrays in ESRI ASCII raster format in the style of Ash3d.
!  This style includes: corner-specification, file-naming, 10-float columns with rows separated
!  by a newline.
!  Full format specification is given at the following web sites:
!   https://help.arcgis.com/en/arcgisdesktop/10.0/help/index.html#/ESRI_ASCII_raster_format/009t0000000z000000/
!   https://en.wikipedia.org/wiki/Esri_grid
!  This format can be post-processed with gmt converting to grid files with
!   gmt grdconvert out.dat=ef out.grd
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine write_2D_ASCII(nx,ny,OutVar,VarMask,Fill_Value,filename_root)

      use global_param,  only  : &
         KM_2_M

      use mesh,          only : &
         dx,dy,de,dn,IsLatLon,latLL,lonLL,xLL,yLL

      use io_data,       only : &
         isFinal_TS,iout3d,WriteTimes

      integer          ,intent(in) :: nx,ny
      real(kind=ip)    ,intent(in) :: OutVar(nx,ny)
      logical          ,intent(in) :: VarMask(nx,ny)
      character(len=6) ,intent(in) :: Fill_Value
      character(len=20),intent(in) :: filename_root

      real(kind=op)  :: OVar(nx,ny)
      real(kind=op)  :: FValue
      integer :: i,j
      character (len=9)  :: cio
      character(len=50)  :: filename_out
      integer            :: iostatus
      character(len=120) :: iomessage
      character(len= 50) :: linebuffer050
      character(len= 80) :: linebuffer080

      do io=1,2;if(VB(io) <= verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine write_2D_ASCII"
      endif;enddo

      read(Fill_Value,*,iostat=iostatus,iomsg=iomessage)FValue
      linebuffer080 = Fill_Value
      linebuffer050 = "Reading FValue from ASCII file"
      if(iostatus /= 0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)

      if(isFinal_TS)then
        cio='____final'
      else
        if (WriteTimes(iout3d) < 10.0_ip) then
          write(cio,1) WriteTimes(iout3d)
1         format('00',f4.2,'hrs')
        elseif (WriteTimes(iout3d) < 100.0_ip) then
          write(cio,2) WriteTimes(iout3d)
2         format('0',f5.2,'hrs')
        else
          write(cio,3) WriteTimes(iout3d)
3         format(f6.2,'hrs')
        endif
      endif
      if(index(filename_root,'outvar') > 0)then
        ! If this subroutine is called with a filename root of 'outvar', then just
        ! write to file 'outvar.dat'
        write(filename_out,*)trim(adjustl(filename_root)),'.dat'
      elseif(index(filename_root,'ArrivalTime') > 0)then
          ! For the special cases of DepositArrivalTime.dat and
          ! CloudArrivalTime.dat
        write(filename_out,*)trim(adjustl(filename_root)),'.dat'
      else
        write(filename_out,*)trim(adjustl(filename_root)),cio,'.dat'
      endif

      open(unit=fid_ascii2dout,file=trim(adjustl(filename_out)), status='replace',action='write',err=2500)

      write(fid_ascii2dout,3000) nx        ! write header values
      write(fid_ascii2dout,3001) ny
      if (IsLatLon) then
        if (lonLL > 180.0_ip)then
          ! GIS software seem to prefer the domain -180->180
          write(fid_ascii2dout,3002) lonLL -360.0_ip
        else
          write(fid_ascii2dout,3002) lonLL
        endif
        write(fid_ascii2dout,3003) latLL
        write(fid_ascii2dout,3004) de,dn
      else
        write(fid_ascii2dout,3002) xLL*KM_2_M    ! convert xLL from km to meters so ArcMap can read it
        write(fid_ascii2dout,3003) yLL*KM_2_M    ! same with yLL
        write(fid_ascii2dout,3004) dx*KM_2_M,dy*KM_2_M    ! and with dx and dy
      endif
      write(fid_ascii2dout,3005)Fill_Value

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

      ! Write out arrays
      do j=ny,1,-1
        write(fid_ascii2dout,3006) (OVar(i,j), i=1,nx)
        write(fid_ascii2dout,*)" "               ! make a blank line between rows
      enddo

      close(fid_ascii2dout)

!     format statements
3000  format('NCOLS ',i5)
3001  format('NROWS ',i5)
3002  format('XLLCORNER ',f15.3)
3003  format('YLLCORNER ',f15.3)
3004  format('CELLSIZE ',2f15.3)
3005  format('NODATA_VALUE ',a6)
!3006  format(10f15.3)               ! Older ASCII output file from Ash3d used this format
3006  format(10f18.6)

      do io=1,2;if(VB(io) <= verbosity_debug1)then
        write(outlog(io),*)"     Exited Subroutine write_2D_ASCII"
      endif;enddo

      return

!     Error traps
2500  do io=1,2;if(VB(io) <= verbosity_error)then
        write(errlog(io),*) 'Error opening output file ASCII_output_file.txt.  Program stopped'
      endif;enddo
      stop 1

      end subroutine write_2D_ASCII

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  write_2D_ASCII_flt
!
!  Called from: topography and some auxilary tools
!  Arguments:
!    nx            = x length of output array OutVar
!    ny            = y length of output array OutVar
!    IsLL          = logical flag for Lon/Lat vs projected
!    x1            = x (or lon) coordinate of Lower-left cell
!    y1            = y (or lat) coordinate of Lower-left cell
!    IsCC          = logical flag for cell-centered vs corner grid registration
!    dx            = cell length in x
!    dy            = cell length in y
!    Fill_Value    = number used for No-data (given as 6-char)
!    OutVar        = 2-d array to be written to ASCII file
!    filename      = name of output file (20 characters)
!
!  Subroutine that writes out 2-D arrays in ESRI ASCII raster format.
!  Full format specification is given at the following web sites:
!   https://help.arcgis.com/en/arcgisdesktop/10.0/help/index.html#/ESRI_ASCII_raster_format/009t0000000z000000/
!   https://en.wikipedia.org/wiki/Esri_grid
!  This format can be post-processed with gmt converting to grid files with
!   gmt grdconvert out.dat=ef out.grd
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine write_2D_ASCII_flt(nx,ny,IsLL,x1,y1,IsCC,dx,dy,Fill_Value,OutVar,filename)

      use global_param,  only  : &
         KM_2_M

      integer          ,intent(in) :: nx,ny
      logical          ,intent(in) :: IsLL
      real(kind=sp)    ,intent(in) :: x1,y1
      logical          ,intent(in) :: IsCC
      real(kind=sp)    ,intent(in) :: dx,dy
      character(len=6) ,intent(in) :: Fill_Value
      real(kind=sp)    ,intent(in) :: OutVar(nx,ny)
      character(len=20),intent(in) :: filename

      real(kind=op)  :: FValue
      integer :: i,j
      character(len=50)  :: filename_out
      integer            :: iostatus
      character(len=120) :: iomessage
      character(len= 50) :: linebuffer050
      character(len= 80) :: linebuffer080

      do io=1,2;if(VB(io) <= verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine write_2D_ASCII_flt"
      endif;enddo

      read(Fill_Value,*,iostat=iostatus,iomsg=iomessage)FValue
      linebuffer080 = Fill_Value
      linebuffer050 = "Reading FValue from ASCII file"
      if(iostatus /= 0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)

      write(filename_out,*)trim(adjustl(filename))

      open(unit=fid_ascii2dout,file=trim(adjustl(filename_out)), status='replace',action='write',err=2500)

      write(fid_ascii2dout,3000) nx        ! write header values
      write(fid_ascii2dout,3001) ny
      if (IsLL) then
        if (x1 > 180.0_ip)then
          ! GIS software seem to prefer the domain -180->180
          if(IsCC)then
            write(fid_ascii2dout,3012) x1 -360.0_sp
          else
            write(fid_ascii2dout,3002) x1 -360.0_sp
          endif
        else
          if(IsCC)then
            write(fid_ascii2dout,3012) x1
          else
            write(fid_ascii2dout,3002) x1
          endif
        endif
        if(IsCC)then
          write(fid_ascii2dout,3013) y1
        else
          write(fid_ascii2dout,3003) y1
        endif
        write(fid_ascii2dout,3004) dx,dy
      else
        if(IsCC)then
          write(fid_ascii2dout,3012) x1*KM_2_M    ! convert xLL from km to meters so ArcMap can read it
          write(fid_ascii2dout,3013) y1*KM_2_M    ! same with yLL
        else
          write(fid_ascii2dout,3002) x1*KM_2_M    ! convert xLL from km to meters so ArcMap can read it
          write(fid_ascii2dout,3003) y1*KM_2_M    ! same with yLL
        endif
        write(fid_ascii2dout,3004) dx*KM_2_M,dy*KM_2_M    ! and with dx and dy
      endif
      write(fid_ascii2dout,3005)Fill_Value

      ! Write out arrays
      do j=ny,1,-1
        write(fid_ascii2dout,*) (OutVar(i,j), i=1,nx)
      enddo

      close(fid_ascii2dout)

!     format statements
3000  format('NCOLS ',i5)
3001  format('NROWS ',i5)
3002  format('XLLCORNER ',f15.3)
3003  format('YLLCORNER ',f15.3)
3012  format('XLLCENTER ',f15.3)
3013  format('YLLCENTER ',f15.3)
3004  format('CELLSIZE ',2f15.3)
3005  format('NODATA_VALUE ',a6)

      do io=1,2;if(VB(io) <= verbosity_debug1)then
        write(outlog(io),*)"     Exited Subroutine write_2D_ASCII_flt"
      endif;enddo

      return

!     Error traps
2500  do io=1,2;if(VB(io) <= verbosity_error)then
        write(errlog(io),*) 'Error opening output file ASCII_output_file.txt.  Program stopped'
      endif;enddo
      stop 1

      end subroutine write_2D_ASCII_flt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  write_2D_ASCII_flt_regular
!
!  Called from: write_2Dmap_PNG_python
!  Arguments:
!    nx            = x length of output array OutVar
!    ny            = y length of output array OutVar
!    IsLL          = logical flag for Lon/Lat vs projected
!    x1            = x (or lon) coordinate of Lower-left cell
!    y1            = y (or lat) coordinate of Lower-left cell
!    IsCC          = logical flag for cell-centered vs corner grid registration
!    dx            = cell length in x
!    dy            = cell length in y
!    Fill_Value    = number used for No-data (given as 6-char)
!    OutVar        = 2-d array to be written to ASCII file
!    filename      = name of output file (20 characters)
!
!  Subroutine that writes out 2-D arrays in ESRI ASCII raster format, resampling onto a regular grid.
!  Full format specification is given at the following web sites:
!   https://help.arcgis.com/en/arcgisdesktop/10.0/help/index.html#/ESRI_ASCII_raster_format/009t0000000z000000/
!   https://en.wikipedia.org/wiki/Esri_grid
!  This format can be post-processed with gmt converting to grid files with
!   gmt grdconvert out.dat=ef out.grd
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine write_2D_ASCII_flt_regular(nx,ny,IsLL,x1,y1,IsCC,dx,dy,Fill_Value,OutVar,filename)

      use global_param,  only  : &
         KM_2_M

      integer          ,intent(in) :: nx,ny
      logical          ,intent(in) :: IsLL
      real(kind=sp)    ,intent(in) :: x1,y1
      logical          ,intent(in) :: IsCC
      real(kind=sp)    ,intent(in) :: dx,dy
      character(len=6) ,intent(in) :: Fill_Value
      real(kind=sp)    ,intent(in) :: OutVar(nx,ny)
      character(len=20),intent(in) :: filename

      real(kind=op)  :: FValue
      integer :: i,j
      character(len=50)  :: filename_out
      integer            :: iostatus
      character(len=120) :: iomessage
      character(len= 50) :: linebuffer050
      character(len= 80) :: linebuffer080
      integer             :: nx_loc, ny_loc
      real(kind=sp)       :: dx_loc,x_width,y_width
      real(kind=sp),dimension(nx) :: x_in  ! cell-center coordinate of input grid
      real(kind=sp),dimension(ny) :: y_in
      real(kind=sp),dimension(:),allocatable    :: x_loc  ! cell-center coordinate of output grid
      real(kind=sp),dimension(:),allocatable    :: y_loc
      real(kind=sp),dimension(:,:),allocatable  :: OutVar_loc
      real(kind=sp)       :: A11,A12,A21,A22
      real(kind=sp)       :: z11,z12,z21,z22
      real(kind=sp)       :: xl,xr,yl,yr
      real(kind=sp)       :: AREA
      integer             :: ix,iy
      real(kind=sp)       :: xnow,ynow,znow

      do io=1,2;if(VB(io) <= verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine write_2D_ASCII_flt_regular"
      endif;enddo

      read(Fill_Value,*,iostat=iostatus,iomsg=iomessage)FValue
      linebuffer080 = Fill_Value
      linebuffer050 = "Reading FValue from ASCII file"
      if(iostatus /= 0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)

      write(filename_out,*)trim(adjustl(filename))

      open(unit=fid_ascii2dout,file=trim(adjustl(filename_out)), status='replace',action='write',err=2500)

      ! Find the smaller of dx,dy
      dx_loc = min(dx,dy)
      do i=1,nx
        x_in(i) = x1+0.5_sp*dx + (i-1)*dx
      enddo
      do j=1,ny
        y_in(j) = y1+0.5_sp*dy + (j-1)*dy
      enddo
      x_width = nx*dx
      y_width = ny*dy
      nx_loc = floor(x_width/dx_loc)
      ny_loc = floor(y_width/dx_loc)
      allocate(x_loc(nx_loc))
      allocate(y_loc(ny_loc))
      allocate(OutVar_loc(nx_loc,ny_loc))
      do i=1,nx_loc
        x_loc(i) = x1+0.5_sp*dx_loc + (i-1)*dx_loc
      enddo
      do j=1,ny_loc
        y_loc(j) = y1+0.5_sp*dx_loc + (j-1)*dx_loc
      enddo

      ! Loop over all points of the new grid and interpolate
      do i=1,nx_loc
        xnow = min(max(x_loc(i),x_in(1)),x_in(nx))  ! make sure the point is within domain of x_in
        ix = floor((xnow-x1)/dx)
        ix = min(max(ix,1),nx-1)                    ! make sure ix is between 1 and nx-1
        xl = x_in(ix)
        xr = x_in(ix+1)
        do j=1,ny_loc
          ynow = min(max(y_loc(j),y_in(1)),y_in(ny))  ! make sure the point is within domain of y_in
          iy = floor((ynow-y1)/dy)
          iy = min(max(iy,1),ny-1)                    ! make sure iy is between 1 and ny-1
          yl = y_in(iy)
          yr = y_in(iy+1)

        ! The bilinear interpolation is the sum of the values in the four corners
        ! that surround the point in question, weighted by the areas of the rectangles
        ! to the upper left, upper right, lower left, and lower right of the point in
        ! question, divided by the area of the entire square in which the point sits.
        A11  = (xnow-  xl) * (ynow -   yl)    !area of lower  left rectangle
        A21  = (xr  -xnow) * (ynow -   yl)    !area of lower right rectangle
        A12  = (xnow-  xl) * (yr   - ynow)    !area of upper  left rectangle
        A22  = (xr  -xnow) * (yr   - ynow)    !area of upper right rectangle
        AREA = (yr  -  yl) * (xr   -   xl)    !total area of node

        ! values at corners of the rectangle
        z11 = OutVar(ix  ,iy  )
        z21 = OutVar(ix+1,iy  )
        z12 = OutVar(ix  ,iy+1)
        z22 = OutVar(ix+1,iy+1)

        ! Thickness at the airport
        znow = (1.0_sp/AREA) * &  !reciprocal of area of large rectangle
           (z22 * A11    + &  !z22 times area in lower left
            z12 * A21    + &  !z12 times area in lower right
            z21 * A12    + &  !z21 times area in upper left
            z11 * A22)        !z11 times area in upper right
        OutVar_loc(i,j) = znow
        enddo
      enddo

      write(fid_ascii2dout,3000) nx_loc        ! write header values
      write(fid_ascii2dout,3001) ny_loc
      if (IsLL) then
        if (x1 > 180.0_ip)then
          ! GIS software seem to prefer the domain -180->180
          if(IsCC)then
            write(fid_ascii2dout,3012) x1 -360.0_sp
          else
            write(fid_ascii2dout,3002) x1 -360.0_sp
          endif
        else
          if(IsCC)then
            write(fid_ascii2dout,3012) x1
          else
            write(fid_ascii2dout,3002) x1
          endif
        endif
        if(IsCC)then
          write(fid_ascii2dout,3013) y1
        else
          write(fid_ascii2dout,3003) y1
        endif
        write(fid_ascii2dout,3004) dx_loc,dx_loc
      else
        if(IsCC)then
          write(fid_ascii2dout,3012) x1*KM_2_M    ! convert xLL from km to meters so ArcMap can read it
          write(fid_ascii2dout,3013) y1*KM_2_M    ! same with yLL
        else
          write(fid_ascii2dout,3002) x1*KM_2_M    ! convert xLL from km to meters so ArcMap can read it
          write(fid_ascii2dout,3003) y1*KM_2_M    ! same with yLL
        endif
        write(fid_ascii2dout,3004) dx_loc*KM_2_M,dx_loc*KM_2_M    ! and with dx and dy
      endif
      write(fid_ascii2dout,3005)Fill_Value

      ! Write out arrays
      do j=ny_loc,1,-1
        write(fid_ascii2dout,*) (OutVar_loc(i,j), i=1,nx_loc)
      enddo

      deallocate(x_loc)
      deallocate(y_loc)
      deallocate(OutVar_loc)

      close(fid_ascii2dout)

!     format statements
3000  format('NCOLS ',i5)
3001  format('NROWS ',i5)
3002  format('XLLCORNER ',f15.3)
3003  format('YLLCORNER ',f15.3)
3012  format('XLLCENTER ',f15.3)
3013  format('YLLCENTER ',f15.3)
3004  format('CELLSIZE ',2f15.3)
3005  format('NODATA_VALUE ',a6)

      do io=1,2;if(VB(io) <= verbosity_debug1)then
        write(outlog(io),*)"     Exited Subroutine write_2D_ASCII_flt_regular"
      endif;enddo

      return

!     Error traps
2500  do io=1,2;if(VB(io) <= verbosity_error)then
        write(errlog(io),*) 'Error opening output file ASCII_output_file.txt.  Program stopped'
      endif;enddo
      stop 1

      end subroutine write_2D_ASCII_flt_regular

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  write_2D_ASCII_csv
!
!  Called from:
!  Arguments:
!    nx            = x length of output array OutVar
!    ny            = y length of output array OutVar
!    x1            = x (or lon) coordinate of Lower-left cell
!    y1            = y (or lat) coordinate of Lower-left cell
!    dx            = cell length in x
!    dy            = cell length in y
!    OutVar        = 2-d array to be written to ASCII file
!    filename      = name of output file (20 characters)
!
!  Subroutine that writes out 2-D data as a comma-separated-values file: x,y,OutVar(:,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine write_2D_ASCII_csv(nx,ny,x1,y1,dx,dy,OutVar,filename)

      use global_param,  only  : &
         KM_2_M

      integer          ,intent(in) :: nx,ny
      real(kind=sp)    ,intent(in) :: x1,y1
      real(kind=sp)    ,intent(in) :: dx,dy
      real(kind=sp)    ,intent(in) :: OutVar(nx,ny)
      character(len=20),intent(in) :: filename

      integer :: i,j
      character(len=50)  :: filename_out
      integer            :: iostatus
      character(len=120) :: iomessage
      character(len= 50) :: linebuffer050
      character(len= 80) :: linebuffer080
      real(kind=sp),dimension(nx) :: x_in  ! cell-center coordinate of input grid
      real(kind=sp),dimension(ny) :: y_in

      do io=1,2;if(VB(io) <= verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine write_2D_ASCII_csv"
      endif;enddo

      write(filename_out,*)trim(adjustl(filename)),'.csv'

      open(unit=fid_ascii2dout,file=trim(adjustl(filename_out)), status='replace',action='write',err=2500)

      do i=1,nx
        x_in(i) = x1+0.5_sp*dx + (i-1)*dx
      enddo
      do j=1,ny
        y_in(j) = y1+0.5_sp*dy + (j-1)*dy
      enddo

      do j=1,ny
        do i=1,nx
          write(fid_ascii2dout,'(f10.4,a3,f10.4,a3,g15.5)')x_in(i),' , ',y_in(j),' , ',OutVar(i,j)
        enddo
      enddo

      close(fid_ascii2dout)

      do io=1,2;if(VB(io) <= verbosity_debug1)then
        write(outlog(io),*)"     Exited Subroutine write_2D_ASCII_csv"
      endif;enddo

      return

!     Error traps
2500  do io=1,2;if(VB(io) <= verbosity_error)then
        write(errlog(io),*) 'Error opening output file ASCII_output_file.txt.  Program stopped'
      endif;enddo
      stop 1

      end subroutine write_2D_ASCII_csv


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  read_2D_ASCII
!
!  Called from: Ash3d_ASCII_check.f90 and Ash3d_PostProc.F90
!  Arguments:
!    filename = full filename (80 characters)
!
!  Subroutine that reads in 2-D arrays in ESRI ASCII raster format and
!  populates A_nx,A_ny,A_XY,A_xll,A_yll,A_dx,A_dy,A_Fill
!  Note that this is not a generic ERSI ASCII reader. This subroutine is
!  designed to read ASCII files written by the subroutine write_2D_ASCII above.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine read_2D_ASCII(filename)

      character(len=80),intent(in) :: filename

      integer           :: i,j
      character(len=080):: linebuffer080
      integer           :: iostatus
      character(len=120):: iomessage
      character(len=20) :: tst_str
      integer           :: substr_pos1
      logical           :: IsAsh3dASCII

      do io=1,2;if(VB(io) <= verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine read_2D_ASCII"
      endif;enddo

      ! Initialization
      IsAsh3dASCII = .true.

      open(unit=fid_ascii2din,file=trim(adjustl(filename)), status='old',action='read',err=2500)

      ! The header of the ESRI ASCII file has 6 lines
      ! NCOLS (ncols)                           int
      ! NROWS (nrows)                           int
      ! XLLCORNER/XLLCENTER (xllcorner)         double
      ! YLLCORNER/YLLCENTER (yllcorner)         double
      ! CELLSIZE (cellsize)                     double (and maybe a second double)
      ! NODATA_VALUE (NODATA_value)             double or int

      read(fid_ascii2din,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      read(linebuffer080(7:),*,iostat=iostatus,iomsg=iomessage)A_nx
      if(iostatus /= 0)then
        ! We might have an empty file
        ! Issue warning and return
        do io=1,2;if(VB(io) <= verbosity_error)then
          write(errlog(io),*) 'Error reading file ',trim(adjustl(filename))
          write(errlog(io),*) 'Check for zero-length file.'
        endif;enddo
        return
      endif
      read(fid_ascii2din,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      read(linebuffer080(7:),*,iostat=iostatus,iomsg=iomessage)A_ny
      read(fid_ascii2din,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      read(linebuffer080(10:),*,iostat=iostatus,iomsg=iomessage)A_xll
      read(fid_ascii2din,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      read(linebuffer080(10:),*,iostat=iostatus,iomsg=iomessage)A_yll
      read(fid_ascii2din,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      read(linebuffer080(10:),*,iostat=iostatus,iomsg=iomessage)A_dx
      ! test if two cell lengths are provided
      read(linebuffer080(10:),*,iostat=iostatus,iomsg=iomessage)A_dx,A_dy
      if(iostatus /= 0)then
        ! only one cell size provided; copy dx to dy
        A_dy = A_dx
      endif

      ! Try to read the nondata value, first as a float
      read(fid_ascii2din,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      tst_str(1:20) = linebuffer080(13:32)
      substr_pos1 = index(tst_str,'.')
      if(substr_pos1 > 0)then
        ! period found, assume data are floats
        read(linebuffer080(13:),*,iostat=iostatus,iomsg=iomessage)A_Fill
        A_IsInt = .false.
      else
        read(linebuffer080(13:),*,iostat=iostatus,iomsg=iomessage)A_Fill_int
        A_IsInt = .true.
      endif

      if(A_IsInt)then
        allocate(A_XY_int(A_nx,A_ny))
        do j=A_ny,1,-1
          ! This format reads the full line of integers
          read(fid_ascii2din,*,err=2600,iostat=iostatus,iomsg=iomessage) (A_XY_int(i,j), i=1,A_nx)
        enddo
      else
        allocate(A_XY(A_nx,A_ny))
        do j=A_ny,1,-1
          ! This format ID directs to read 10 floats at a time, matching the Ash3d ASCII output format
          read(fid_ascii2din,3006,iostat=iostatus,iomsg=iomessage) (A_XY(i,j), i=1,A_nx)
          if(iostatus /= 0)then
            ! If there is an error, then the ESRI ASCII file might not be an Ash3d-generated file.
            ! Try reading the whole line of A_nx values without the blank-line separator
            IsAsh3dASCII = .false.
            exit
          endif
          read(fid_ascii2din,*,iostat=iostatus,iomsg=iomessage) (A_XY(i,j), i=1,A_nx)
          if(iostatus /= 0)then
            IsAsh3dASCII = .false.
            exit
          endif
        enddo
      endif

      if(.not.IsAsh3dASCII)then
        ! Run this part only if the above loop fails (non-Ash3d-type of ASCII file)
        ! First rewind and zip past the 6 header lines
        rewind(fid_ascii2din)
        read(fid_ascii2din,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        read(fid_ascii2din,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        read(fid_ascii2din,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        read(fid_ascii2din,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        read(fid_ascii2din,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        read(fid_ascii2din,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        ! read the raster file from the top down as per the ESRI ASCII specification
        do j=A_ny,1,-1
          read(fid_ascii2din,*,err=2600,iostat=iostatus,iomsg=iomessage) (A_XY(i,j), i=1,A_nx)
        enddo
      endif

      close(fid_ascii2din)

!     format statements
!3006  format(10f15.3)               ! Older ASCII output file from Ash3d used this format
3006  format(10f18.6)

      do io=1,2;if(VB(io) <= verbosity_debug1)then
        write(outlog(io),*)"     Exited Subroutine read_2D_ASCII"
      endif;enddo

      return

!     Error traps
2500  do io=1,2;if(VB(io) <= verbosity_error)then
        write(errlog(io),*) 'Error opening ASCII file. Program stopped'
      endif;enddo
      stop 1

2600  do io=1,2;if(VB(io) <= verbosity_error)then
        write(errlog(io),*) 'Error reading from ASCII file.'
      endif;enddo
      stop 1

      end subroutine read_2D_ASCII

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  write_3D_ASCII
!
!  Called from: output_results
!  Arguments:
!    cio = time string to be inserted into filename; either '________final'
!          or yyyymmddhh.h
!
!  Subroutine that writes out 3-D arrays in ESRI ASCII raster format
!  Format specification is given at the following web sites:
!   https://help.arcgis.com/en/arcgisdesktop/10.0/help/index.html#/ESRI_ASCII_raster_format/009t0000000z000000/
!   https://en.wikipedia.org/wiki/Esri_grid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine write_3D_ASCII(cio,nx,ny,nz,ashcon_tot)

      use mesh,          only : &
         lon_cc_pd,lat_cc_pd,IsLatLon,&
         x_cc_pd,y_cc_pd,z_cc_pd

      character(len=13) ,intent(in) :: cio
      integer           ,intent(in) :: nx
      integer           ,intent(in) :: ny
      integer           ,intent(in) :: nz
      real(kind=op)     ,intent(in) :: ashcon_tot(nx,ny,nz)

      integer :: i,j,k
      character(len=32) :: DepOutfileName

      do io=1,2;if(VB(io) <= verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine write_3D_ASCII"
      endif;enddo

      ! Output data in ASCII format
      DepOutfileName='3d_tephra_fall_'//cio//'.dat'
      open(unit=fid_ascii3dout,file=DepOutfileName,status='replace',action='write')
      write(fid_ascii3dout,*)'VARIABLES = "X","Y","Z","AshConc"'
      write(fid_ascii3dout,3000) 'ZONE I = ',nx,' J = ',ny,' K = ',nz

      do k=1,nz
        do j=1,ny
          do i=1,nx
            if (IsLatLon) then
              write(fid_ascii3dout,'(3(4x,f20.3),g20.8)') &
                lon_cc_pd(i), lat_cc_pd(j), z_cc_pd(k), ashcon_tot(i,j,k)
            else
              write(fid_ascii3dout,'(3(4x,f20.3),g20.8)') &
                x_cc_pd(i), y_cc_pd(j), z_cc_pd(k), ashcon_tot(i,j,k)
            endif
          enddo
        enddo
      enddo

      close(fid_ascii3dout)


      do io=1,2;if(VB(io) <= verbosity_debug1)then
        write(outlog(io),*)"     Exited Subroutine write_3D_ASCII"
      endif;enddo

      return

!     format statements
3000  format(a9,i5,a5,i5,a5,i5)

      end subroutine write_3D_ASCII

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  read_3D_ASCII
!
!  Called from: Ash3d_PostProc.F90
!  Arguments:
!    filename = root name of file (20 characters)
!
!  Subroutine that reads in 3-D arrays in ESRI ASCII raster format and
!  populates A_nx,A_ny,A_nz,A_XYZ,A_xll,A_yll,A_zll,A_dx,A_dy,A_dz,A_Fill
!  Note that this is not a generic ERSI ASCII reader. This subroutine is
!  designed to read ASCII files written by the subroutine write_3D_ASCII above.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine read_3D_ASCII(filename)

      character(len=80),intent(in) :: filename

      integer :: i,j,k
      integer :: iostatus
      character(len=120) :: iomessage
      character(len= 50) :: linebuffer050
      character(len=130) :: linebuffer130
      real(kind=ip) :: value1,value2,value3

      do io=1,2;if(VB(io) <= verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine read_3D_ASCII"
      endif;enddo

      open(unit=fid_ascii3din,file=trim(adjustl(filename)), status='old',action='read',err=2500)

      ! Read the header lines
      read(fid_ascii3din,'(a130)',iostat=iostatus,iomsg=iomessage)linebuffer130
      linebuffer050 = "Reading line from ASCII file"
      if(iostatus /= 0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer130(1:80),iomessage)
      read(fid_ascii3din,3000,err=2600,iostat=iostatus,iomsg=iomessage) A_nx,A_ny,A_nz
      if(iostatus /= 0)then
        ! We might have an empty file
        ! Issue warning and return
        do io=1,2;if(VB(io) <= verbosity_error)then
          write(errlog(io),*) 'Error reading file ',trim(adjustl(filename))
          write(errlog(io),*) 'Check for zero-length file.'
        endif;enddo
        return
      endif

      if(.not.allocated(A_XYZ)) allocate(A_XYZ(1:A_nx,1:A_ny,1:A_nz))
      do k=1,A_nz
        do j=1,A_ny
          do i=1,A_nx
            read(fid_ascii3din,*,iostat=iostatus,iomsg=iomessage)linebuffer130
            read(linebuffer130,'(3(4x,f20.3),g20.8)',iostat=iostatus,iomsg=iomessage) &
                   value1,value2,value3,A_XYZ(i,j,k)
            linebuffer050 = "Reading line from ASCII file"
            if(iostatus /= 0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer130(1:80),iomessage)

            if(i == 1.and.j == 1.and.k == 1)then
              A_xll = value1
              A_yll = value2
              A_zll = value3
            endif
            if(i == 2) A_dx = value1 - A_xll
          enddo
          if(j == 2) A_dy = value2 - A_yll
        enddo
        if(k == 2) A_dz = value3 - A_zll
      enddo

      close(fid_ascii3din)

!     format statements
3000  format(9x,i5,5x,i5,5x,i5)

      do io=1,2;if(VB(io) <= verbosity_debug1)then
        write(outlog(io),*)"     Exited Subroutine read_3D_ASCII"
      endif;enddo

      return

!     Error traps
2500  do io=1,2;if(VB(io) <= verbosity_error)then
        write(errlog(io),*) 'Error opening ASCII file. Program stopped'
      endif;enddo
      stop 1

2600  do io=1,2;if(VB(io) <= verbosity_error)then
        write(errlog(io),*) 'Error reading from ASCII file.'
      endif;enddo
      stop 1

      end subroutine read_3D_ASCII

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Write_PointData_Airports_ASCII
!
!  Called from: output_results and Ash3d_PostProc.f90
!  Arguments:
!    none
!
!  Subroutine that writes the arrival time at airports to an ASCII file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Write_PointData_Airports_ASCII

      use global_param,  only : &
         M_2_MM,UseCalcFallVel

      use io_data,       only : &
         infile,WriteGSD,WriteAirportFile_ASCII,VolcanoName

      use mesh,          only : &
         nsmax

       use solution,      only : &
         DepositGranularity

      use time_data,     only : &
         time,SimStartHour,OutputOffset,BaseYear,useLeap,RunStartMinute,&
         RunStartYear,RunStartMonth,RunStartDay,RunStartHr

      use Airports,      only : &
         Airport_AshArrivalTime,Airport_CloudArrivalTime, &
         Airport_thickness,Airport_AshDuration, &
         Airport_CloudDuration,nairports,Airport_CloudArrived,&
         Airport_AshArrived,Airport_depRate,Airport_i,Airport_j,&
         Airport_Longitude,Airport_Latitude,Airport_Name,&
         bilinear_thickness

      use Output_Vars,   only : &
         DepositThickness,DEPRATE_THRESH, &
         CloudLoad,CLOUDLOAD_THRESH

      use Tephra,        only : &
         n_gs_max,Tephra_gsdiam,Tephra_rho_m,Tephra_v_s

      use Source,        only : &
         neruptions,e_Duration,e_Volume,e_PlumeHeight

      integer             :: i
      integer             :: nWrittenOut
      character (len=20)  :: yyyymmddhh_ash, yyyymmddhh_cloud
      character (len=1)   :: cloud_morethan, deposit_morethan      !equals ">" if cloud is still overhead or ash is still falling
      character (len=13)  :: nwsthickness                ! nwp ashfall  terms (trace, minor, substantial, heavy, severe)
      integer             :: isize
      real(kind=ip)       :: longitude_now

      INTERFACE
        character (len=20) function HS_xmltime(HoursSince,byear,useLeaps)
          implicit none
          !implicit none (type, external)
          integer        ,parameter   :: dp        = 8  ! double precision
          real(kind=dp)  ,intent(in)  :: HoursSince
          integer        ,intent(in)  :: byear
          logical        ,intent(in)  :: useLeaps
        end function HS_xmltime
      END INTERFACE

      do io=1,2;if(VB(io) <= verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine Write_PointData_Airports_ASCII"
      endif;enddo

      ! Write values out to ASCII file
      if (WriteAirportFile_ASCII) then

        ! Write out source parameters for airports file
        open(unit=fid_asharrive,file='ash_arrivaltimes_airports.txt',status='replace',action='write',err=2000)
        write(fid_asharrive,98)  infile, RunStartYear,RunStartMonth,RunStartDay,RunStartHr, &
                              RunStartMinute, VolcanoName  !write infile, simulation time
        do i=1,neruptions  ! write source parameters
          write(fid_asharrive,99) i, HS_xmltime(SimStartHour+OutputOffset,&
                                           BaseYear,useLeap), &
                           e_Duration(i), e_PlumeHeight(i), e_PlumeHeight(i)*3280.8_ip, e_Volume(i)
        enddo
        write(fid_asharrive,995)
        if (WriteGSD) then
          ! if we're writing out grain sizes.
          if (UseCalcFallVel) then
            ! if fall velocity is calculated
            write(fid_asharrive,100) (Tephra_gsdiam(isize)*M_2_MM, isize=1,n_gs_max)
            write(fid_asharrive,101)
            write(fid_asharrive,102) (Tephra_rho_m(isize), isize=1,n_gs_max)
            write(fid_asharrive,103)
          else
            ! if fall velocity is specified
            write(fid_asharrive,110) (Tephra_v_s(isize), isize=1,n_gs_max)
          endif
        else
          ! if not writing out grain sizes
          write(fid_asharrive,1)
        endif
        nWrittenOut = 0
        do i=1,nairports                      ! write out the airports that are hit.
          if (Airport_AshArrived(i).or.Airport_CloudArrived(i)) then
            ! rank ash thickness in NWS rank
            Airport_thickness(i) = bilinear_thickness(i,DepositThickness)  ! interpolate to find thickness
            if (Airport_thickness(i) <= 0.7935_ip) then         ! <1/32" thickness
              nwsthickness="trace or less"
            elseif (Airport_thickness(i) <= 6.35_ip) then       ! <=1/4"
              nwsthickness="minor"
            elseif (Airport_thickness(i) <= 25.4_ip) then       ! <=1"
              nwsthickness="substantial"
            elseif (Airport_thickness(i) <= 101.6_ip) then      ! <=4"
              nwsthickness="heavy"
            else
              nwsthickness="severe"
            endif

            ! get yyyymmddhh of arrival
            yyyymmddhh_cloud = HS_xmltime(Airport_CloudArrivalTime(i)+SimStartHour+OutputOffset,&
                                          BaseYear,useLeap)
            if (Airport_AshArrived(i)) then
              yyyymmddhh_ash = HS_xmltime(Airport_AshArrivaltime(i)+SimStartHour+OutputOffset,&
                                          BaseYear,useLeap)
            else
              yyyymmddhh_ash = '0000-00-00T00:00:00Z'
            endif

            ! See whether cloud is still overhead, or whether ash is still
            ! falling
            if((Airport_AshArrived(i)).and.(Airport_depRate(i) > DEPRATE_THRESH)) then
              Airport_AshDuration(i) = time-Airport_AshArrivalTime(i)
              deposit_morethan = '>'
            else
              deposit_morethan = ' '
            endif
            if (CloudLoad(Airport_i(i),Airport_j(i)) > CLOUDLOAD_THRESH)then
              Airport_CloudDuration(i) = time-Airport_CloudArrivalTime(i)
              cloud_morethan = '>'
            else
              cloud_morethan = ' '
            endif

            if (Airport_Longitude(i) > 180.0_ip) then
              longitude_now = Airport_Longitude(i)-360.0_ip
            else
              longitude_now = Airport_Longitude(i)
            endif
            if (WriteGSD) then
              write(fid_asharrive,20) Airport_Name(i), Airport_Latitude(i), longitude_now, &
                       yyyymmddhh_cloud, Airport_CloudArrivalTime(i), cloud_morethan, Airport_CloudDuration(i), &
                       yyyymmddhh_ash, Airport_AshArrivalTime(i), deposit_morethan, Airport_AshDuration(i), &
                       Airport_thickness(i), nwsthickness, &
                       ((DepositGranularity(Airport_i(i),Airport_j(i),isize)/ &
                         sum(DepositGranularity(Airport_i(i),Airport_j(i),:))),isize=1,nsmax)  ! mass fraction of size i
            else
              write(fid_asharrive,2) Airport_Name(i), Airport_Latitude(i), longitude_now, &
                       yyyymmddhh_cloud, Airport_CloudArrivalTime(i), cloud_morethan, Airport_CloudDuration(i), &
                       yyyymmddhh_ash, Airport_AshArrivalTime(i), deposit_morethan, Airport_AshDuration(i), &
                       Airport_thickness(i), nwsthickness
            endif
            nWrittenOut = nWrittenOut + 1
          endif
        enddo
        if (nWrittenOut == 0) write(fid_asharrive,3)    ! if no airports are hit, say it in the file.
        write(fid_asharrive,120)                        ! write footnotes & caveats
        close(fid_asharrive)
      endif

      do io=1,2;if(VB(io) <= verbosity_debug1)then
        write(outlog(io),*)"     Exited Subroutine Write_PointData_Airports_ASCII"
      endif;enddo

      return

!     Error traps
2000  do io=1,2;if(VB(io) <= verbosity_error)then
        write(errlog(io),*)  'Error opening ash_arrivaltimes_airports.txt.  Program stopped.'
      endif;enddo
      stop 1

!     Format statements
1     format('---------------------------------------------------------------------------------------------------', &
             '----------------------------------------------------------------',/, &
             '                      LOCATION                          |                  CLOUD                  |', &
             '                               DEPOSIT                         |',/, &
             '                                                        |                                         |', &
             '                                                               |',/, &
             '                                                        |                                         |', &
             '                                                               |',/, &
             '(Airport code & ) Place name         Latitude Longitude |        Cloud Arrival Time      Duration |', &
             '      Deposit Arrival Time        Duration      Thickness      |',/, &
             '                                                        |                        hrs after        |', &
             '                          hrs after                            |',/, &
             '                                                        |      date/time UTC    start      hrs    |', &
             '       date/time UTC      start     hrs      mm   NWS rank     |')
98    format('ARRIVAL TIME OF ASH IN AREA MODELED BY ASH3D',/, &  !if WriteGSD=.true. and CalcFallVel=.true.
             'simulation using input file: ',A130,/,                        &
             'Model run date: ',i4,'.',i2.2,'.',i2,', time UTC: ',i4,':',i2.2,//, &
             '---------------------------------------------------------------------------------',/, &
             '             SOURCE PARAMETERS FOR SIMULATION OF: ',a30,'|',/, &
             'Pulse       start time                Duration       Plume height      volume   |',/, &
             '               UTC                      hrs          km     feet         km3    |')
99    format(i3,4x,a20,f16.2,f13.1,f9.0,f13.5,'  |')
995   format('---------------------------------------------------------------------------------',/)
100   format('---------------------------------------------------------------------------------------------------', &
             '---------------------------------------------------------------',/, &
             '                      LOCATION                          |                  CLOUD                  |', &
             '                               DEPOSIT                         |  GRAIN SIZE, MM',/, &
             '                                                        |                                         |', &
             '                                                               |',30f12.4)                                    !grain size, mm
101   format('                                                        |                                         |', &
             '                                                               |  DENSITY, KG/M3')
102   format('(Airport code & ) Place name         Latitude Longitude |        Cloud Arrival Time      Duration |', &
             '      Deposit Arrival Time        Duration      Thickness      |',30f12.1)                                    !density (kg/m3)
103   format('                                                        |                      hrs after          |', &
             '                                                        |',/, &
             '                                                        |      date/time UTC    start      hrs    |', &
             '       date/time UTC       hrs      hrs      mm   NWS rank     |  MASS FRACTION AT THE GIVEN GRAIN SIZE')

110   format('---------------------------------------------------------------------------------------------------', &
             '----------------------------------------------------------------',/, &
             '                      LOCATION                          |                  CLOUD                  |', &
             '                               DEPOSIT                         |', &
             '    MASS PER UNIT AREA AT THE GIVEN FALL SPEED (M/S)',/, &
             '                                                        |                                         |', &
             '                                                               |',/, &
             '                                                        |                                         |', &
             '                                                               |',/, &
             '(Airport code & ) Place name         Latitude Longitude |        Cloud Arrival Time      Duration |', &
             '      Deposit Arrival Time        Duration      Thickness      |',/, &
             '                                                        |                                         |', &
             '                                                               |',/, &
             '                                                        |      date/time UTC     hrs       hrs    |', &
             '       date/time UTC       hrs      hrs      mm   NWS rank     |',30f12.4)
120   format(//, &
             '---------------------------------------------------------------------------------------------------', &
             '----------------------------------------------------------------',//, &
             'NOTES ON ITEMS IN THIS TABLE:',/, &
             'LOCATION: If the location is an airport, the first three letters are the ICAO airport code',/, &
             'CLOUD DATA: Cloud arrival time is given in hours after the eruption start and the date and time', &
                          ' in UTC,in the format yyyy-mm-ddThh:mm:ssZ).  Duration is the',/, &
             '  number of hours during which the cloud is overhead (or, if there is a break in the cloud, the ', &
                          'time from the clouds first arrival until it last passes).',/, &
             '  A character ">" in front of the duration indicates that the cloud was still overhead at the end ', &
                          'of the simulation. The vertically integrated cloud mass must',/, &
             '  exceed a threshold value of 0.2 tonnes per square kilometer to be considered in these calculations.',/, &
             'DEPOSIT DATA:  The deposit arrival time is given in hours since eruption start or in the date and time ', &
                          'UTC as formatted in the cloud output (The cloud and ',/, &
             '   deposit may arrive at different times). "Deposit arrival time" is the time of arrival of the deposit ', &
                          'at a thickness exceeding 0.01 mm (0.0004 inches).',/, &
             '   An arrival time of -9999.00 hrs indicates the deposit did not arrive at this location.  Deposit ', &
                          'duration is the time period (hrs) over which the deposit',/, &
             '   was falling at a rate exceeding 0.01 mm/hr.  A ">" character before this number indicates that the ', &
                          'deposit was still falling at the end of the simulation.',/, &
             '   The thickness of the deposit is given in millimeters (left column) and as ranked according to the ', &
                          'following system devised by the U.S. National Weather Service:',/, &
             '     NWS Rank         Thickness',/, &
             '                        up to',/, &
             '                     mm      in.',/, &
             '     trace           0.8     1/32"',/, &
             '     minor           6.3     1/4"',/, &
             '     substantial    25.4     1"',/, &
             '     heavy          100      4"',/, &
             '     severe        >100     >4"',//, &
             'NOTE: This table is the estimate at time of issuance: changing conditions at the volcano may require ', &
                         'updating the forecast.')
2     format(a35,2f10.4,' |',2x,a20,f7.2,3x,a1,f5.2,'   |',2x,a20,f9.2,3x,a1,f5.2,f8.2,2x,a13,'|')
20    format(a35,2f10.4,' |',2x,a20,f7.2,3x,a1,f5.2,'   |',2x,a20,f9.2,3x,a1,f5.2,f8.2,2x,a13,'|',30e12.4)
3     format(/,'No airports affected by ash')

      end subroutine Write_PointData_Airports_ASCII

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module Ash3d_ASCII_IO

!##############################################################################
