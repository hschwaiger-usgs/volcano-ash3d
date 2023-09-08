!##############################################################################
      module citywriter

      use precis_param,  only : &
        ip,op,sp,dp

      use io_units

      implicit none

      !private

      contains

!##############################################################################

      subroutine citylist(outCode,inlonLL,inlonUR,inlatLL,inlatUR,maxcities, &
                          CityLon_out,CityLat_out,CityName_out)

      use global_param,  only : &
        DirDelim

      use io_data,       only : &
         Ash3dHome

!     subroutine that reads from a list of cities and figures out which ones to include
!     on a map.

      implicit none

      integer      ,intent(in) :: outCode ! 0 for list only, 1 for GMT, 2 for gnuplot
      real(kind=ip),intent(in) :: inlonLL
      real(kind=ip),intent(in) :: inlonUR
      real(kind=ip),intent(in) :: inlatLL
      real(kind=ip),intent(in) :: inlatUR
      integer      ,intent(in) :: maxcities

      real(kind=ip)    ,dimension(maxcities),intent(out) :: CityLon_out
      real(kind=ip)    ,dimension(maxcities),intent(out) :: CityLat_out
      character(len=26),dimension(maxcities),intent(out) :: CityName_out

      integer            :: Iostatus = 1
      integer            :: i, ncities, nread
      integer            :: resolution                              ! # of cells in x and y
      character(len=26)  :: CityName
      character(len=133) :: inputline
      real(kind=ip)      :: lonLL,lonUR,latLL,latUR
      real(kind=ip)      :: CityLat, CityLon
      real(kind=ip)      :: dlat, dlon, cell_width, cell_height
      real(kind=ip)      :: minspace_x, minspace_y
      logical            :: IsOkay                     ! true if city is not near any others
      logical            :: IsThere,IsThere2

      character(len=130)             :: CityMasterFile

      CityName_out = ''           ! set default values
      CityLon_out  = 0.0_ip
      CityLat_out  = 0.0_ip
      resolution   = 100          ! number of cells in width & height

      ! All city longitudes are between -180 and 180 degrees.
      ! Make sure the requested computational domain is in the same range.
      if (inlonLL.gt.180.0_ip)then
        lonLL = inlonLL-360.0_ip
      else
        lonLL = inlonLL
      endif
      if (inlonUR.gt.180.0_ip)then
        lonUR = inlonUR-360.0_ip
      else
        lonUR = inlonUR
      endif
      ! if the model domain wraps across the prime meridian add 360 to longitude
      if (inlonLL.gt.inlonUR)lonUR = inlonUR + 360.0_ip
      latLL = inlatLL
      latUR = inlatUR

      dlat = latUR - latLL
      dlon = lonUR - lonLL
      cell_width = dlon/resolution
      cell_height = dlat/resolution
      minspace_x  = 3.0_ip*cell_width
      minspace_y  = 3.0_ip*cell_height

      nread   = 0
      ncities = 0

      CityMasterFile = 'world_cities.txt'
      inquire( file=trim(adjustl(CityMasterFile)), exist=IsThere )
      if(.not.IsThere)then
        CityMasterFile = trim(Ash3dHome) // &
                          DirDelim // 'share' // &
                          DirDelim // 'post_proc' // &
                          DirDelim // 'world_cities.txt'
        inquire( file=trim(adjustl(CityMasterFile)), exist=IsThere2 )
        if(.not.IsThere2)then
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(outlog(io),*)"ERROR: Could not find file: ",trim(adjustl(CityMasterFile))
            write(outlog(io),*)"       Skipping cities"
          endif;enddo
          ncities = 0
          return
        endif
      endif
      open(unit=fid_cities,file=trim(adjustl(CityMasterFile)),status='old',action='read')

      ! skip the first line
      read(fid_cities,*)

      do while ((ncities.lt.maxcities).and.(Iostatus.ge.0))
        read(fid_cities,'(a133)',iostat=Iostatus) inputline
        read(inputline,2) CityLon, CityLat, CityName
2       format(f16.4,f15.4,a26)
        if ((CityLon.gt.lonLL).and.(CityLon.lt.lonUR).and. &
            (CityLat.gt.latLL).and.(CityLat.lt.latUR)) then
          ! Make sure this city is not near any others
          IsOkay=.true.
          call space_checker(maxcities,CityLon_out,CityLat_out,ncities, & 
                             CityLon,CityLat, &
                             minspace_x,minspace_y,IsOkay)
          if (IsOkay) then
            ncities = ncities+1
            CityName_out(ncities) = CityName
            CityLon_out(ncities)  = CityLon
            CityLat_out(ncities)  = CityLat
          endif
          ! if the model domain crosses over the prime meridian
        elseif ((CityLon+360.0_ip.gt.lonLL).and.(CityLon+360.0_ip.lt.lonUR).and. &
                (CityLat.gt.latLL).and.(CityLat.lt.latUR)) then
          ! Make sure this city is not near any others
          IsOkay=.true.
          call space_checker(maxcities,CityLon_out,CityLat_out,ncities, & 
                             CityLon,CityLat, &
                             minspace_x,minspace_y,IsOkay)
          if (IsOkay) then
            ncities = ncities+1
            CityName_out(ncities) = CityName
            CityLon_out(ncities)  = CityLon
            CityLat_out(ncities)  = CityLat
          endif
        endif
        nread=nread+1
      enddo
      close(fid_cities)
      if(outCode.gt.0)then
        if(ncities.gt.0) then
          open(unit=fid_citiesxy,file='cities.xy',status='replace',action='write')
          if(outCode.eq.1)then
            ! for gmt
            do i=1,ncities
              write(fid_citiesxy,4) CityLon_out(i),CityLat_out(i),CityName_out(i)
4             format(2f10.4,'  10  0  9  BL    ',a26)
            enddo
          elseif(outCode.eq.2)then
            do i=1,ncities
              ! for gnuplot
              write(fid_citiesxy,*) real(CityLon_out(i),kind=4),real(CityLat_out(i),kind=4),&
                          '"',trim(adjustl(CityName_out(i))),'"'
!3             format(2f10.4,1x,a26)
            enddo
          else
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"outCode not recognized. No output file written"
            endif;enddo
          endif
          close(fid_citiesxy)
        endif
      endif

      end subroutine citylist
         
!***************************************************************************************

      subroutine space_checker(maxcities,CityLon_out,CityLat_out,ncities, & 
                                  CityLon,CityLat, &
                                  minspace_x,minspace_y,IsOkay)

      implicit none

      integer      ,intent(in)    :: maxcities
      real(kind=ip),intent(in)    :: CityLon_out(maxcities)
      real(kind=ip),intent(in)    :: CityLat_out(maxcities)
      integer      ,intent(in)    :: ncities
      real(kind=ip),intent(in)    :: CityLon
      real(kind=ip),intent(in)    :: CityLat
      real(kind=ip),intent(in)    :: minspace_x
      real(kind=ip),intent(in)    :: minspace_y
      logical      ,intent(inout) :: IsOkay                 ! true if city is not near any others

      integer            :: icity

      do icity=1,ncities
        if ((abs(CityLon_out(icity)-CityLon).lt.minspace_x).and. &
            (abs(CityLat_out(icity)-CityLat).lt.minspace_y)) then
          IsOkay=.false.
          exit
        endif
      enddo

      return

      end subroutine space_checker

      end module
!##############################################################################

