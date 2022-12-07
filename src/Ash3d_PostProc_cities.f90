      subroutine citylist(outCode,inlonLL,inlonUR,inlatLL,inlatUR,maxcities, &
                          CityLon_out,CityLat_out,CityName_out)

      use precis_param

      use global_param,  only : &
        DirDelim

      use io_data,       only : &
         Ash3dHome

!     subroutine that reads from a list of cities and figures out which ones to include
!     on a map.

      implicit none

      integer     ,intent(in) :: outCode ! 0 for list only, 1 for GMT, 2 for gnuplot
      real(kind=8),intent(in) :: inlonLL
      real(kind=8),intent(in) :: inlonUR
      real(kind=8),intent(in) :: inlatLL
      real(kind=8),intent(in) :: inlatUR
      integer     ,intent(in) :: maxcities

      real(kind=8),dimension(maxcities),intent(out) :: CityLon_out
      real(kind=8),dimension(maxcities),intent(out) :: CityLat_out
      character(len=26),dimension(maxcities),intent(out) :: CityName_out

      integer            :: iostatus = 1
      integer            :: i, nargs, ncities, nmax, nread
      integer            :: status
      integer            :: resolution                              !# of cells in x and y
      character(len=26)  :: CityName
      character(len=133) :: inputline
      character(len=9)   :: lonLL_char, lonUR_char, latLL_char, latUR_char
      real(kind=8)       :: lonLL,lonUR,latLL,latUR
      real(kind=8)       :: CityLat, CityLon
      real(kind=8)       :: dlat, dlon, cell_width, cell_height
      real(kind=8)       :: minspace_x, minspace_y
      logical            :: IsOkay                 !true if city is not near any others
      logical            :: IsThere,IsThere2

      character(len=130)             :: CityMasterFile

      CityName_out = ''           !set default values
      CityLon_out  = 0.0_8
      CityLat_out  = 0.0_8
      resolution   = 100          !number of cells in width & height

      ! make sure everything between -180 and 180 degrees.
      if (inlonLL.gt.180.0_8) lonLL=inlonLL-360.0_8
      if (inlonUR.gt.180.0_8) lonUR=inlonUR-360.0_8
      ! if the model domain wraps across the prime meridian add 360 to longitude
      if (inlonLL.gt.inlonUR)lonUR = inlonUR + 360.0_8
      latLL = inlatLL
      latUR = inlatUR

      dlat = latUR - latLL
      dlon = lonUR - lonLL
      cell_width = dlon/resolution
      cell_height = dlat/resolution
      minspace_x  = 3.0_8*cell_width
      minspace_y  = 3.0_8*cell_height

      nread   = 0
      ncities = 0

      CityMasterFile = 'world_cities.txt'
      inquire( file=trim(adjustl(CityMasterFile)), exist=IsThere )
      if(.not.IsThere)then
        write(6,*)"WARNING: Could not find file: ",trim(adjustl(CityMasterFile)),&
                  "in CWD."
        write(6,*)"         Trying in default install path: ",trim(Ash3dHome)
        CityMasterFile = trim(Ash3dHome) // &
                          DirDelim // 'share' // &
                          DirDelim // 'post_proc' // &
                          DirDelim // 'world_cities.txt'
        inquire( file=trim(adjustl(CityMasterFile)), exist=IsThere2 )
        if(.not.IsThere2)then
          write(6,*)"ERROR: Could not find file: ",trim(adjustl(CityMasterFile))
          write(6,*)"       Skipping cities"
          ncities = 0
          return
        endif
      endif
      open(unit=12,file=trim(adjustl(CityMasterFile)))

      read(12,*)                                     !skip the first line
      do while ((ncities.lt.maxcities).and.(iostatus.ge.0))
        read(12,'(a133)',IOSTAT=iostatus) inputline
        read(inputline,2) CityLon, CityLat, CityName
2       format(f16.4,f15.4,a26)
        if ((CityLon.gt.lonLL).and.(CityLon.lt.lonUR).and. &
            (CityLat.gt.latLL).and.(CityLat.lt.latUR)) then
          ! Make sure this city is not near any others
          IsOkay=.true.
          call space_checker(maxcities,CityLon_out,CityLat_out,CityName_out,ncities, & 
                             CityLon,CityLat, &
                             minspace_x,minspace_y,IsOkay)
          if (IsOkay) then
            ncities = ncities+1
            CityName_out(ncities) = CityName
            CityLon_out(ncities)  = CityLon
            CityLat_out(ncities)  = CityLat
          endif
          ! if the model domain crosses over the prime meridian
        elseif ((CityLon+360.0_8.gt.lonLL).and.(CityLon+360.0_8.lt.lonUR).and. &
                (CityLat.gt.latLL).and.(CityLat.lt.latUR)) then
          ! Make sure this city is not near any others
          IsOkay=.true.
          call space_checker(maxcities,CityLon_out,CityLat_out,CityName_out,ncities, & 
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
      close(12)

      if(outCode.gt.0)then
        if(ncities.gt.0) then
          open(unit=13,file='cities.xy')
          if(outCode.eq.1)then
            ! for gmt
            do i=1,ncities
              write(13,4) CityLon_out(i),CityLat_out(i),CityName_out(i)
4             format(2f10.4,'  10  0  9  BL    ',a26)
            enddo
          elseif(outCode.eq.2)then
            do i=1,ncities
              ! for gnuplot
              write(13,3) CityLon_out(i),CityLat_out(i),CityName_out(i)
              write(13,*) real(CityLon_out(i),kind=4),real(CityLat_out(i),kind=4),&
                          '"',trim(adjustl(CityName_out(i))),'"'
3             format(2f10.4,1x,a26)
            enddo
          else
            write(*,*)"outCode not recognized. No output file written"
          endif
          close(13)
        endif
      endif

!      write(6,*) 'wrote ',ncities,' to cities.xy'
!
!      close(13)

      end subroutine citylist
         
!***************************************************************************************

      subroutine space_checker(maxcities,CityLon_out,CityLat_out,CityName_out,ncities, & 
                                  CityLon,CityLat, &
                                  minspace_x,minspace_y,IsOkay)
      implicit none

      integer ,intent(in) :: maxcities
      integer            :: icity, ncities
      real(kind=8)       :: CityLat, CityLat_out(maxcities), CityLon, CityLon_out(maxcities)
      character(len=26)  :: CityName_out(maxcities)
      !character(len=1)   :: answer
      real(kind=8)       :: minspace_x, minspace_y
      logical            :: IsOkay                 !true if city is not near any others

!      write(6,*) 'compare with:'
      do icity=1,ncities
!        write(6,6) CityName_out(icity), CityLon_out(icity), CityLat_out(icity)
!6       format(a26,2f10.4)
        if ((abs(CityLon_out(icity)-CityLon).lt.minspace_x).and. &
            (abs(CityLat_out(icity)-CityLat).lt.minspace_y)) then
!          write(6,*) 'too close'
          IsOkay=.false.
          exit
        endif
      enddo

      !write(6,*) 'continue?'
      !read(5,'(a1)') answer
      !if (answer.eq.'n') stop

      return

      end subroutine
