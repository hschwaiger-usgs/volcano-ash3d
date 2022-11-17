!##############################################################################
!
!    write_2Dprof_PNG_plplot
!
!##############################################################################

      subroutine write_2Dprof_PNG_plplot(vprof_ID)

      use precis_param

      use global_param,  only : &
         KG_2_MG,KM3_2_M3

      use mesh,          only : &
         nzmax,z_cc_pd

      use Output_Vars,   only : &
         pr_ash,CLOUDCON_THRESH

      use time_data,     only : &
         ntmax,time_native

      use io_data,       only : &
         Site_vprofile,x_vprofile,y_vprofile,cdf_b3l1,VolcanoName

      use Source,        only : &
         neruptions,e_Volume,e_Duration,e_StartTime,e_PlumeHeight

      use time_data,     only : &
         os_time_log,BaseYear,useLeap

      use plplot
      use iso_c_binding, only: c_ptr, c_loc, c_f_pointer

      implicit none

      integer,intent(in) :: vprof_ID

      !integer,parameter :: NLEVEL     = 11
      integer,parameter :: NLEVEL     = 20
      integer,parameter :: NUM_AXES   = 1
      integer,parameter :: NUM_LABELS = 1

      character(len=14) :: dp_pngfile
      character(len=26) :: coord_str
      character(len=76) :: title_str
      character(len=80) :: outstring
      integer :: k,i
      integer :: ioerr,iw,iwf

      ! PLPLOT variables
      real(kind=plflt)  :: tmin
      real(kind=plflt)  :: tmax
      real(kind=plflt)  :: zmin
      real(kind=plflt)  :: zmax
      real(kind=plflt)  :: cmin
      real(kind=plflt)  :: cmax
      real(kind=plflt), dimension(:),   allocatable :: t, z
      real(kind=plflt), dimension(:,:), allocatable :: conc
      real(kind=plflt), dimension(:),   allocatable :: shedge
      integer           :: cont_color
      real(kind=plflt)  :: fill_width, cont_width
      real(kind=plflt)  :: colorbar_width, colorbar_height
      real(kind=plflt)  :: dy_newline
      !integer           :: NUM_AXES, NUM_LABELS
      !parameter(NUM_AXES=1, NUM_LABELS=1)
      character(len=20) :: axis_opts(NUM_AXES)
      integer           :: num_values(NUM_AXES)
      real(kind=plflt)  :: values(NUM_AXES,NLEVEL+1)
      real(kind=plflt)  :: axis_ticks(NUM_AXES)
      integer           :: axis_subticks(NUM_AXES)
      character(len=100):: labels(NUM_LABELS)
      integer           :: label_opts(NUM_LABELS)
      character(len=1)  :: defined
      real(kind=plflt)  :: cloudcon_thresh_mgm3

      real(kind=plflt)   :: tr(6)

      real(kind=plflt)   :: clevel(1)

      INTERFACE
        character (len=20) function HS_xmltime(HoursSince,byear,useLeaps)
          real(kind=8)              :: HoursSince
          integer                   :: byear
          logical                   :: useLeaps
        end function HS_xmltime
      END INTERFACE

      cloudcon_thresh_mgm3 = CLOUDCON_THRESH * KG_2_MG / KM3_2_M3 !convert from kg/km3 to mg/m3

      clevel(1) = cloudcon_thresh_mgm3

      write(dp_pngfile,54) vprof_ID,".png"
 54   format('plplt_',i4.4,a4)

      tmin=real(0,kind=plflt)
      tmax=real(ceiling(time_native(ntmax)),kind=plflt)
      zmin=real(0,kind=plflt)
      zmax=real(z_cc_pd(nzmax),kind=plflt)
      cmin=real(0,kind=plflt)
      cmax=real(maxval(pr_ash(:,:,vprof_ID)),kind=plflt)

      tr = (/ tmax/real(ntmax-1,kind=plflt), 0.0_plflt, 0.0_plflt, &
              0.0_plflt, zmax/real(nzmax-1,kind=plflt), 0.0_plflt /)

      allocate(t(ntmax))
      allocate(z(nzmax))
      allocate(conc(ntmax,nzmax))

      t = real(time_native(1:ntmax),kind=plflt)
      z = real(z_cc_pd(1:nzmax),kind=plflt)
      do i=1,ntmax
        do k=1,nzmax
          conc(i,k) = pr_ash(k,i,vprof_ID)
        enddo
      enddo

      allocate(shedge(NLEVEL+1))
      ! Here we linearly interpolate color levels to the min/max of the data
      do i = 1, NLEVEL+1
        shedge(i) = cmin + (cmax - cmin) * real(i-1,kind=plflt) / real(NLEVEL,kind=plflt)
      enddo
      ! Here we hard-wire the color levels to the same as the kml plot
      !shedge(:) = (/ 0.0_plflt,  0.1_plflt,   0.3_plflt,   1.0_plflt, 2.0_plflt, &
      !              10.0_plflt, 30.0_plflt, 100.0_plflt, 300.0_plflt, & 
      !              1000.0_plflt, 3000.0_plflt, 10000.0_plflt /)
      fill_width = 2
      cont_color = 0
      cont_width = 0
      axis_opts(1) = 'bcvtm'
      axis_ticks(1) = 0.0_plflt
      axis_subticks(1) = 0
      label_opts(1) = PL_COLORBAR_LABEL_RIGHT
      labels(1) = 'Ash conc. mg/m3'

      write(coord_str,101)x_vprofile(vprof_ID),y_vprofile(vprof_ID)
 101  format(' (lon=',f7.2,', lat=',f6.2,')')
      write(title_str,*)trim(adjustl(Site_vprofile(vprof_ID))),coord_str

      ! Set up for plplot
      call plsdev("pngcairo")      ! Set output device (png, pdf, etc.)
      call plsfnam ( dp_pngfile )  ! Set output filename

      call plsetopt("geometry","854x603, 854x603")  ! Set image size
      !call plsetopt("bg","FFFFFF")                  ! Set background color to white
      !-----------------------------------
      !call plspal0('cmap0_black_on_white.pal')
      !call plspal1('cmap1_gray.pal',1)
      !-----------------------------------
      !call plspal0('cmap0_black_on_white.pal')
      !call plspal1('cmap1_blue_red.pal',1)
      !call plscmap0n(3)
      !call plscmap1l
      !-----------------------------------
      call plspal0('cmap0_black_on_white.pal')
      call plspal1('cmap1_blue_yellow.pal',1)
      call plscmap0n(3)

      ! Initialize plplot
      call plinit()

        ! pladv: Advance the (sub-)page
      call pladv(0)
        ! plvpor: Specify viewport using normalized subpage coordinates
      call plvpor(0.12_plflt, 0.75_plflt, 0.3_plflt, 0.8_plflt)
        ! plwind: Specify window 
      call plwind(tmin,tmax,zmin,zmax)
        ! plpsty: Select area fill pattern (0 is for solid)
      call plpsty(0)

      ! Now plot the data
      ! Call used with plplot5.10.0
      call plshades(conc(:ntmax,:nzmax), defined, &
        tmin,tmax,zmin,zmax, &
        shedge, fill_width, &
        cont_color, cont_width )

      ! Smaller text, scale by second argument
      call plschr( 0.0_plflt, 0.5_plflt )
      ! Small ticks on the vertical axis
      call plsmaj( 0.0_plflt, 0.5_plflt )
      call plsmin( 0.0_plflt, 0.5_plflt )
      call plcol0(1)
      call plbox('bcnst', 0.0_plflt, 0, 'bcnstv', 0.0_plflt, 0)
      call plcol0(2)
      call plschr(0.0_plflt,1.0_plflt)  ! Change font scale
      call pllab("Time (hours after eruption)", "Height (km)", trim(adjustl(title_str)))

      num_values(1) = NLEVEL + 1;
      values(1,:)   = shedge;
      call plcolorbar( colorbar_width, colorbar_height, &  ! these are output values
            PL_COLORBAR_SHADE ,  &
            !ior(PL_COLORBAR_SHADE, PL_COLORBAR_SHADE_LABEL), & ! sets options
            0, &                                               ! sets position
            0.015_plflt, 0.1_plflt, 0.0375_plflt, 0.8_plflt, &  ! x,y, dimesions
            0, 1, 1, &                                         ! bg_color,bb_color,bb_style
            0.0_plflt, shedge(NLEVEL+1), & ! low/high 
            cont_color, cont_width, &
            label_opts, labels, &
            axis_opts, &
            axis_ticks, axis_subticks, &
            num_values, values )

      ! This plots a contour line at the threshold level
      !call plcont(conc,1,ntmax,1,nzmax,clevel, tr)

      ! Now add the annotation box
        ! pladv: Advance the (sub-)page
      call pladv(1)
        ! plvpor: Specify viewport using normalized subpage coordinates
      call plvpor(0.05_plflt, 0.35_plflt, 0.05_plflt, 0.2_plflt)
        ! plwind: Specify window 
      call plwind(0.0_plflt, 1.0_plflt, 0.0_plflt, 1.0_plflt )
        ! plbox: Draw a box with axes, etc
        !       (xopt, xtick, nxsub, yopt, ytick, nysub)
      !call plbox('bc', 0.0_plflt, 0, 'bc', 0.0_plflt, 0 )
        !  plschr: Set character size
        !          (def, scale)
      call plschr( 0.0_plflt, 0.7_plflt )
      dy_newline = 0.13_plflt
        ! plptex: Write text inside the viewport
      write(outstring,*)"Volcano: ",trim(adjustl(VolcanoName))
      call plptex(0.02_plflt, 1.0_plflt-1.0_plflt*dy_newline, 1.0_plflt, 0.0_plflt, &
                  0.0_plflt, outstring )
      write(outstring,*)"Run Date: ",trim(adjustl(os_time_log))
      call plptex(0.02_plflt, 1.0_plflt-2.0_plflt*dy_newline, 1.0_plflt, 0.0_plflt, &
                  0.0_plflt, outstring )
      read(cdf_b3l1,*,iostat=ioerr) iw,iwf
      write(outstring,*)"Windfile: ",iwf
      call plptex(0.02_plflt, 1.0_plflt-3.0_plflt*dy_newline, 1.0_plflt, 0.0_plflt, &
                  0.0_plflt, outstring )

        ! pladv: Advance the (sub-)page
      call pladv(2)
        ! plvpor: Specify viewport using normalized subpage coordinates
      call plvpor(0.4_plflt, 0.75_plflt, 0.05_plflt, 0.2_plflt)
        ! plwind: Specify window 
      call plwind(0.0_plflt, 1.0_plflt, 0.0_plflt, 1.0_plflt )
        ! plbox: Draw a box with axes, etc
        !       (xopt, xtick, nxsub, yopt, ytick, nysub)
      !call plbox('bc', 0.0_plflt, 0, 'bc', 0.0_plflt, 0 )

      write(outstring,*)"Erup. Start Time: ",HS_xmltime(e_StartTime(1),BaseYear,useLeap)
      call plptex(0.02_plflt, 1.0_plflt-1.0_plflt*dy_newline, 1.0_plflt, 0.0_plflt, &
                  0.0_plflt, outstring )
      write(outstring,111)e_PlumeHeight(1)
 111  format(' Erup. Plume Height: ',f7.2,' km')
      call plptex(0.02_plflt, 1.0_plflt-2.0_plflt*dy_newline, 1.0_plflt, 0.0_plflt, &
                  0.0_plflt, outstring )
      write(outstring,121)e_Duration(1)
 121  format(' Erup. Duration: ',f7.2,' hours')
      call plptex(0.02_plflt, 1.0_plflt-3.0_plflt*dy_newline, 1.0_plflt, 0.0_plflt, &
                  0.0_plflt, outstring )
      write(outstring,131)e_Volume(1)
 131  format(' Erup. Volume: ',f10.5,' km3 (DRE)')
      call plptex(0.02_plflt, 1.0_plflt-4.0_plflt*dy_newline, 1.0_plflt, 0.0_plflt, &
                  0.0_plflt, outstring )


      ! Close PLplot library
      call plend

      end subroutine write_2Dprof_PNG_plplot

!##############################################################################
!
!    write_DepPOI_TS_PNG_plplot
!
!##############################################################################

      subroutine write_DepPOI_TS_PNG_plplot(pt_indx)

      use precis_param

      use Airports,      only : &
         nairports,Airport_Code,Airport_Name,Airport_x,Airport_y,&
         Airport_Latitude,Airport_Longitude,Airport_Thickness_TS

      use io_data,       only : &
         nWriteTimes,WriteTimes,VolcanoName

      use time_data,     only : &
         Simtime_in_hours

      use plplot

      implicit none

      integer :: pt_indx,i

      real(kind=8) :: ymaxpl
      character(len=14) :: dp_pngfile
      integer,save      :: plot_index = 0

      real(kind=plflt) :: xmin
      real(kind=plflt) :: xmax
      real(kind=plflt) :: ymin
      real(kind=plflt) :: ymax
      real(kind=plflt), dimension(:), allocatable :: x, y, x0, y0
      integer(kind=4) :: r1 = 0
      integer(kind=4) :: g1 = 0
      integer(kind=4) :: b1 = 0
      integer(kind=4) :: r2 = 136
      integer(kind=4) :: g2 = 136
      integer(kind=4) :: b2 = 136

      if(Airport_Thickness_TS(pt_indx,nWriteTimes).lt.0.01_ip)then
        !write(*,*)"No deposit at ",pt_indx,Airport_Name(pt_indx)
        return
      else
        plot_index = plot_index + 1
        !write(*,*)"Processing ",pt_indx,plot_index,Airport_Name(pt_indx),'--'
      endif

      write(dp_pngfile,55) plot_index,".png"
 55   format('plplt_',i4.4,a4)

      if(Airport_Thickness_TS(plot_index,nWriteTimes).lt.0.01)then
        ymaxpl = 1.0
      elseif(Airport_Thickness_TS(plot_index,nWriteTimes).lt.1.0)then
        ymaxpl = 1.0
      elseif(Airport_Thickness_TS(plot_index,nWriteTimes).lt.5.0)then
        ymaxpl = 5.0
      elseif(Airport_Thickness_TS(plot_index,nWriteTimes).lt.25.0)then
        ymaxpl = 25.0
      else
        ymaxpl = 100.0
      endif

      xmin=real(0,kind=plflt)
      xmax=real(ceiling(Simtime_in_hours),kind=plflt)
      ymin=real(0,kind=plflt)
      ymax=real(ymaxpl,kind=plflt)

      allocate(x(nWriteTimes))
      allocate(y(nWriteTimes))
      allocate(x0(nWriteTimes+1))
      allocate(y0(nWriteTimes+1))

      x = real(WriteTimes,kind=plflt)
      y = real(Airport_Thickness_TS(pt_indx,1:nWriteTimes),kind=plflt)
      x0(1:nWriteTimes)=x(1:nWriteTimes)
      x0(nWriteTimes+1)=x(nWriteTimes)
      y0(1:nWriteTimes)=y(1:nWriteTimes)
      y0(nWriteTimes+1)=0.0_plflt

      ! Set up for plplot
      call plsdev("pngcairo")      ! Set output device (png, pdf, etc.)
      call plsfnam ( dp_pngfile )  ! Set output filename

      call plsetopt("geometry","400x300, 400x300")  ! Set image size
      call plsetopt("bg","FFFFFF")                  ! Set background color to white
      ! Initialize plplot
      call plinit()

      ! Create a labelled box to hold the plot.
      call plscolbg(r1,g1,b1) ! Set color back to black
      call plcol0(0)
      call plschr(0.0_plflt,1.7_plflt)  ! Chage font scale
      call plenv( xmin, xmax, ymin, ymax, 0, 0 )
      call pllab( "Time (hours after eruption)", "Deposit Thickeness (mm)", Airport_Name(pt_indx))

      call plscolbg(r2,g2,b2) ! Set pen color to grey
      call plcol0(0)
      call plline( x , y ) ! Simple line plot
      call plfill( x0(1:nWriteTimes+1), y0(1:nWriteTimes+1) )

      ! Close PLplot library
      call plend

      end subroutine write_DepPOI_TS_PNG_plplot
