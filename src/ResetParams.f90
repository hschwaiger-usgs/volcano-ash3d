!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Subroutine input_data_ResetParams()
!
!  Called from: 
!  Arguments:
!    none
!
!  This subroutine
!
!  This subroutine is called from Ash3d.F90 after Read_Control_File only if the
!  input file has a block with the keyword OPTMOD=RESETPARAMS.
!  An example block with all the variables availalbe to reset is given here. Not
!  all variables need to be listed; only those changed from the defaults shown below.
!
!***********************
!# Reset parameters
!***********************
!OPTMOD=RESETPARAMS
! MagmaDensity         = 3500.0
! DepositDensity       = 1300.0
! LAM_GS_THRESH        = 250.0
! AIRBORNE_THRESH      = 1.0e-3
! GRAV                 = 9.81
! RAD_EARTH            = 6371.229
! CFL                  = 0.80
! DT_MIN               = 1.0e-5
! DT_MAX               = 1.0
! ZPADDING             = 1.3
! DEPO_THRESH          = 1.0e-1
! DEPRATE_THRESH       = 1.0e-2
! CLOUDCON_THRESH      = 2.0e-1
! CLOUDCON_GRID_THRESH = 2.0e-1
! CLOUDLOAD_THRESH     = 1.0e-2
! THICKNESS_THRESH     = 1.0e-3
! DBZ_THRESH           = -2.0e+1
! VelMod_umb           = 1
! lambda_umb           = 0.2
! N_BV_umb             = 0.02
! k_entrainment_umb    = 0.1
! SuzK_umb             = 12.0
! useMoistureVars      = F
! useVz_rhoG           = T
! cdf_institution      = USGS
! cdf_run_class        = Analysis
! cdf_url              = https://vsc-ash.wr.usgs.gov/ash3d-gui
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine input_data_ResetParams

      use precis_param

      use io_units

      use global_param,  only : &
         nmods,GRAV,CFL,DT_MIN,DT_MAX,RAD_EARTH,EPS_SMALL,&
         useMoistureVars, useVz_rhoG

      use mesh,          only : &
         ZPADDING

      use io_data,       only : &
         infile,cdf_institution,cdf_run_class,cdf_url

      use Tephra,        only : &
         MagmaDensity,DepositDensity,LAM_GS_THRESH,AIRBORNE_THRESH

      use Output_Vars,   only : &
         DEPO_THRESH,DEPRATE_THRESH,CLOUDCON_THRESH,CLOUDLOAD_THRESH,&
         THICKNESS_THRESH,DBZ_THRESH,CLOUDCON_GRID_THRESH

      use Source_Umbrella, only : &
         k_entrainment_umb,lambda_umb,N_BV_umb,SuzK_umb ,&
         VelMod_umb

      integer, parameter :: MAXPARAMS = 50

      character(len=80) :: linebuffer080
      character :: testkey
      integer :: ios
      character(len=20) :: mod_name
      integer :: substr_pos
      integer :: iparam
      character(len=20),dimension(MAXPARAMS) :: pname
      real(kind=ip)    ,dimension(MAXPARAMS) :: pvalue
      character(len=50),dimension(MAXPARAMS) :: pvalue_str
      integer :: i

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"    Searching for OPTMOD=RESETPARAMS"
      endif;enddo
      nmods = 0
      open(unit=10,file=infile,status='old',err=1900)

      read(10,'(a80)',iostat=ios)linebuffer080
      do while(ios.eq.0)
        read(10,'(a80)',iostat=ios)linebuffer080
        substr_pos = index(linebuffer080,'OPTMOD')
        if(substr_pos.eq.1)then
          ! found an optional module
          !  Parse for the keyword
            write(*,*)"Found optmod"
          read(linebuffer080,1104)mod_name
          if(trim(adjustl(mod_name)).eq.'RESETPARAMS')then
            exit
          endif
        endif
1104    format(7x,a20)
      enddo

      !write(*,*)"Now parsing RESETPARAMS block"
      
      read(10,'(a80)',iostat=ios)linebuffer080
      read(linebuffer080,*)testkey
      iparam = 0
      do while(ios.eq.0.and. &
               testkey.ne.'#'.and.testkey.ne.'*')
        iparam = iparam + 1
        substr_pos = index(linebuffer080,'=')
        pname(iparam)=trim(adjustl(linebuffer080(1:substr_pos-1)))
        ! first try to read this parameter as a real value
        read(linebuffer080(substr_pos+1:80),*,iostat=ios)pvalue(iparam)
        if(ios.ne.0)then
          ! If reading a floating point value fails, then try to read as
          ! a string
          read(linebuffer080(substr_pos+1:50),*)pvalue_str(iparam)
        endif
        read(10,'(a80)',iostat=ios)linebuffer080
        read(linebuffer080,*)testkey
      enddo

      ! We've read all the parameters to reset, now loop through the list
      ! and reset values, error-check, etc.

      do i = 1,iparam
        if (pname(i).eq.'MagmaDensity') then
          ! error-checking
          if (pvalue(i).le.0.0_ip)then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: MagmaDensity must be > 0"
            endif;enddo
            stop 1
          elseif (pvalue(i).gt.10000.0_ip)then
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"WARNING: MagmaDensity seems high."
              write(outlog(io),*)"         Units should be kg/m3"
            endif;enddo
          endif
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"  Resetting MagmaDensity from ",MagmaDensity,&
                              "to ",pvalue(i)
          endif;enddo
          MagmaDensity = pvalue(i)
        elseif (pname(i).eq.'DepositDensity') then
          ! error-checking
          if (pvalue(i).le.0.0_ip)then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: DepositDensity must be > 0"
            endif;enddo
            stop 1
          elseif (pvalue(i).gt.10000.0_ip)then
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"WARNING: DepositDensity seems high."
              write(outlog(io),*)"         Units should be kg/m2"
            endif;enddo
          endif
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"  Resetting DepositDensity from ",DepositDensity,&
                              "to ",pvalue(i)
          endif;enddo 
          DepositDensity = pvalue(i)
        elseif (pname(i).eq.'LAM_GS_THRESH') then
          ! error-checking
          if (pvalue(i).le.0.0_ip)then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: LAM_GS_THRESH must be > 0"
            endif;enddo
            stop 1
          elseif (pvalue(i).gt.1000.0_ip)then
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"WARNING: LAM_GS_THRESH seems high."
              write(outlog(io),*)"         Units should be m"
            endif;enddo
          endif
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"  Resetting LAM_GS_THRESH from ",LAM_GS_THRESH,&
                              "to ",pvalue(i)
          endif;enddo
          LAM_GS_THRESH = pvalue(i)
        elseif (pname(i).eq.'AIRBORNE_THRESH') then
          ! error-checking
          if (pvalue(i).le.0.0_ip)then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: AIRBORNE_THRESH must be > 0"
            endif;enddo
            stop 1
          elseif (pvalue(i).gt.1.0_ip)then
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"WARNING: AIRBORNE_THRESH seems high."
              write(outlog(io),*)"         Units should be kg"
            endif;enddo
          endif
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"  Resetting AIRBORNE_THRESH from ",AIRBORNE_THRESH,&
                              "to ",pvalue(i)
          endif;enddo
          AIRBORNE_THRESH = pvalue(i)
        elseif (pname(i).eq.'GRAV') then
          ! error-checking
          if (pvalue(i).le.0.0_ip)then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: GRAV must be > 0"
            endif;enddo
            stop 1
          elseif (pvalue(i).gt.20.0_ip)then
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"WARNING: GRAV seems high."
              write(outlog(io),*)"         Units should be m/s2"
            endif;enddo
          endif
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"  Resetting GRAV from ",GRAV,&
                              "to ",pvalue(i)
          endif;enddo
          GRAV = pvalue(i)
        elseif (pname(i).eq.'RAD_EARTH') then
          ! error-checking
          if (pvalue(i).le.0.0_ip)then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: RAD_EARTH must be > 0"
            endif;enddo
            stop 1
          elseif (pvalue(i).gt.10000.0_ip)then
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"WARNING: RAD_EARTH seems high."
              write(outlog(io),*)"         Units should be km"
            endif;enddo
          endif
          do io=1,2;if(VB(io).le.verbosity_info)then          
            write(outlog(io),*)"  Resetting RAD_EARTH from ",RAD_EARTH,&
                              "to ",pvalue(i)
          endif;enddo
          RAD_EARTH = pvalue(i)
        elseif (pname(i).eq.'CFL') then
          ! error-checking
          if (pvalue(i).le.0.0_ip)then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: CFL must be > 0"
            endif;enddo
            stop 1
          elseif (pvalue(i).ge.1.0_ip)then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: CFL must be < 1"
            endif;enddo
            stop 1
          endif
          do io=1,2;if(VB(io).le.verbosity_info)then            
            write(outlog(io),*)"  Resetting CFL from ",CFL,&
                              "to ",pvalue(i)
          endif;enddo 
          CFL = pvalue(i)
        elseif (pname(i).eq.'DT_MIN') then
          ! error-checking
          if (pvalue(i).le.0.0_ip)then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: DT_MIN must be > 0"
            endif;enddo
            stop 1
          elseif (pvalue(i).gt.1.0_ip)then
            do io=1,2;if(VB(io).le.verbosity_info)then            
              write(outlog(io),*)"WARNING: DT_MIN seems high."
              write(outlog(io),*)"         Units should be hours"
            endif;enddo
          endif
          do io=1,2;if(VB(io).le.verbosity_info)then         
            write(outlog(io),*)"  Resetting DT_MIN from ",DT_MIN,&
                              "to ",pvalue(i)
          endif;enddo 
          DT_MIN = pvalue(i)

        elseif (pname(i).eq.'DT_MAX') then
          ! error-checking
          if (pvalue(i).le.0.0_ip)then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: DT_MAX must be > 0"
            endif;enddo
            stop 1
          elseif (pvalue(i).le.DT_MIN)then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: DT_MAX must be > DT_MIN"
            endif;enddo
            stop 1
          elseif (pvalue(i).gt.10.0_ip)then
            do io=1,2;if(VB(io).le.verbosity_info)then         
              write(outlog(io),*)"WARNING: DT_MAX seems high."
              write(outlog(io),*)"         Units should be hours"
            endif;enddo
          endif
          do io=1,2;if(VB(io).le.verbosity_info)then         
            write(outlog(io),*)"  Resetting DT_MAX from ",DT_MAX,&
                              "to ",pvalue(i)
          endif;enddo
          DT_MAX = pvalue(i)
        elseif (pname(i).eq.'ZPADDING') then
          ! error-checking
          if (pvalue(i).le.1.0_ip)then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: ZPADDING must be > 1"
            endif;enddo
            stop 1
          elseif (pvalue(i).gt.5.0_ip)then
            do io=1,2;if(VB(io).le.verbosity_info)then         
              write(outlog(io),*)"WARNING: ZPADDING seems high."
              write(outlog(io),*)"         This is the factor times the plume height (~1.3)"
            endif;enddo
          endif
          do io=1,2;if(VB(io).le.verbosity_info)then         
            write(outlog(io),*)"  Resetting ZPADDING from ",ZPADDING,&
                              "to ",pvalue(i)
          endif;enddo
          ZPADDING = pvalue(i)
        elseif (pname(i).eq.'DEPO_THRESH') then
          ! error-checking
          if (pvalue(i).le.0.0_ip)then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: DEPO_THRESH must be > 0"
            endif;enddo
            stop 1
          elseif (pvalue(i).gt.5.0_ip)then
            do io=1,2;if(VB(io).le.verbosity_info)then         
              write(outlog(io),*)"WARNING: DEPO_THRESH seems high."
              write(outlog(io),*)&
                "         This is the threshold to track deposit (mm)"
            endif;enddo
          endif
          do io=1,2;if(VB(io).le.verbosity_info)then         
            write(outlog(io),*)"  Resetting DEPO_THRESH from ",DEPO_THRESH,&
                              "to ",pvalue(i)
          endif;enddo
          DEPO_THRESH = pvalue(i)
        elseif (pname(i).eq.'DEPRATE_THRESH') then
          ! error-checking
          if (pvalue(i).le.0.0_ip)then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: DEPRATE_THRESH must be > 0"
            endif;enddo
            stop 1
          elseif (pvalue(i).gt.1.0_ip)then
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"WARNING: seems high."
              write(outlog(io),*)&
                "         This is the threshold to track deposit rate (mm/hr)"
            endif;enddo
          endif
          do io=1,2;if(VB(io).le.verbosity_info)then         
            write(outlog(io),*)"  Resetting DEPRATE_THRESH from ",DEPRATE_THRESH,&
                              "to ",pvalue(i)
          endif;enddo
          DEPRATE_THRESH = pvalue(i)
        elseif (pname(i).eq.'CLOUDCON_THRESH') then
          ! error-checking
          if (pvalue(i).le.0.0_ip)then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: CLOUDCON_THRESH must be > 0"
            endif;enddo
            stop 1
          elseif (pvalue(i).gt.5.0_ip)then
            do io=1,2;if(VB(io).le.verbosity_info)then         
              write(outlog(io),*)"WARNING: CLOUDCON_THRESH seems high."
              write(outlog(io),*)&
                "         This is the threshold to track cloud concentration (t/km3)"
            endif;enddo
          endif
          do io=1,2;if(VB(io).le.verbosity_info)then         
            write(outlog(io),*)"  Resetting CLOUDCON_THRESH from ",CLOUDCON_THRESH,&
                              "to ",pvalue(i)
          endif;enddo
          CLOUDCON_THRESH = pvalue(i)
        elseif (pname(i).eq.'CLOUDCON_GRID_THRESH') then
          ! error-checking
          if (pvalue(i).le.0.0_ip)then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: CLOUDCON_GRID_THRESH must be > 0"
            endif;enddo
            stop 1
          elseif (pvalue(i).gt.1.0e-1_ip)then
            do io=1,2;if(VB(io).le.verbosity_info)then         
              write(outlog(io),*)"WARNING: CLOUDCON_GRID_THRESH seems high."
              write(outlog(io),*)&
                "         This is the concentration threshold to identify regions to calculate (t/km3)"
             endif;enddo
          endif
          do io=1,2;if(VB(io).le.verbosity_info)then         
            write(outlog(io),*)"  Resetting CLOUDCON_GRID_THRESH from ",CLOUDCON_GRID_THRESH,&
                              "to ",pvalue(i)
          endif;enddo
          CLOUDCON_GRID_THRESH = pvalue(i)
        elseif (pname(i).eq.'CLOUDLOAD_THRESH') then
          ! error-checking
          if (pvalue(i).le.0.0_ip)then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: CLOUDLOAD_THRESH must be > 0"
            endif;enddo
            stop 1
          elseif (pvalue(i).gt.5.0_ip)then
            do io=1,2;if(VB(io).le.verbosity_info)then         
              write(outlog(io),*)"WARNING: CLOUDLOAD_THRESH seems high."
              write(outlog(io),*)&
                "         This is the threshold to track cloud load (t/km2)"
            endif;enddo
          endif
          do io=1,2;if(VB(io).le.verbosity_info)then         
            write(outlog(io),*)"  Resetting CLOUDLOAD_THRESH from ",CLOUDLOAD_THRESH,&
                              "to ",pvalue(i)
          endif;enddo
          CLOUDLOAD_THRESH = pvalue(i)
        elseif (pname(i).eq.'THICKNESS_THRESH') then
          ! error-checking
          if (pvalue(i).le.0.0_ip)then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: THICKNESS_THRESH must be > 0"
            endif;enddo
            stop 1
          elseif (pvalue(i).gt.5.0_ip)then
            do io=1,2;if(VB(io).le.verbosity_info)then         
              write(outlog(io),*)"WARNING: THICKNESS_THRESH seems high."
              write(outlog(io),*)&
                "         This is the threshold to track deposit (mm)"
            endif;enddo
          endif
          do io=1,2;if(VB(io).le.verbosity_info)then         
            write(outlog(io),*)"  Resetting THICKNESS_THRESH from ",THICKNESS_THRESH,&
                              "to ",pvalue(i)
          endif;enddo
          THICKNESS_THRESH = pvalue(i)
        elseif (pname(i).eq.'DBZ_THRESH') then
          ! error-checking
          if (pvalue(i).ge.0.0_ip)then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: DBZ_THRESH must be < 0"
            endif;enddo
            stop 1
          elseif (pvalue(i).lt.-1000.0_ip)then
            do io=1,2;if(VB(io).le.verbosity_info)then         
              write(outlog(io),*)"WARNING: DBZ_THRESH seems low."
              write(outlog(io),*)&
                "         This is the threshold to reflectivity (db)"
            endif;enddo
          endif
          do io=1,2;if(VB(io).le.verbosity_info)then         
            write(outlog(io),*)"  Resetting DBZ_THRESH from ",DBZ_THRESH,&
                              "to ",pvalue(i)
          endif;enddo
          DBZ_THRESH = pvalue(i)
        elseif (pname(i).eq.'VelMod_umb') then
          if(abs(pvalue(i)-1.0_ip).lt.EPS_SMALL)then
            VelMod_umb = 1
          elseif(abs(pvalue(i)-2.0_ip).lt.EPS_SMALL)then
            VelMod_umb = 2
          else
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: VelMod_umb not recognized; must be 1 or 2"
            endif;enddo
          endif
        elseif (pname(i).eq.'lambda_umb') then
          ! error-checking
          if (pvalue(i).le.0.0_ip)then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: lambda_umb must be > 0"
            endif;enddo
            stop 1
          elseif (pvalue(i).gt.100.0_ip)then
            do io=1,2;if(VB(io).le.verbosity_info)then         
              write(outlog(io),*)"WARNING: lambda_umb seems high."
              write(outlog(io),*)&
                "         This is the Shape factor for the umbrella cloud"
            endif;enddo
          endif
          do io=1,2;if(VB(io).le.verbosity_info)then         
            write(outlog(io),*)"  Resetting lambda_umb from ",lambda_umb,&
                              "to ",pvalue(i)
          endif;enddo
          lambda_umb = pvalue(i)
        elseif (pname(i).eq.'N_BV_umb') then
          ! error-checking
          if (pvalue(i).le.0.0_ip)then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR:  must be > 0"
            endif;enddo
            stop 1
          elseif (pvalue(i).gt.10.0_ip)then
            do io=1,2;if(VB(io).le.verbosity_info)then         
              write(outlog(io),*)"WARNING: N_BV_umb seems high."
              write(outlog(io),*)&
                "         This is the Brunt-Vaisala frequency (1/s) for umbrella sources"
            endif;enddo
          endif
          do io=1,2;if(VB(io).le.verbosity_info)then         
            write(outlog(io),*)"  Resetting N_BV_umb from ",N_BV_umb,&
                              "to ",pvalue(i)
          endif;enddo
          N_BV_umb = pvalue(i)
        elseif (pname(i).eq.'k_entrainment_umb') then
          ! error-checking
          if (pvalue(i).le.0.0_ip)then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: k_entrainment_umb must be > 0"
            endif;enddo
            stop 1
          elseif (pvalue(i).gt.5.0_ip)then
            do io=1,2;if(VB(io).le.verbosity_info)then         
              write(outlog(io),*)"WARNING: k_entrainment_umb seems high."
              write(outlog(io),*)&
                "         This is the umbrella entrainment coefficient"
            endif;enddo
          endif
          do io=1,2;if(VB(io).le.verbosity_info)then         
            write(outlog(io),*)"  Resetting k_entrainment_umb from ",k_entrainment_umb,&
                              "to ",pvalue(i)
          endif;enddo
          k_entrainment_umb = pvalue(i)
        elseif (pname(i).eq.'SuzK_umb') then
          ! error-checking
          if (pvalue(i).le.0.0_ip)then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: SuzK_umb must be > 0"
            endif;enddo
            stop 1
          elseif (pvalue(i).gt.20.0_ip)then
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"WARNING: SuzK_umb seems high."
              write(outlog(io),*)&
                "         This is the Suzuki parameter for the umbrella"
            endif;enddo
          endif
          do io=1,2;if(VB(io).le.verbosity_info)then         
            write(outlog(io),*)"  Resetting SuzK_umb from ",SuzK_umb,&
                              "to ",pvalue(i)
          endif;enddo
          SuzK_umb = pvalue(i)
        elseif (pname(i).eq.'useMoistureVars') then
          ! error-checking
          ! This should either be 0 for .false. or 1 for .true.
          ! but we read the value as a float
          if(abs(pvalue(i)).lt.EPS_SMALL)then
            useMoistureVars = .false.
          elseif(abs(pvalue(i)-1.0_ip).lt.EPS_SMALL)then
            useMoistureVars = .true.
          else
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: useMoistureVars should be:"
              write(errlog(io),*)"         0 for false"
              write(errlog(io),*)"         1 for true"
              write(errlog(io),*)"     Value read = ",pvalue(i)
              stop 1
            endif;enddo
          endif
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"  Resetting useMoistureVars to ",useMoistureVars
          endif;enddo
        elseif (pname(i).eq.'useVz_rhoG') then
          ! error-checking
          ! This should either be 0 for .false. or 1 for .true.
          ! but we read the value as a float
          if(abs(pvalue(i)).lt.EPS_SMALL)then
            useVz_rhoG = .false.
          elseif(abs(pvalue(i)-1.0_ip).lt.EPS_SMALL)then
            useVz_rhoG = .true.
          else
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: useVz_rhoG should be:"
              write(errlog(io),*)"         0 for false"
              write(errlog(io),*)"         1 for true"
              write(errlog(io),*)"     Value read = ",pvalue(i)
              stop 1
            endif;enddo
          endif
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"  Resetting useVz_rhoG to ",useVz_rhoG
          endif;enddo
        elseif (pname(i).eq.'cdf_institution') then
          cdf_institution = trim(adjustl(pvalue_str(i)))
        elseif (pname(i).eq.'cdf_run_class') then
          if(abs(pvalue(i)-1.0_ip).lt.EPS_SMALL)then
            cdf_url = 'Analysis'
          elseif(abs(pvalue(i)-2.0_ip).lt.EPS_SMALL)then
            cdf_url = 'Hypothetical'
          elseif(abs(pvalue(i)-3.0_ip).lt.EPS_SMALL)then
            cdf_url = 'Forecast'
          else
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: run_class not recognized; must be 1,2, or 3"
            endif;enddo
          endif
        elseif (pname(i).eq.'cdf_url') then
          cdf_url = trim(adjustl(pvalue_str(i)))
        else
          do io=1,2;if(VB(io).le.verbosity_info)then         
            write(outlog(io),*)"Found unknown parameter/value: ", &
                                pname(i),pvalue(i)
            write(outlog(io),*)"No action taken."
          endif;enddo
        endif
      enddo

      close(10)

      return

1900  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)  'error: cannot find input file: ',infile
        write(errlog(io),*)  'Program stopped'
      endif;enddo
      stop 1

      end subroutine input_data_ResetParams

