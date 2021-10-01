!******************************************************************************

      subroutine input_data_ResetParams

      use precis_param

      use io_units

      use global_param,  only : &
         nmods,GRAV,CFL,DT_MIN,DT_MAX,RAD_EARTH

      use mesh,          only : &
         ZPADDING

      use io_data,       only : &
         infile

      use Tephra,        only : &
         MagmaDensity,DepositDensity,LAM_GS_THRESH,AIRBORNE_THRESH

      use Output_Vars,   only : &
         DEPO_THRESH,DEPRATE_THRESH,CLOUDCON_THRESH,CLOUDLOAD_THRESH,&
         THICKNESS_THRESH,DBZ_THRESH

      integer, parameter :: MAXPARAMS = 50

      character(len=80)  :: linebuffer
      character :: testkey
      integer :: ios
      character(len=20) :: mod_name
      integer :: substr_pos
      integer :: iparam
      character(len=20),dimension(MAXPARAMS) :: pname
      real(kind=ip),dimension(MAXPARAMS)     :: pvalue
      integer :: i

      write(global_info,*)"    Searching for OPTMOD=RESETPARAMS"
      nmods = 0
      open(unit=10,file=infile,status='old',err=1900)

      read(10,'(a80)',iostat=ios)linebuffer
      do while(ios.eq.0)
        read(10,'(a80)',iostat=ios)linebuffer
        substr_pos = index(linebuffer,'OPTMOD')
        if(substr_pos.eq.1)then
          ! found an optional module
          !  Parse for the keyword
            write(*,*)"Found optmod"
          read(linebuffer,1104)mod_name
          if(adjustl(trim(mod_name)).eq.'RESETPARAMS')then
            exit
          endif
        endif
1104    format(7x,a20)
      enddo

      write(*,*)"Now parsing RESETPARAMS block"
      
      read(10,'(a80)',iostat=ios)linebuffer
      read(linebuffer,*)testkey
      iparam = 0
      do while(ios.eq.0.and. &
               testkey.ne.'#'.and.testkey.ne.'*')
        iparam = iparam + 1
        substr_pos = index(linebuffer,'=')
        pname(iparam)=adjustl(trim(linebuffer(1:substr_pos-1)))
        read(linebuffer(substr_pos+1:80),*)pvalue(iparam)
        read(10,'(a80)',iostat=ios)linebuffer
        read(linebuffer,*)testkey
      enddo

      ! We've read all the parameters to reset, now loop through the list
      ! and reset values, error-check, etc.

      do i = 1,iparam
        if (pname(i).eq.'MagmaDensity') then
          ! error-checking
          if (pvalue(i).le.0.0_ip)then
            write(global_error,*)"ERROR: MagmaDensity must be > 0"
            stop 0
          elseif (pvalue(i).gt.10000.0_ip)then
            write(global_error,*)"WARNING: MagmaDensity seems high."
            write(global_error,*)"         Units should be kg/m3"
          endif
          write(global_info,*)"  Resetting MagmaDensity from ",MagmaDensity,&
                              "to ",pvalue(i)
          MagmaDensity = pvalue(i)
        elseif (pname(i).eq.'DepositDensity') then
          ! error-checking
          if (pvalue(i).le.0.0_ip)then
            write(global_error,*)"ERROR: DepositDensity must be > 0"
            stop 0
          elseif (pvalue(i).gt.10000.0_ip)then
            write(global_error,*)"WARNING: DepositDensity seems high."
            write(global_error,*)"         Units should be kg/m2"
          endif
          write(global_info,*)"  Resetting DepositDensity from ",DepositDensity,&
                              "to ",pvalue(i)
          DepositDensity = pvalue(i)
        elseif (pname(i).eq.'LAM_GS_THRESH') then
          ! error-checking
          if (pvalue(i).le.0.0_ip)then
            write(global_error,*)"ERROR: LAM_GS_THRESH must be > 0"
            stop 0
          elseif (pvalue(i).gt.1000.0_ip)then
            write(global_error,*)"WARNING: LAM_GS_THRESH seems high."
            write(global_error,*)"         Units should be m"
          endif
          write(global_info,*)"  Resetting LAM_GS_THRESH from ",LAM_GS_THRESH,&
                              "to ",pvalue(i)
          LAM_GS_THRESH = pvalue(i)
        elseif (pname(i).eq.'AIRBORNE_THRESH') then
          ! error-checking
          if (pvalue(i).le.0.0_ip)then
            write(global_error,*)"ERROR: AIRBORNE_THRESH must be > 0"
            stop 0
          elseif (pvalue(i).gt.1.0_ip)then
            write(global_error,*)"WARNING: AIRBORNE_THRESH seems high."
            write(global_error,*)"         Units should be kg"
          endif
          write(global_info,*)"  Resetting AIRBORNE_THRESH from ",AIRBORNE_THRESH,&
                              "to ",pvalue(i)
          AIRBORNE_THRESH = pvalue(i)
        elseif (pname(i).eq.'GRAV') then
          ! error-checking
          if (pvalue(i).le.0.0_ip)then
            write(global_error,*)"ERROR: GRAV must be > 0"
            stop 0
          elseif (pvalue(i).gt.20.0_ip)then
            write(global_error,*)"WARNING: GRAV seems high."
            write(global_error,*)"         Units should be m/s2"
          endif
          write(global_info,*)"  Resetting GRAV from ",GRAV,&
                              "to ",pvalue(i)
          GRAV = pvalue(i)
        elseif (pname(i).eq.'RAD_EARTH') then
          ! error-checking
          if (pvalue(i).le.0.0_ip)then
            write(global_error,*)"ERROR: RAD_EARTH must be > 0"
            stop 0
          elseif (pvalue(i).gt.10000.0_ip)then
            write(global_error,*)"WARNING: RAD_EARTH seems high."
            write(global_error,*)"         Units should be km"
          endif
          write(global_info,*)"  Resetting RAD_EARTH from ",RAD_EARTH,&
                              "to ",pvalue(i)
          RAD_EARTH = pvalue(i)
        elseif (pname(i).eq.'CFL') then

          ! error-checking
          if (pvalue(i).le.0.0_ip)then
            write(global_error,*)"ERROR: CFL must be > 0"
            stop 0
          elseif (pvalue(i).ge.1.0_ip)then
            write(global_error,*)"ERROR: CFL must be < 1"
            stop 0
          endif
          write(global_info,*)"  Resetting CFL from ",CFL,&
                              "to ",pvalue(i)
          CFL = pvalue(i)
        elseif (pname(i).eq.'DT_MIN') then
          ! error-checking
          if (pvalue(i).le.0.0_ip)then
            write(global_error,*)"ERROR: DT_MIN must be > 0"
            stop 0
          elseif (pvalue(i).gt.1.0_ip)then
            write(global_error,*)"WARNING: DT_MIN seems high."
            write(global_error,*)"         Units should be hours"
          endif
          write(global_info,*)"  Resetting DT_MIN from ",DT_MIN,&
                              "to ",pvalue(i)
          DT_MIN = pvalue(i)

        elseif (pname(i).eq.'DT_MAX') then
          ! error-checking
          if (pvalue(i).le.0.0_ip)then
            write(global_error,*)"ERROR: DT_MAX must be > 0"
            stop 0
          elseif (pvalue(i).le.DT_MIN)then
            write(global_error,*)"ERROR: DT_MAX must be > DT_MIN"
            stop 0
          elseif (pvalue(i).gt.10.0_ip)then
            write(global_error,*)"WARNING: DT_MAX seems high."
            write(global_error,*)"         Units should be hours"
          endif
          write(global_info,*)"  Resetting DT_MAX from ",DT_MAX,&
                              "to ",pvalue(i)
          DT_MAX = pvalue(i)
        elseif (pname(i).eq.'ZPADDING') then
          ! error-checking
          if (pvalue(i).le.1.0_ip)then
            write(global_error,*)"ERROR: ZPADDING must be > 1"
            stop 0
          elseif (pvalue(i).gt.5.0_ip)then
            write(global_error,*)"WARNING: ZPADDING seems high."
            write(global_error,*)"         This is the factor times the plume height (~1.3)"
          endif
          write(global_info,*)"  Resetting ZPADDING from ",ZPADDING,&
                              "to ",pvalue(i)
          ZPADDING = pvalue(i)
        elseif (pname(i).eq.'DEPO_THRESH') then
          ! error-checking
          if (pvalue(i).le.0.0_ip)then
            write(global_error,*)"ERROR: DEPO_THRESH must be > 0"
            stop 0
          elseif (pvalue(i).gt.5.0_ip)then
            write(global_error,*)"WARNING: DEPO_THRESH seems high."
            write(global_error,*)&
              "         This is the threshold to track deposit (mm)"
          endif
          write(global_info,*)"  Resetting DEPO_THRESH from ",DEPO_THRESH,&
                              "to ",pvalue(i)
          DEPO_THRESH = pvalue(i)
        elseif (pname(i).eq.'DEPRATE_THRESH') then
          ! error-checking
          if (pvalue(i).le.0.0_ip)then
            write(global_error,*)"ERROR: DEPRATE_THRESH must be > 0"
            stop 0
          elseif (pvalue(i).gt.1.0_ip)then
            write(global_error,*)"WARNING: seems high."
            write(global_error,*)&
              "         This is the threshold to track deposit rate (mm/hr)"
          endif
          write(global_info,*)"  Resetting DEPRATE_THRESH from ",DEPRATE_THRESH,&
                              "to ",pvalue(i)
          DEPRATE_THRESH = pvalue(i)
        elseif (pname(i).eq.'CLOUDCON_THRESH') then
          ! error-checking
          if (pvalue(i).le.0.0_ip)then
            write(global_error,*)"ERROR: CLOUDCON_THRESH must be > 0"
            stop 0
          elseif (pvalue(i).gt.5.0_ip)then
            write(global_error,*)"WARNING: CLOUDCON_THRESH seems high."
            write(global_error,*)&
              "         This is the threshold to track cloud concentration (t/km3)"
          endif
          write(global_info,*)"  Resetting CLOUDCON_THRESH from ",CLOUDCON_THRESH,&
                              "to ",pvalue(i)
          CLOUDCON_THRESH = pvalue(i)
        elseif (pname(i).eq.'CLOUDCON_GRID_THRESH') then
          ! error-checking
          if (pvalue(i).le.0.0_ip)then
            write(global_error,*)"ERROR: CLOUDCON_GRID_THRESH must be > 0"
            stop 0
          elseif (pvalue(i).gt.1.0e-1_ip)then
            write(global_error,*)"WARNING: CLOUDCON_GRID_THRESH seems high."
            write(global_error,*)&
              "         This is the concentration threshold to identify regions to calculate (t/km3)"
          endif
          write(global_info,*)"  Resetting CLOUDCON_GRID_THRESH from ",CLOUDCON_GRID_THRESH,&
                              "to ",pvalue(i)
          CLOUDCON_GRID_THRESH = pvalue(i)
        elseif (pname(i).eq.'CLOUDLOAD_THRESH') then
          ! error-checking
          if (pvalue(i).le.0.0_ip)then
            write(global_error,*)"ERROR: CLOUDLOAD_THRESH must be > 0"
            stop 0
          elseif (pvalue(i).gt.5.0_ip)then
            write(global_error,*)"WARNING: CLOUDLOAD_THRESH seems high."
            write(global_error,*)&
              "         This is the threshold to track cloud load (t/km2)"
          endif
          write(global_info,*)"  Resetting CLOUDLOAD_THRESH from ",CLOUDLOAD_THRESH,&
                              "to ",pvalue(i)
          CLOUDLOAD_THRESH = pvalue(i)
        elseif (pname(i).eq.'THICKNESS_THRESH') then
          ! error-checking
          if (pvalue(i).le.0.0_ip)then
            write(global_error,*)"ERROR: THICKNESS_THRESH must be > 0"
            stop 0
          elseif (pvalue(i).gt.5.0_ip)then
            write(global_error,*)"WARNING: THICKNESS_THRESH seems high."
            write(global_error,*)&
              "         This is the threshold to track deposit (mm)"
          endif
          write(global_info,*)"  Resetting THICKNESS_THRESH from ",THICKNESS_THRESH,&
                              "to ",pvalue(i)
          THICKNESS_THRESH = pvalue(i)

        elseif (pname(i).eq.'DBZ_THRESH') then
          ! error-checking
          if (pvalue(i).ge.0.0_ip)then
            write(global_error,*)"ERROR: DBZ_THRESH must be < 0"
            stop 0
          elseif (pvalue(i).lt.-1000.0_ip)then
            write(global_error,*)"WARNING: DBZ_THRESH seems low."
            write(global_error,*)&
              "         This is the threshold to reflectivity (db)"
          endif
          write(global_info,*)"  Resetting DBZ_THRESH from ",DBZ_THRESH,&
                              "to ",pvalue(i)
          DBZ_THRESH = pvalue(i)
        elseif (pname(i).eq.'lambda') then
          ! error-checking
          if (pvalue(i).le.0.0_ip)then
            write(global_error,*)"ERROR: lambda must be > 0"
            stop 0
          elseif (pvalue(i).gt.100.0_ip)then
            write(global_error,*)"WARNING: lambda seems high."
            write(global_error,*)&
              "         This is the Suzuki parameter for the umbrella"
          endif
          !write(global_info,*)"  Resetting lambda from ",lambda,&
          !                    "to ",pvalue(i)
          !lambda = pvalue(i)
        elseif (pname(i).eq.'N_BV') then
          ! error-checking
          if (pvalue(i).le.0.0_ip)then
            write(global_error,*)"ERROR:  must be > 0"
            stop 0
          elseif (pvalue(i).gt.10.0_ip)then
            write(global_error,*)"WARNING: N_BV seems high."
            write(global_error,*)&
              "         This is the Brunt-Vaisala frequency (1/s)"
          endif
          !write(global_info,*)"  Resetting N_BV from ",N_BV,&
          !                    "to ",pvalue(i)
          !N_BV = pvalue(i)
        elseif (pname(i).eq.'k_entrainment') then
          ! error-checking
          if (pvalue(i).le.0.0_ip)then
            write(global_error,*)"ERROR: k_entrainment must be > 0"
            stop 0
          elseif (pvalue(i).gt.5.0_ip)then
            write(global_error,*)"WARNING: k_entrainment seems high."
            write(global_error,*)&
              "         This is the umbrella entrainment coefficient"
          endif
          !write(global_info,*)"  Resetting k_entrainment from ",k_entrainment,&
          !                    "to ",pvalue(i)
          !k_entrainment = pvalue(i)
        else
          write(global_info,*)"Found unknown parameter/value: ", &
                              pname(i),pvalue(i)
          write(global_info,*)"No action taken."
        endif
      enddo

      close(10)

      return

1900  write(global_info,*)  'error: cannot find input file: ',infile
      write(global_info,*)  'Program stopped'
      write(global_log,*)  'error: cannot find input file: ',infile
      write(global_log,*)  'Program stopped'
      stop 1

      end subroutine input_data_ResetParams

