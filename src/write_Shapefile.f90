
      subroutine write_ShapeFile_Polyline(iprod,itime)

      ! input should be the contours, the output product, and the time step

      use precis_param

      use global_param,  only : &
         IsLitEnd,IsLinux,IsWindows,IsMacOS

!      use global_param,  only : &
!         IsLinux,IsWindows,IsMacOS,DirPrefix,DirDelim

      use Output_Vars,   only : &
         ContourDataX,ContourDataY,ContourDataNcurves,ContourDataNpoints,&
         ContourLev,nConLev

      use io_data,       only : &
         WriteTimes,VolcanoName

      use io_data,       only : &
         cdf_institution,cdf_run_class,cdf_url

      implicit none

      integer,intent(in) :: iprod
      integer,intent(in) :: itime

      character(len=8)  :: ov_fileroot = "testfile"
      character(len=12) :: ov_mainfile
      character(len=12) :: ov_indxfile
      character(len=12) :: ov_dbasfile
      character(len=12) :: ov_projfile
      character(len=12) :: ov_zipfile
      character(len=4)  :: ov_mainext = ".shp"
      character(len=4)  :: ov_indxext = ".shx"
      character(len=4)  :: ov_dbasext = ".dbf"
      character(len=4)  :: ov_projext = ".prj"
      character(len=4)  :: ov_zipext  = ".zip"
      integer           :: ov_mainID  = 22
      integer           :: ov_indxID  = 23
      integer           :: ov_dbasID  = 24
      integer           :: ov_projID  = 25

      character(len=40) :: title_plot
      character(len=40) :: plot_variable
      character(len=15) :: plot_units

      logical           :: debugmode = .false.

      integer                                  :: nrec      ! num of records (e.g. contour levels)
      integer(kind=4),dimension(:),allocatable :: NumParts
      integer(kind=4),dimension(:),allocatable :: NumPoints ! total points
      integer :: i,irec,ipart,ipt
      integer,dimension(:),allocatable :: reclen
      integer :: offset

      integer(kind=4) :: file_code
      integer(kind=4) :: tmp4
      integer(kind=4) :: file_length
      integer(kind=4) :: version
      integer(kind=4) :: shape_type
      real(kind=8),dimension(:),allocatable :: xmin
      real(kind=8),dimension(:),allocatable :: xmax
      real(kind=8),dimension(:),allocatable :: ymin
      real(kind=8),dimension(:),allocatable :: ymax
      real(kind=8)    :: zmin
      real(kind=8)    :: zmax
      real(kind=8)    :: mmin
      real(kind=8)    :: mmax

      ! Variables needed for the dBASE file
      integer(kind=1)  :: DBASE_zero     = 0
      integer(kind=1)  :: DBASE_v
      integer(kind=1)  :: DBASE_yy
      integer(kind=1)  :: DBASE_mm
      integer(kind=1)  :: DBASE_dd
      integer(kind=4)  :: DBASE_nrec      ! 1
      integer(kind=2)  :: DBASE_headlen   ! 129
      integer(kind=2)  :: DBASE_reclen    ! 115
      integer(kind=1)  :: DBASE_transflag ! 0
      integer(kind=1)  :: DBASE_cryptflag ! 0
      integer(kind=1)  :: DBASE_mdxflag   ! 0
      integer(kind=1)  :: DBASE_langID    ! 87

      character(len=11):: DBASE_FieldName ! 'name   '
      character(len=1) :: DBASE_FieldTyp
      integer(kind=1)  :: DBASE_FieldLen
      integer(kind=1)  :: DBASE_FieldDesTerm
      integer          :: fldlen
      character(len=1) :: DBASE_RecStart = ' '
      integer(kind=1)  :: DBASE_EOF      = 26

      integer          :: nattr          = 15
      character(len=10):: DBASE_TableRecData01  ! Organizaion
      character(len=42):: DBASE_TableRecData02  ! Volcano
      character(len=20):: DBASE_TableRecData03  ! Run date
      character(len= 5):: DBASE_TableRecData04  ! iwindformat
      character(len=20):: DBASE_TableRecData05  ! Run class
      character(len=20):: DBASE_TableRecData06  ! Erup. Start Time
      character(len=20):: DBASE_TableRecData07  ! Erup. Plume Height
      character(len=20):: DBASE_TableRecData08  ! Erup. Duration
      character(len=20):: DBASE_TableRecData09  ! Erup. Vol
      character(len=80):: DBASE_TableRecData10  ! URL
      character(len=24):: DBASE_TableRecData11  ! Variable
      character(len=24):: DBASE_TableRecData12  ! Value of contour level
      character(len=10):: DBASE_TableRecData13  ! unit for level
      character(len=10):: DBASE_TableRecData14  ! index
      character(len=20):: DBASE_TableRecData15  ! Time of data

      logical                 :: IsThere
      character(len=71):: zipcom

      INTERFACE
        subroutine writeShapFileFieldDesArr(ov_dbasID,fldlen,DBASE_FieldName,&
                                            DBASE_FieldTyp,DBASE_FieldLen)
          integer,intent(in)               :: ov_dbasID
          integer,intent(in)               :: fldlen
          character(len=fldlen),intent(in) :: DBASE_FieldName
          character(len=1),intent(in)      :: DBASE_FieldTyp
          integer(kind=1),intent(in)       :: DBASE_FieldLen
        end subroutine writeShapFileFieldDesArr
        integer(kind=2) function BigEnd_2int(isLit,r)
          logical         :: isLit
          integer(kind=2) :: r
        end function BigEnd_2int
        integer(kind=4) function BigEnd_4int(isLit,r)
          logical         :: isLit
          integer(kind=4) :: r
        end function BigEnd_4int
        integer(kind=2) function LitEnd_2int(isLit,r)
          logical         :: isLit
          integer(kind=2) :: r
        end function LitEnd_2int
        integer(kind=4) function LitEnd_4int(isLit,r)
          logical         :: isLit
          integer(kind=4) :: r
        end function LitEnd_4int
        real(kind=8) function LitEnd_8real(isLit,r)
          logical         :: isLit
          real(kind=8)    :: r
        end function LitEnd_8real
      END INTERFACE

      if(iprod.eq.3)then       ! deposit at specified times (mm)
        write(title_plot,'(a22,f5.2,a6)')': Deposit Thickness t=',WriteTimes(itime),' hours'
        plot_variable = 'Deposit Thickness'
        ov_fileroot = 'depothik'
        plot_units = 'mm'
      elseif(iprod.eq.4)then   ! deposit at specified times (inches)
        write(title_plot,'(a22,f5.2,a6)')': Deposit Thickness t=',WriteTimes(itime),' hours'
        plot_variable = 'Deposit Thickness'
        ov_fileroot = 'depothik'
        plot_units  = 'in'
      elseif(iprod.eq.5)then       ! deposit at final time (mm)
        title_plot = ': Final Deposit Thickness'
        plot_variable = 'Final Deposit Thickness'
        ov_fileroot = 'depothik'
        plot_units = 'mm'
      elseif(iprod.eq.6)then   ! deposit at final time (inches)
        title_plot = ': Final Deposit Thickness'
        plot_variable = 'Final Deposit Thickness'
        ov_fileroot = 'depothik'
        plot_units = 'in'
      elseif(iprod.eq.7)then   ! ashfall arrival time (hours)
        write(title_plot,'(a22)')': Ashfall arrival time'
        plot_variable = 'Ashfall arrival time'
        ov_fileroot = 'DepAvlTm'
        plot_units = 'hours'
      elseif(iprod.eq.8)then   ! ashfall arrival at airports/POI (mm)
        write(*,*)"ERROR: No map PNG output option for airport arrival time data."
        write(*,*)"       Should not be in write_2Dmap_PNG_dislin"
        stop 1
      elseif(iprod.eq.9)then   ! ash-cloud concentration
        write(title_plot,'(a28,f5.2,a6)')': Ash-cloud concentration t=',WriteTimes(itime),' hours'
        plot_variable ='Ash-cloud concentration'
        ov_fileroot = 'AshCdCon'
        plot_units = 'mg/m3'
      elseif(iprod.eq.10)then   ! ash-cloud height
        write(title_plot,'(a21,f5.2,a6)')': Ash-cloud height t=',WriteTimes(itime),' hours'
        plot_variable ='Ash-cloud height'
        ov_fileroot = 'AshCdHgt'
        plot_units = 'km'
      elseif(iprod.eq.11)then   ! ash-cloud bottom
        write(title_plot,'(a21,f5.2,a6)')': Ash-cloud bottom t=',WriteTimes(itime),' hours'
        plot_variable ='Ash-cloud bottom'
        ov_fileroot = 'AshCdBot'
        plot_units = 'km'
      elseif(iprod.eq.12)then   ! ash-cloud load
        write(title_plot,'(a19,f5.2,a6)')': Ash-cloud load t=',WriteTimes(itime),' hours'
        plot_variable ='Ash-cloud load'
        ov_fileroot = 'AshCdLod'
        plot_units = 'T/km2'
      elseif(iprod.eq.13)then  ! radar reflectivity
        write(title_plot,'(a26,f5.2,a6)')': Ash-cloud radar refl. t=',WriteTimes(itime),' hours'
        plot_variable ='Ash-cloud radar refl.'
        ov_fileroot = 'AshClRad'
        plot_units = 'dBz'
      elseif(iprod.eq.14)then   ! ashcloud arrival time (hours)
        write(title_plot,'(a24)')': Ash-cloud arrival time'
        plot_variable ='Ash-cloud arrival time'
        ov_fileroot = 'AshAvlTm'
        plot_units = 'hours'
      elseif(iprod.eq.15)then   ! topography
        write(title_plot,'(a12)')': Topography'
        plot_variable ='Topography'
        ov_fileroot = 'Topogrph'
        plot_units = 'km'
      elseif(iprod.eq.16)then   ! profile plots
        write(*,*)"ERROR: No map PNG output option for vertical profile data."
        write(*,*)"       Should not be in write_2Dmap_PNG_dislin"
        stop 1
      else
        write(*,*)"ERROR: unexpected variable"
        stop 1
      endif

      ! Get number of records; number of actual contour levels used
      nrec=0
      do i = 1,nConLev
        if(ContourDataNcurves(i).gt.0)nrec=nrec+1
      enddo

      ! We allocate variables here to hold the contour data mainly so ensure that
      ! the variables we write to the shapefile have exactly the kind values expected.
      allocate(NumParts(nrec))
      allocate(NumPoints(nrec))
      allocate(reclen(nrec))
      allocate(xmin(nrec)); xmin(:) =  1.0e8_8
      allocate(xmax(nrec)); xmax(:) = -1.0e8_8
      allocate(ymin(nrec)); ymin(:) =  1.0e8_8
      allocate(ymax(nrec)); ymax(:) = -1.0e8_8

      do irec = 1,nrec  ! This is the loop of the layers (records)
        NumParts(irec) = ContourDataNcurves(irec)
        NumPoints(irec) = 0  ! Initialize the number of points for this record
        ! Loop through all the parts (individual curves) for this level and sum points
        do ipart = 1,ContourDataNcurves(irec)
          NumPoints(irec) = NumPoints(irec) + ContourDataNpoints(irec,ipart)
        enddo
        ! Set the min/max for each record, as well as the global min/max and
        ! the length of each polyline record
        do ipart = 1,ContourDataNcurves(irec)
          do ipt = 1,ContourDataNpoints(irec,ipart)
            if(ContourDataX(irec,ipart,ipt).lt.xmin(irec))then
              xmin(irec) = ContourDataX(irec,ipart,ipt)
            endif
            if(ContourDataX(irec,ipart,ipt).gt.xmax(irec))then
              xmax(irec) = ContourDataX(irec,ipart,ipt)
            endif
            if(ContourDataY(irec,ipart,ipt).lt.ymin(irec))then
              ymin(irec) = ContourDataY(irec,ipart,ipt)
            endif
            if(ContourDataY(irec,ipart,ipt).gt.ymax(irec))then
              ymax(irec) = ContourDataY(irec,ipart,ipt)
            endif
          enddo
        enddo

        ! Calculate the lengths of each record
        ! Polyline content header length
        reclen(irec) = 1*4 + & ! Shape Type : 1-Integer : Little
                       4*8 + & ! Box        : 4-Double  : Little
                       1*4 + & ! NumParts   : 1-Integer : Little
                       1*4 + & ! NumPoints  : 1-Integer : Little
                       NumParts(irec)*4 + &  ! Parts      : NumParts-Integers : Little
                       NumPoints(irec)*2*8  ! Points     : 2-Double          : Little
        reclen(irec) = reclen(irec) / 2     ! total record length in 16-bit words
      enddo
      zmin = 0.0_8
      zmax = 0.0_8
      mmin = 0.0_8
      mmax = 0.0_8

      ! Note: all data in shapefiles are either
      !        integer : Signed 32-bit integer (4 bytes)
      !        double  : Signed 64-bit IEEE double-precision floating point number (8 bytes)
      ov_mainfile = trim(adjustl(ov_fileroot)) // ov_mainext
      ov_indxfile = trim(adjustl(ov_fileroot)) // ov_indxext
      ov_dbasfile = trim(adjustl(ov_fileroot)) // ov_dbasext
      ov_projfile = trim(adjustl(ov_fileroot)) // ov_projext
      ov_zipfile  = trim(adjustl(ov_fileroot)) // ov_zipext

      write(*,*)"Writing shapefile"

      ! First, writing the main .shp file which contains the contour (polyline) data.
      ! The specification has a mess of little-endian and big-endian writes with non-integer
      ! record offsets (e.g. 8-byte reals written to addres byte-36), so we open the file with
      ! access='stream' and use the functions [Lit,Big}End_[4int,8real] to write the correct
      ! bits.  For all the details, see
      ! https://www.esri.com/content/dam/esrisites/sitecore-archive/Files/Pdfs/library/whitepapers/pdfs/shapefile.pdf
      !  All writes of real*8 will use transfer() to write the correct bits in 4-byte packets
      open(ov_mainID, file=trim(adjustl(ov_mainfile)), access='stream', form='unformatted', status='replace')
      ! Note: The index file (.shx) has the same head as the .shp file (except the file length field)
      !       So duplicate all the writes
      open(ov_indxID, file=trim(adjustl(ov_indxfile)), access='stream', form='unformatted', status='replace')

      ! File header is 100 bytes
      file_code = 9994_4
      tmp4      = 0_4
      ! shp file
      write(ov_mainID)BigEnd_4int(IsLitEnd,file_code)   !  1-  4 : FH Byte 0 File Code 9994 Integer Big
      write(ov_mainID)BigEnd_4int(IsLitEnd,tmp4)        !  5-  8 : FH Byte 4 Unused 0 Integer Big
      write(ov_mainID)BigEnd_4int(IsLitEnd,tmp4)        !  9- 12 : FH Byte 8 Unused 0 Integer Big
      write(ov_mainID)BigEnd_4int(IsLitEnd,tmp4)        ! 13- 16 : FH Byte 12 Unused 0 Integer Big
      write(ov_mainID)BigEnd_4int(IsLitEnd,tmp4)        ! 17- 20 : FH Byte 16 Unused 0 Integer Big
      write(ov_mainID)BigEnd_4int(IsLitEnd,tmp4)        ! 21- 24 : FH Byte 20 Unused 0 Integer Big
      ! shx file
      write(ov_indxID)BigEnd_4int(IsLitEnd,file_code)
      write(ov_indxID)BigEnd_4int(IsLitEnd,tmp4)
      write(ov_indxID)BigEnd_4int(IsLitEnd,tmp4)
      write(ov_indxID)BigEnd_4int(IsLitEnd,tmp4)
      write(ov_indxID)BigEnd_4int(IsLitEnd,tmp4)
      write(ov_indxID)BigEnd_4int(IsLitEnd,tmp4)

      if(debugmode)then
        write(*,*)file_code,   " 1-  4 : FH Byte 0 File Code 9994 Integer Big"
        write(*,*)tmp4,        " 5-  8 : FH Byte 4 Unused 0 Integer Big"
        write(*,*)tmp4,        " 9- 12 : FH Byte 8 Unused 0 Integer Big"
        write(*,*)tmp4,        "13- 16 : FH Byte 12 Unused 0 Integer Big"
        write(*,*)tmp4,        "17- 20 : FH Byte 16 Unused 0 Integer Big"
        write(*,*)tmp4,        "21- 24 : FH Byte 20 Unused 0 Integer Big"
      endif

      !              -- File Header length
      !              |     -- Record header length for all records
      !              |     |        -- Sum of all length of polyline record contents
      !              |     |        |    
      file_length = 50 + (4*nrec) + sum(reclen(1:nrec))  ! total file length in 16-bit words
      write(ov_mainID)BigEnd_4int(IsLitEnd,file_length)          ! 25-28 : FH Byte 24 File Length Integer Big
      if(debugmode)write(*,*)file_length,"25-28 : FH Byte 24 File Length Integer Big"
      ! Recalculating length for the index file
      file_length = 50 + (4*nrec)
      write(ov_indxID)BigEnd_4int(IsLitEnd,file_length)

      ! Here's where we start writing as little-endian
      version = 1000
      write(ov_mainID)LitEnd_4int(IsLitEnd,version)             ! 29-32 : FH Byte 28 Version 1000 Integer Little
      if(debugmode)write(*,*)version,"29-32 : FH Byte 28 Version 1000 Integer Little"
      write(ov_indxID)LitEnd_4int(IsLitEnd,version)
      shape_type = 3
      write(ov_mainID)LitEnd_4int(IsLitEnd,shape_type)          ! 33-36 : FH Shape Type Integer Little
      if(debugmode)write(*,*)shape_type,"33-36 : FH Shape Type Integer Little"
      write(ov_indxID)LitEnd_4int(IsLitEnd,shape_type)          !         (use code 3 for polyline)

      ! Start of real values
      write(ov_mainID)LitEnd_8real(IsLitEnd,minval(xmin(1:nrec)))    ! 37- 44 : FH Bounding Box Xmin Double Little
      write(ov_mainID)LitEnd_8real(IsLitEnd,minval(ymin(1:nrec)))    ! 45- 52 : FH Bounding Box Ymin Double Little
      write(ov_mainID)LitEnd_8real(IsLitEnd,maxval(xmax(1:nrec)))    ! 53- 60 : FH Bounding Box Xmax Double Little
      write(ov_mainID)LitEnd_8real(IsLitEnd,maxval(ymax(1:nrec)))    ! 61- 68 : FH Bounding Box Ymax Double Little
      write(ov_indxID)LitEnd_8real(IsLitEnd,minval(xmin(1:nrec)))
      write(ov_indxID)LitEnd_8real(IsLitEnd,minval(ymin(1:nrec)))
      write(ov_indxID)LitEnd_8real(IsLitEnd,maxval(xmax(1:nrec)))
      write(ov_indxID)LitEnd_8real(IsLitEnd,maxval(ymax(1:nrec)))
      if(debugmode)then
        write(*,*)minval(xmin(1:nrec)),"37- 44 : FH Bounding Box Xmin Double Little"
        write(*,*)minval(ymin(1:nrec)),"45- 52 : FH Bounding Box Ymin Double Little"
        write(*,*)maxval(xmax(1:nrec)),"53- 60 : FH Bounding Box Xmax Double Little"
        write(*,*)maxval(ymax(1:nrec)),"61- 68 : FH Bounding Box Ymax Double Little"
      endif

      write(ov_mainID)LitEnd_8real(IsLitEnd,zmin)    ! 69- 76 : FH Bounding Box Zmin Double Little
      write(ov_mainID)LitEnd_8real(IsLitEnd,zmax)    ! 77- 84 : FH Bounding Box Zmax Double Little
      write(ov_mainID)LitEnd_8real(IsLitEnd,mmin)    ! 85- 92 : FH Bounding Box Mmin Double Little
      write(ov_mainID)LitEnd_8real(IsLitEnd,mmax)    ! 94-100 : FH Bounding Box Mmax Double Little
      write(ov_indxID)LitEnd_8real(IsLitEnd,zmin)
      write(ov_indxID)LitEnd_8real(IsLitEnd,zmax)
      write(ov_indxID)LitEnd_8real(IsLitEnd,mmin)
      write(ov_indxID)LitEnd_8real(IsLitEnd,mmax)

      if(debugmode)then
        write(*,*)zmin,"69- 76 : FH Bounding Box Zmin Double Little"
        write(*,*)zmax,"77- 84 : FH Bounding Box Zmax Double Little"
        write(*,*)mmin,"85- 92 : FH Bounding Box Mmin Double Little"
        write(*,*)mmax,"94-100 : FH Bounding Box Mmax Double Little"
      endif
      ! End of Main File Header

      ! Loop through all the records
      offset = 50
      do irec=1,nrec
        ! Each record contains
        !  Byte 0 Record Number Record Number Integer Big
        !  Byte 4 Content Length Content Length Integer Big
        !   Note: content length is the number of 16-bit words
        write(ov_mainID)BigEnd_4int(IsLitEnd,irec)
        if(debugmode)write(*,*)irec,"Byte 0 Record Number Record Number Integer Big"
        write(ov_mainID)BigEnd_4int(IsLitEnd,reclen(irec))
        if(debugmode)write(*,*)reclen(irec),"Byte 4 Content Length Content Length Integer Big"

        !  PolyLine Record Contents
        !  Byte 0 Shape Type 3 Integer 1 Little
        !  Byte 4 Box Box Double 4 Little
        !  Byte 36 NumParts NumParts Integer 1 Little
        !  Byte 40 NumPoints NumPoints Integer 1 Little
        !  Byte 44 Parts Parts Integer NumParts Little
        !  Byte X Points Points Point NumPoints Little
        write(ov_mainID)LitEnd_4int(IsLitEnd,shape_type)
        write(ov_mainID)LitEnd_8real(IsLitEnd,xmin(irec))
        write(ov_mainID)LitEnd_8real(IsLitEnd,ymin(irec))
        write(ov_mainID)LitEnd_8real(IsLitEnd,xmax(irec))
        write(ov_mainID)LitEnd_8real(IsLitEnd,ymax(irec))
        write(ov_mainID)LitEnd_4int(IsLitEnd,NumParts(irec))
        write(ov_mainID)LitEnd_4int(IsLitEnd,NumPoints(irec))

        if(debugmode)then
          write(*,*)shape_type,"Byte 0 Shape Type 3 Integer 1 Little"
          write(*,*)xmin(irec),"Byte 4 xmin Double Little"
          write(*,*)ymin(irec),"Byte 12 ymin Double Little"
          write(*,*)xmax(irec),"Byte 20 xmax Double Little"
          write(*,*)ymax(irec),"Byte 28 ymax Double Little"
          write(*,*)NumParts(irec),"Byte 36 NumParts NumParts Integer 1 Little"
          write(*,*)NumPoints(irec),"Byte 40 NumPoints NumPoints Integer 1 Little"
        endif
        ! Address of the start of the first part is 0
        write(ov_mainID)LitEnd_4int(IsLitEnd,0_4) ! An array of length NumParts with each value
                                                  ! the address (zero-offset) of the start of the part
        do i=1,NumParts(irec)-1
          tmp4 = sum(ContourDataNpoints(irec,1:i)) !-1
          write(ov_mainID)LitEnd_4int(IsLitEnd,tmp4)
        enddo

        do ipart=1,ContourDataNcurves(irec)
          do i=1,ContourDataNpoints(irec,ipart)
            write(ov_mainID)LitEnd_8real(IsLitEnd,ContourDataX(irec,ipart,i))
            write(ov_mainID)LitEnd_8real(IsLitEnd,ContourDataY(irec,ipart,i))
          enddo ! points in part
        enddo ! ipart

        ! The offset of the first record is just 50, but for subsequent offsets, we need
        ! to add the length of the previous record (as well as the record header).
        if(irec.ne.1)offset = offset + 4 + reclen(irec-1)
        write(ov_indxID)BigEnd_4int(IsLitEnd,offset)
        write(ov_indxID)BigEnd_4int(IsLitEnd,reclen(irec))

      enddo
      ! Close shp main file
      close(ov_mainID)
      close(ov_indxID)

      ! Now write the dbase file
      ! Full for dbase 7 specifications at
      ! https://www.dbase.com/Knowledgebase/INT/db7_file_fmt.htm
      ! For dbase 5 see
      ! https://www.oocities.org/geoff_wass/dBASE/GaryWhite/dBASE/FAQ/qformt.htm Sec. D
      ! Shapefile additional requirements:
      !  1. name consistent with shp, but with dbf extension
      !  2. table must contain one record per shape feature
      !  3. The record order must be the same as the order of shape features in the main file
      !  4. The year value in the dBASE header must be the year since 1900.
      ! Note: We will write dBASE V – MS-Windows (Level 5) format since we know that works.
      open(ov_dbasID, file=trim(adjustl(ov_dbasfile)), access='stream', form='unformatted', status='replace')

      ! Populate each of the TableRecData fields with dummy values so we can get lengths to put
      ! in the header
      ! HFS KLUDGE
      cdf_institution="USGS"
      cdf_run_class="Analysis"
      cdf_url="https://vsc-ash.wr.usgs.gov/ash3d-gui"

      write(DBASE_TableRecData01,*)trim(adjustl(cdf_institution))
      write(DBASE_TableRecData02,*)trim(adjustl(VolcanoName))
      write(DBASE_TableRecData03,'(a20)')'1800-01-01T00:00:00Z'
      write(DBASE_TableRecData04,'(i2)')25
      write(DBASE_TableRecData05,*)trim(adjustl(cdf_run_class))
      write(DBASE_TableRecData06,'(a20)')'1800-01-01T00:00:00Z'  ! start time
      write(DBASE_TableRecData07,'(f10.3)')40.0 ! plume height
      write(DBASE_TableRecData08,'(f10.3)')24.0 ! duration
      write(DBASE_TableRecData09,'(e10.3)')40.0 ! erup.vol
      write(DBASE_TableRecData10,*)trim(adjustl(cdf_url))
      write(DBASE_TableRecData11,*)trim(adjustl(plot_variable))
      write(DBASE_TableRecData12,'(a24)')"       0.000000000000000"
      write(DBASE_TableRecData13,*)trim(adjustl(plot_units))
      write(DBASE_TableRecData14,'(a10)')"         0"
      write(DBASE_TableRecData15,'(a20)')'1800-01-01T00:00:00Z'

      DBASE_zero = 0
      DBASE_v  = 3
      DBASE_yy =95
      DBASE_mm =07
      DBASE_dd =26

      DBASE_nrec      = nrec
      DBASE_headlen   = int(32,kind=2) +     &   ! Table File Header length
                        int(nattr*32,kind=2) &   ! length of field descriptor (attributes)
                        + int(1,kind=2)
      DBASE_reclen    = len(DBASE_TableRecData01) + &
                        len(DBASE_TableRecData02) + &
                        len(DBASE_TableRecData03) + &
                        len(DBASE_TableRecData04) + &
                        len(DBASE_TableRecData05) + &
                        len(DBASE_TableRecData06) + &
                        len(DBASE_TableRecData07) + &
                        len(DBASE_TableRecData08) + &
                        len(DBASE_TableRecData09) + &
                        len(DBASE_TableRecData10) + &
                        len(DBASE_TableRecData11) + &
                        len(DBASE_TableRecData12) + &
                        len(DBASE_TableRecData13) + &
                        len(DBASE_TableRecData14) + &
                        len(DBASE_TableRecData15) &
                        + int(1,kind=2)
      DBASE_transflag = 0
      DBASE_cryptflag = 0
      DBASE_mdxflag   = 0
      DBASE_langID    = 87

      DBASE_FieldDesTerm = 13

      ! 0   :1 byte   : version number
      write(ov_dbasID)DBASE_v
      ! 1-3 :3 bytes  : Date of last update; in YYMMDD format. 
      write(ov_dbasID)DBASE_yy
      write(ov_dbasID)DBASE_mm
      write(ov_dbasID)DBASE_dd
      !4-7  : 4 bytes : Number of records in the table. 
      write(ov_dbasID)LitEnd_4int(IsLitEnd,DBASE_nrec)
      !8-9  : 2 bytes : Number of bytes in the header. 
      write(ov_dbasID)LitEnd_2int(IsLitEnd,DBASE_headlen)
      !10-11: 2 bytes : Number of bytes in the record.
      write(ov_dbasID)LitEnd_2int(IsLitEnd,DBASE_reclen)
      !12-13: 2 bytes : Reserved; filled with zeros.
      write(ov_dbasID)DBASE_zero
      write(ov_dbasID)DBASE_zero
      !14   : 1 byte  : Flag indicating incomplete dBASE transaction.
      write(ov_dbasID)DBASE_transflag
      !15   : 1 byte  : Encryption flag.
      write(ov_dbasID)DBASE_cryptflag
      !16-27: 12 bytes: Reserved for multi-user processing.
      do i=1,12
        write(ov_dbasID)DBASE_zero
      enddo
      !28   : 1 byte  : Production MDX flag; x00 stored in this byte if a no .mdx exists
      write(ov_dbasID)DBASE_mdxflag
      !29   : 1-byte  : Language driver ID
      write(ov_dbasID)DBASE_langID
      !30-31: 1 byte  : Reserved; filled with zeros.
      write(ov_dbasID)DBASE_zero
      write(ov_dbasID)DBASE_zero

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !  FIELDS
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      DBASE_FieldName ='ORG'
      fldlen=len(trim(adjustl(DBASE_FieldName)))
      DBASE_FieldTyp = 'C'
      DBASE_FieldLen=len(DBASE_TableRecData01)
      call writeShapFileFieldDesArr(ov_dbasID,fldlen,DBASE_FieldName,&
                                    DBASE_FieldTyp,DBASE_FieldLen)
      DBASE_FieldName ='VOLC'
      fldlen=len(trim(adjustl(DBASE_FieldName)))
      DBASE_FieldTyp = 'C'
      DBASE_FieldLen=len(DBASE_TableRecData02)
      call writeShapFileFieldDesArr(ov_dbasID,fldlen,DBASE_FieldName,&
                                    DBASE_FieldTyp,DBASE_FieldLen)
      DBASE_FieldName ='RUN DATE'
      fldlen=len(trim(adjustl(DBASE_FieldName)))
      DBASE_FieldTyp = 'C'
      DBASE_FieldLen=len(DBASE_TableRecData03)
      call writeShapFileFieldDesArr(ov_dbasID,fldlen,DBASE_FieldName,&
                                    DBASE_FieldTyp,DBASE_FieldLen)
      DBASE_FieldName ='WINDFRMT'
      fldlen=len(trim(adjustl(DBASE_FieldName)))
      DBASE_FieldTyp = 'C'
      DBASE_FieldLen=len(DBASE_TableRecData04)
      call writeShapFileFieldDesArr(ov_dbasID,fldlen,DBASE_FieldName,&
                                    DBASE_FieldTyp,DBASE_FieldLen)
      DBASE_FieldName ='RUN CLASS'
      fldlen=len(trim(adjustl(DBASE_FieldName)))
      DBASE_FieldTyp = 'C'
      DBASE_FieldLen=len(DBASE_TableRecData05)
      call writeShapFileFieldDesArr(ov_dbasID,fldlen,DBASE_FieldName,&
                                    DBASE_FieldTyp,DBASE_FieldLen)
      DBASE_FieldName ='E_STIME'  ! 1800-01-01T00:00:00Z
      fldlen=len(trim(adjustl(DBASE_FieldName)))
      DBASE_FieldTyp = 'C'
      DBASE_FieldLen=len(DBASE_TableRecData06)
      call writeShapFileFieldDesArr(ov_dbasID,fldlen,DBASE_FieldName,&
                                    DBASE_FieldTyp,DBASE_FieldLen)
      DBASE_FieldName ='E_PLMH'   ! 40.0 km
      fldlen=len(trim(adjustl(DBASE_FieldName)))
      DBASE_FieldTyp = 'C'
      DBASE_FieldLen=len(DBASE_TableRecData07)
      call writeShapFileFieldDesArr(ov_dbasID,fldlen,DBASE_FieldName,&
                                    DBASE_FieldTyp,DBASE_FieldLen)
      DBASE_FieldName ='E_DUR'     ! 24.0 hours
      fldlen=len(trim(adjustl(DBASE_FieldName)))
      DBASE_FieldTyp = 'C'
      DBASE_FieldLen=len(DBASE_TableRecData08)
      call writeShapFileFieldDesArr(ov_dbasID,fldlen,DBASE_FieldName,&
                                    DBASE_FieldTyp,DBASE_FieldLen)
      DBASE_FieldName ='E_VOL'
      fldlen=len(trim(adjustl(DBASE_FieldName)))
      DBASE_FieldTyp = 'C'
      DBASE_FieldLen=len(DBASE_TableRecData09)
      call writeShapFileFieldDesArr(ov_dbasID,fldlen,DBASE_FieldName,&
                                    DBASE_FieldTyp,DBASE_FieldLen)
      DBASE_FieldName ='URL'
      fldlen=len(trim(adjustl(DBASE_FieldName)))
      DBASE_FieldTyp = 'C'
      DBASE_FieldLen=len(DBASE_TableRecData10)
      call writeShapFileFieldDesArr(ov_dbasID,fldlen,DBASE_FieldName,&
                                    DBASE_FieldTyp,DBASE_FieldLen)
      DBASE_FieldName ='VAR'      ! variable
      fldlen=len(trim(adjustl(DBASE_FieldName)))
      DBASE_FieldTyp = 'C'
      DBASE_FieldLen=len(DBASE_TableRecData11)
      call writeShapFileFieldDesArr(ov_dbasID,fldlen,DBASE_FieldName,&
                                    DBASE_FieldTyp,DBASE_FieldLen)
      DBASE_FieldName ='VALUE'    ! contour level
      fldlen=len(trim(adjustl(DBASE_FieldName)))
      DBASE_FieldTyp = 'N'
      DBASE_FieldLen=len(DBASE_TableRecData12)
      call writeShapFileFieldDesArr(ov_dbasID,fldlen,DBASE_FieldName,&
                                    DBASE_FieldTyp,DBASE_FieldLen)
      DBASE_FieldName ='UNITS'    ! units for level
      fldlen=len(trim(adjustl(DBASE_FieldName)))
      DBASE_FieldTyp = 'C'
      DBASE_FieldLen=len(DBASE_TableRecData13)
      call writeShapFileFieldDesArr(ov_dbasID,fldlen,DBASE_FieldName,&
                                    DBASE_FieldTyp,DBASE_FieldLen)
      DBASE_FieldName ='INDEX'    ! index of contour level
      fldlen=len(trim(adjustl(DBASE_FieldName)))
      DBASE_FieldTyp = 'N'
      DBASE_FieldLen=len(DBASE_TableRecData14)
      call writeShapFileFieldDesArr(ov_dbasID,fldlen,DBASE_FieldName,&
                                    DBASE_FieldTyp,DBASE_FieldLen)
      DBASE_FieldName ='TIME'
      fldlen=len(trim(adjustl(DBASE_FieldName)))
      DBASE_FieldTyp = 'C'
      DBASE_FieldLen=len(DBASE_TableRecData15)
      call writeShapFileFieldDesArr(ov_dbasID,fldlen,DBASE_FieldName,&
                                    DBASE_FieldTyp,DBASE_FieldLen)


      DBASE_FieldTyp = 'N'

      write(ov_dbasID)DBASE_FieldDesTerm

      do irec=1,nrec
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !  FIELDS
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(DBASE_TableRecData12,'(2x,f22.15)')ContourLev(irec)
        write(DBASE_TableRecData14,'(i10)')irec-1

        ! Start of records
        write(ov_dbasID)DBASE_RecStart
        write(ov_dbasID)adjustl(DBASE_TableRecData01)
        write(ov_dbasID)DBASE_TableRecData02
        write(ov_dbasID)DBASE_TableRecData03
        write(ov_dbasID)DBASE_TableRecData04
        write(ov_dbasID)DBASE_TableRecData05
        write(ov_dbasID)DBASE_TableRecData06
        write(ov_dbasID)DBASE_TableRecData07
        write(ov_dbasID)DBASE_TableRecData08
        write(ov_dbasID)DBASE_TableRecData09
        write(ov_dbasID)DBASE_TableRecData10
        write(ov_dbasID)DBASE_TableRecData11
        write(ov_dbasID)DBASE_TableRecData12
        write(ov_dbasID)DBASE_TableRecData13
        write(ov_dbasID)DBASE_TableRecData14
        write(ov_dbasID)DBASE_TableRecData15
      enddo
      write(ov_dbasID)DBASE_EOF

      close(ov_dbasID)

      open(ov_projID, file=trim(adjustl(ov_projfile)), access='stream', form='unformatted', status='replace')
      write(ov_projID)'GEOGCS["GCS_WGS_1984",'
      write(ov_projID)'DATUM["D_WGS_1984",'
      write(ov_projID)'SPHEROID["WGS_1984",'
      write(ov_projID)'6378137,298.257223563]],'
      write(ov_projID)'PRIMEM["Greenwich",0],'
      write(ov_projID)'UNIT["Degree",0.017453292519943295]]'
      close(ov_projID)

      ! Now zip it up if we can
      ! Test if zip is installed
      if(IsLinux.or.IsMacOS)then
        inquire( file=trim(adjustl('/usr/bin/zip')), exist=IsThere)
      elseif(IsWindows)then
        ! this is a placeholder for now
        IsThere = .false.
      else
        IsThere = .false.
      endif
      if(IsThere)then
        write(*,*)ov_zipfile
        write(*,*)ov_projfile
        write(*,*)ov_indxfile
        write(*,*)ov_mainfile
        write(*,*)ov_dbasfile
        write(zipcom,212)'zip -r',ov_zipfile,ov_projfile,&
                                  ov_indxfile,ov_mainfile,ov_dbasfile
 212    format(a6,1x,a12,1x,a12,1x,a12,1x,a12,1x,a12)
        call execute_command_line(zipcom)
      endif

      end subroutine write_ShapeFile_Polyline

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Subroutine writeShapFileFieldDesArr writes the Field Description to the
!  Field Descriptor Array
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine writeShapFileFieldDesArr(ov_dbasID,fldlen,DBASE_FieldName,&
                                   DBASE_FieldTyp,DBASE_FieldLen)

      implicit none

      integer,intent(in)               :: ov_dbasID
      integer,intent(in)               :: fldlen
      character(len=fldlen),intent(in) :: DBASE_FieldName
      character(len=1),intent(in)      :: DBASE_FieldTyp
      integer(kind=1),intent(in)       :: DBASE_FieldLen

      integer(kind=1)  :: DBASE_FieldDecCount
      integer(kind=1)  :: DBASE_FieldSetFieldFlag
      integer(kind=1)  :: DBASE_FieldWrkArID
      integer(kind=1)  :: DBASE_zero     = 0
      integer          :: i

      !  Now the fields followed by the field terminator
      !32-n : 32 bytes: Field descriptor array
      !  0-10  : 11 bytes : Field name in ASCII (zero-filled). 
      !DBASE_FieldName ='name'
      !fldlen=len(trim(adjustl(DBASE_FieldName)))
      write(ov_dbasID)trim(adjustl(DBASE_FieldName))
      do i=fldlen+1,11
        write(ov_dbasID)DBASE_zero
      enddo
      !  11    : 1 byte   : Field type in ASCII (C, D, L, M, or N).
      !DBASE_FieldTyp = 'C'
      write(ov_dbasID)DBASE_FieldTyp
      ! 12-15 : 4 bytes  : Field data address (address is set in memory; not useful on disk). 
      do i=1,4
        write(ov_dbasID)DBASE_zero
      enddo
      ! 16    : 1 byte   : Field length in binary.
      !DBASE_FieldLen=len(DBASE_TableRecData01)   ! Should be 80
      write(ov_dbasID)DBASE_FieldLen
      ! 17    : 1 byte   : Field decimal count in binary. 
      DBASE_FieldDecCount = 0
      write(ov_dbasID)DBASE_FieldDecCount
      ! 18-19 : 2 bytes  : Reserved for dBASE III PLUS on a LAN. 
      do i=1,2
        write(ov_dbasID)DBASE_zero
      enddo
      ! 20    : 1 byte   : Work area ID.
      DBASE_FieldWrkArID = 0
      write(ov_dbasID)DBASE_FieldWrkArID
      ! 21-22 : 2 bytes  : Reserved for dBASE III PLUS on a LAN. 
      do i=1,2
        write(ov_dbasID)DBASE_zero
      enddo
      ! 23    : 1 byte   : SET FIELDS flag.
      DBASE_FieldSetFieldFlag = 0
      write(ov_dbasID)DBASE_FieldSetFieldFlag
      ! 24-32 : 8 bytes  : Reserved bytes
      do i=1,8
        write(ov_dbasID)DBASE_zero
      enddo

      end subroutine writeShapFileFieldDesArr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Subroutine CHECK_ENDIAN checks if the local system uses Big-Endian
!  or Little-Endian byte ordering.  Returns the logical value
!  IsLitEnd = .true. if the system is Little-Endian, .false. otherwise.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      subroutine check_endian(IsLitEnd)
!
!      implicit none
!
!      logical          :: IsLitEnd
!      integer(kind=1)  :: ii(2)
!      integer(kind=2)  :: s
!
!      equivalence (s,ii)
!
!      s = 1
!      IF (ii(1).eq.1)THEN
!        !write(*,*)'System is Little-Endian'
!        IsLitEnd = .true.
!      ELSEIF(ii(2).eq.1)THEN
!        !write(*,*)'System is Big-Endian'
!        IsLitEnd = .false.
!      ELSE
!        write(*,*)'ERROR: cannot figure out endian-ness!'
!        stop
!      ENDIF
!
!      end subroutine check_endian

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Subroutine SWAP_8REAL swaps the byte order of a 8-byte real
!  Subroutine SWAP_4REAL swaps the byte order of a 4-byte real
!  Subroutine SWAP_2INT  swaps the byte order of a 2-byte integer
!  Subroutine SWAP_4INT  swaps the byte order of a 4-byte integer
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      subroutine swap_8real(r)
!
!      implicit none
!
!      integer(kind=1) :: ii(8), jj(8)
!      real(kind=8)    :: r, s, t
!
!      equivalence (s,ii)
!      equivalence (t,jj)
!
!      s = r
!      jj(1) = ii(8)
!      jj(2) = ii(7)
!      jj(3) = ii(6)
!      jj(4) = ii(5)
!      jj(5) = ii(4)
!      jj(6) = ii(3)
!      jj(7) = ii(2)
!      jj(8) = ii(1)
!      r = t
!      end subroutine swap_8real
!
!
!      subroutine swap_4real(r)
!
!      implicit none
!
!      integer(kind=1) :: ii(4), jj(4)
!      real(kind=4)    :: r, s, t
!
!      equivalence (s,ii)
!      equivalence (t,jj)
!
!      s = r
!      jj(1) = ii(4)
!      jj(2) = ii(3)
!      jj(3) = ii(2)
!      jj(4) = ii(1)
!      r = t
!      end subroutine swap_4real
!
!      subroutine swap_2int(r)
!
!      implicit none
!
!      integer(kind=1) :: ii(2), jj(2)
!      integer(kind=2) :: r, s, t
!
!      equivalence (s,ii)
!      equivalence (t,jj)
!
!      s = r
!      jj(1) = ii(2)
!      jj(2) = ii(1)
!      r = t
!      end subroutine swap_2int
!
!      subroutine swap_4int(r)
!
!      implicit none
!
!      integer(kind=1) :: ii(4), jj(4)
!      integer(kind=4) :: r, s, t
!
!      equivalence (s,ii)
!      equivalence (t,jj)
!      s = r
!      jj(1) = ii(4)
!      jj(2) = ii(3)
!      jj(3) = ii(2)
!      jj(4) = ii(1)
!      r = t
!      end subroutine swap_4int

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function BigEnd_2int(isLit,r)

      implicit none

      integer(kind=2) :: BigEnd_2int
      logical         :: isLit
      integer(kind=2) :: r

      integer(kind=2) :: s  = 0
      integer(kind=2) :: t  = 0
      integer         :: bl = 8  ! bit length of the move

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

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      function BigEnd_4int(isLit,r)

      implicit none

      integer(kind=4) :: BigEnd_4int
      logical         :: isLit
      integer(kind=4) :: r

      integer(kind=4) :: s  = 0
      integer(kind=4) :: t  = 0
      integer         :: bl = 8  ! bit length of the move

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

      end function BigEnd_4int

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function LitEnd_2int(isLit,r)

      implicit none

      integer(kind=2) :: LitEnd_2int
      logical         :: isLit
      integer(kind=2) :: r

      integer(kind=2) :: s  = 0
      integer(kind=2) :: t  = 0
      integer         :: bl = 8  ! bit length of the move

!      integer(kind=1) :: ii(2), jj(2)

!      equivalence (s,ii)
!      equivalence (t,jj)

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

      end function LitEnd_2int

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function LitEnd_4int(isLit,r)

      implicit none

      integer(kind=4) :: LitEnd_4int
      logical         :: isLit
      integer(kind=4) :: r

      integer(kind=4) :: s  = 0
      integer(kind=4) :: t  = 0
      integer         :: bl = 8  ! bit length of the move

!      integer(kind=1) :: ii(4), jj(4)

!      equivalence (s,ii)
!      equivalence (t,jj)

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

      end function LitEnd_4int

!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      function BigEnd_8real(isLit,r)
!
!      implicit none
!
!      real(kind=8)    :: BigEnd_8real
!      logical         :: isLit
!      real(kind=8)    :: r
!
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
!
!      end function BigEnd_8real
!
!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function LitEnd_8real(isLit,r)

      implicit none

      real(kind=8)    :: LitEnd_8real
      logical         :: isLit
      real(kind=8)    :: r

      integer(kind=8) :: s  = 0
      integer(kind=8) :: t  = 0
      integer         :: bl = 8  ! bit length of the move

!      integer(kind=1) :: ii(8), jj(8)
!
!      equivalence (s,ii)
!      equivalence (t,jj)

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

      end function LitEnd_8real

