
      subroutine write_ShapeFile_Polyline(iprod,itime)

      ! input should be the contours, the output product, and the time step

      use precis_param

      use global_param,  only : &
        IsLitEnd

      use Output_Vars,   only : &
         ContourDataX,ContourDataY,ContourDataNcurves,ContourDataNpoints,&
         Contour_MaxCurves,Contour_MaxPoints,ContourLev,Con_N

      use io_data,       only : &
         nWriteTimes,WriteTimes,VolcanoName

      implicit none

      integer,intent(in) :: iprod
      integer,intent(in) :: itime

      character(len=8)  :: ov_fileroot = "testfile"
      character(len=22) :: ov_mainfile
      character(len=22) :: ov_indxfile
      character(len=22) :: ov_dbasfile
      character(len=22) :: ov_projfile
      character(len=4)  :: ov_mainext = ".shp"
      character(len=4)  :: ov_indxext = ".shx"
      character(len=4)  :: ov_dbasext = ".dbf"
      character(len=4)  :: ov_projext = ".prj"
      integer           :: ov_mainID  = 22
      integer           :: ov_indxID  = 23
      integer           :: ov_dbasID  = 24
      integer           :: ov_projID  = 25

      character(len=40) :: title_plot
      character(len=15) :: title_units

      logical           :: debugmode = .false.

      character(len=60) :: linebuffer60
      integer                                  :: nrec      ! num of records (e.g. contour levels)
      integer(kind=4),dimension(:),allocatable :: NumPoints ! total points
      integer(kind=4) :: NumParts_irec  ! number of parts in level
      integer(kind=4) :: NumPoints_irec ! sum total points of all curves (parts) in level
      real(kind=8),dimension(:,:,:),allocatable :: contx
      real(kind=8),dimension(:,:,:),allocatable :: conty
      real(kind=8),dimension(:,:,:),allocatable :: contz
      integer :: Iostatus
      integer :: i,j,irec,ipart,ipt
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
      integer(kind=4) :: dumint(2)

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
      integer(kind=1)  :: DBASE_FieldDecCount
      integer(kind=1)  :: DBASE_FieldWrkArID
      integer(kind=1)  :: DBASE_FieldDesTerm
      integer(kind=1)  :: DBASE_FieldSetFieldFlag
      integer          :: fldlen
      character(len=1) :: DBASE_RecStart = ' '
      integer(kind=1)  :: DBASE_EOF      = 26
      character(len=80):: DBASE_TableRecDataName
      character(len=24):: DBASE_TableRecDataValue
      character(len=10):: DBASE_TableRecDataIndex

      INTERFACE
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
        ov_fileroot = 'depothik'
        title_units = '(mm)'
      elseif(iprod.eq.4)then   ! deposit at specified times (inches)
        write(title_plot,'(a22,f5.2,a6)')': Deposit Thickness t=',WriteTimes(itime),' hours'
        ov_fileroot = 'depothik'
        title_units  = '(in)'
      elseif(iprod.eq.5)then       ! deposit at final time (mm)
        title_plot = ': Final Deposit Thickness'
        ov_fileroot = 'depothik'
        title_units = '(mm)'
      elseif(iprod.eq.6)then   ! deposit at final time (inches)
        title_plot = ': Final Deposit Thickness'
        ov_fileroot = 'depothik'
        title_units = '(in)'
      elseif(iprod.eq.7)then   ! ashfall arrival time (hours)
        write(title_plot,'(a22)')': Ashfall arrival time'
        ov_fileroot = 'DepAvlTm'
        title_units = '(hours)'
      elseif(iprod.eq.8)then   ! ashfall arrival at airports/POI (mm)
        write(*,*)"ERROR: No map PNG output option for airport arrival time data."
        write(*,*)"       Should not be in write_2Dmap_PNG_dislin"
        stop 1
      elseif(iprod.eq.9)then   ! ash-cloud concentration
        write(title_plot,'(a28,f5.2,a6)')': Ash-cloud concentration t=',WriteTimes(itime),' hours'
        ov_fileroot = 'AshCdCon'
        title_units = '(mg/m3)'
      elseif(iprod.eq.10)then   ! ash-cloud height
        write(title_plot,'(a21,f5.2,a6)')': Ash-cloud height t=',WriteTimes(itime),' hours'
        ov_fileroot = 'AshCdHgt'
        title_units = '(km)'
      elseif(iprod.eq.11)then   ! ash-cloud bottom
        write(title_plot,'(a21,f5.2,a6)')': Ash-cloud bottom t=',WriteTimes(itime),' hours'
        ov_fileroot = 'AshCdBot'
        title_units = '(km)'
      elseif(iprod.eq.12)then   ! ash-cloud load
        write(title_plot,'(a19,f5.2,a6)')': Ash-cloud load t=',WriteTimes(itime),' hours'
        ov_fileroot = 'AshCdLod'
        title_units = '(T/km2)'
      elseif(iprod.eq.13)then  ! radar reflectivity
        write(title_plot,'(a26,f5.2,a6)')': Ash-cloud radar refl. t=',WriteTimes(itime),' hours'
        ov_fileroot = 'AshClRad'
        title_units = '(dBz)'
      elseif(iprod.eq.14)then   ! ashcloud arrival time (hours)
        write(title_plot,'(a24)')': Ash-cloud arrival time'
        ov_fileroot = 'AshAvlTm'
        title_units = '(hours)'
      elseif(iprod.eq.15)then   ! topography
        write(title_plot,'(a12)')': Topography'
        ov_fileroot = 'Topogrph'
        title_units = '(km)'
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
      do i = 1,Con_N
        if(ContourDataNcurves(i).gt.0)nrec=nrec+1
      enddo

      ! We allocate variables here to hold the contour data mainly so ensure that
      ! the variables we write to the shapefile have exactly the kind values expected.
      allocate(NumPoints(nrec))
      allocate(reclen(nrec))
      allocate(xmin(nrec)); xmin(:) =  1.0e8_8
      allocate(xmax(nrec)); xmax(:) = -1.0e8_8
      allocate(ymin(nrec)); ymin(:) =  1.0e8_8
      allocate(ymax(nrec)); ymax(:) = -1.0e8_8

      do irec = 1,nrec  ! This is the loop of the layers (records)
        NumParts_irec = ContourDataNcurves(irec)
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
                       NumParts_irec*4 + &  ! Parts      : NumParts-Integers : Little
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

        NumParts_irec = ContourDataNcurves(irec)
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
        write(ov_mainID)LitEnd_4int(IsLitEnd,NumParts_irec)
        write(ov_mainID)LitEnd_4int(IsLitEnd,NumPoints(irec))

        if(debugmode)then
          write(*,*)shape_type,"Byte 0 Shape Type 3 Integer 1 Little"
          write(*,*)xmin(irec),"Byte 4 xmin Double Little"
          write(*,*)ymin(irec),"Byte 12 ymin Double Little"
          write(*,*)xmax(irec),"Byte 20 xmax Double Little"
          write(*,*)ymax(irec),"Byte 28 ymax Double Little"
          write(*,*)NumParts_irec,"Byte 36 NumParts NumParts Integer 1 Little"
          write(*,*)NumPoints(irec),"Byte 40 NumPoints NumPoints Integer 1 Little"
        endif
        ! Address of the start of the first part is 0
        write(ov_mainID)LitEnd_4int(IsLitEnd,0_4) ! An array of length NumParts with each value
                                                  ! the address (zero-offset) of the start of the part
        do i=1,NumParts_irec-1
          tmp4 = sum(ContourDataNpoints(irec,1:i)) !-1
          write(ov_mainID)LitEnd_4int(IsLitEnd,tmp4)
        enddo

        do ipart=1,ContourDataNcurves(irec)
          do i=1,ContourDataNpoints(irec,ipart)
            !if(2==1)then
            !  write(ov_mainID)LitEnd_8real(IsLitEnd,contx(irec,ipart,i))
            !  write(ov_mainID)LitEnd_8real(IsLitEnd,conty(irec,ipart,i))
            !else
              write(ov_mainID)LitEnd_8real(IsLitEnd,ContourDataX(irec,ipart,i))
              write(ov_mainID)LitEnd_8real(IsLitEnd,ContourDataY(irec,ipart,i))
            !endif
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
      ! Full specifications at
      ! https://www.dbase.com/Knowledgebase/INT/db7_file_fmt.htm
      ! For dbase 5 see https://www.oocities.org/geoff_wass/dBASE/GaryWhite/dBASE/FAQ/qformt.htm
      ! Shapefile additional requirements:
      !  1. name consistent with shp, but with dbf extension
      !  2. table must contain one record per shape feature
      !  3. The record order must be the same as the order of shape features in the main file
      !  4. The year value in the dBASE header must be the year since 1900.
      ! Note: We will write dBASE V – MS-Windows (Level 5) format since we know that works.
      open(ov_dbasID, file=trim(adjustl(ov_dbasfile)), access='stream', form='unformatted', status='replace')

      write(DBASE_TableRecDataName,*)trim(adjustl(VolcanoName)),title_plot,title_units
      DBASE_TableRecDataValue =  "       0.000000000000000 "
      DBASE_TableRecDataIndex =  "         0"

      DBASE_zero = 0
      DBASE_v  = 3
      DBASE_yy =95
      DBASE_mm =07
      DBASE_dd =26

      DBASE_nrec      = nrec
      DBASE_headlen   = 32 + (32 + 32 + 32) + 1
      DBASE_reclen    = (len(DBASE_TableRecDataName) + &
                         len(DBASE_TableRecDataValue)+ &
                         len(DBASE_TableRecDataIndex)) +1
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

      !  Now the fields followed by the field terminator
      !32-n : 32 bytes: Field descriptor array
      !  0-10  : 11 bytes : Field name in ASCII (zero-filled). 
      DBASE_FieldName ='name'
      fldlen=len(trim(adjustl(DBASE_FieldName)))
      write(ov_dbasID)trim(adjustl(DBASE_FieldName))
      do i=fldlen+1,11
        write(ov_dbasID)DBASE_zero
      enddo
      !  11    : 1 byte   : Field type in ASCII (C, D, L, M, or N).
      DBASE_FieldTyp = 'C'
      write(ov_dbasID)DBASE_FieldTyp
      ! 12-15 : 4 bytes  : Field data address (address is set in memory; not useful on disk). 
      do i=1,4
        write(ov_dbasID)DBASE_zero
      enddo
      ! 16    : 1 byte   : Field length in binary.
      DBASE_FieldLen=len(DBASE_TableRecDataName)   ! Should be 80
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
  
      ! Next field for the value of the first record
      !  0-10  : 11 bytes : Field name in ASCII (zero-filled). 
      DBASE_FieldName ='value'
      fldlen=len(trim(adjustl(DBASE_FieldName)))
      write(ov_dbasID)trim(adjustl(DBASE_FieldName))
      do i=fldlen+1,11
        write(ov_dbasID)DBASE_zero
      enddo
      !  11    : 1 byte   : Field type in ASCII (C, D, L, M, or N).
      DBASE_FieldTyp = 'N'
      write(ov_dbasID)DBASE_FieldTyp
      ! 12-15 : 4 bytes  : Field data address (address is set in memory; not useful on disk). 
      do i=1,4
        write(ov_dbasID)DBASE_zero
      enddo
      ! 16    : 1 byte   : Field length in binary.
      DBASE_FieldLen=len(DBASE_TableRecDataValue)  ! Should be 24
      write(ov_dbasID)DBASE_FieldLen
      ! 17    : 1 byte   : Field decimal count in binary. 
      DBASE_FieldDecCount = 15  ! number of digits past decimal??
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
  
      ! Next field for the index of the first record
      !  0-10  : 11 bytes : Field name in ASCII (zero-filled). 
      DBASE_FieldName  ='index'
      fldlen=len(trim(adjustl(DBASE_FieldName)))
      write(ov_dbasID)trim(adjustl(DBASE_FieldName))
      do i=fldlen+1,11
        write(ov_dbasID)DBASE_zero
      enddo
      !  11    : 1 byte   : Field type in ASCII (C, D, L, M, or N).
      DBASE_FieldTyp = 'N'
      write(ov_dbasID)DBASE_FieldTyp
      ! 12-15 : 4 bytes  : Field data address (address is set in memory; not useful on disk). 
      do i=1,4
        write(ov_dbasID)DBASE_zero
      enddo
      ! 16    : 1 byte   : Field length in binary.
      DBASE_FieldLen=len(DBASE_TableRecDataIndex)  ! Should be 10
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

      write(ov_dbasID)DBASE_FieldDesTerm

      do irec=1,nrec
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !  FIELDS
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(DBASE_TableRecDataValue,'(2x,f22.15)')ContourLev(irec)
        write(DBASE_TableRecDataIndex,'(i10)')irec-1

        ! Start of records
        write(ov_dbasID)DBASE_RecStart
        write(ov_dbasID)adjustl(DBASE_TableRecDataName)
        write(ov_dbasID)DBASE_TableRecDataValue
        write(ov_dbasID)DBASE_TableRecDataIndex
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

      end subroutine write_ShapeFile_Polyline

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

      integer(kind=2) :: ii(2), jj(2)
      integer(kind=2) :: s, t

      equivalence (s,ii)
      equivalence (t,jj)

      if(isLit)then
        ! We need r to be big-endian, but the system is little-endian
        ! swap the bytes
        s = r
        jj(1) = ii(2)
        jj(2) = ii(1)
        BigEnd_2int = t
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

      integer(kind=1) :: ii(4), jj(4)
      integer(kind=4) :: s, t

      equivalence (s,ii)
      equivalence (t,jj)

      if(isLit)then
        ! We need r to be big-endian, but the system is little-endian
        ! swap the bytes
        s = r
        jj(1) = ii(4)
        jj(2) = ii(3)
        jj(3) = ii(2)
        jj(4) = ii(1)
        BigEnd_4int = t
      else
        BigEnd_4int = r
      endif

      end function BigEnd_4int

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function LitEnd_2int(isLit,r)

      implicit none

      integer(kind=2) :: LitEnd_2int
      logical         :: isLit
      integer(kind=2) :: r

      integer(kind=1) :: ii(2), jj(2)
      integer(kind=2) :: s, t

      equivalence (s,ii)
      equivalence (t,jj)

      if(isLit)then
        LitEnd_2int = r
      else
        ! We need r to be little-endian, but the system is big-endian
        ! swap the bytes
        s = r
        jj(1) = ii(2)
        jj(2) = ii(1)
        LitEnd_2int = t
      endif

      end function LitEnd_2int

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function LitEnd_4int(isLit,r)

      implicit none

      integer(kind=4) :: LitEnd_4int
      logical         :: isLit
      integer(kind=4) :: r

      integer(kind=1) :: ii(4), jj(4)
      integer(kind=4) :: s, t

      equivalence (s,ii)
      equivalence (t,jj)

      if(isLit)then
        LitEnd_4int = r
      else
        ! We need r to be little-endian, but the system is big-endian
        ! swap the bytes
        s = r
        jj(1) = ii(4)
        jj(2) = ii(3)
        jj(3) = ii(2)
        jj(4) = ii(1)
        LitEnd_4int = t
      endif

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

      integer(kind=1) :: ii(8), jj(8)
      real(kind=8)    :: s, t

      equivalence (s,ii)
      equivalence (t,jj)

      if(isLit)then
        LitEnd_8real = r
      else
        ! We need r to be little-endian, but the system is big-endian
        ! swap the bytes
        s = r
        jj(1) = ii(8)
        jj(2) = ii(7)
        jj(3) = ii(6)
        jj(4) = ii(5)
        jj(5) = ii(4)
        jj(6) = ii(3)
        jj(7) = ii(2)
        jj(8) = ii(1)
        LitEnd_8real = t
      endif

      end function LitEnd_8real

