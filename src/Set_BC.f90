      subroutine Set_BC

      ! Set the Boundary conditions

      use precis_param

      use mesh,          only : &
         nxmax,nymax,nzmax,nsmax,ts0,IsPeriodic

      use solution,      only : &
         vx_pd,vy_pd,vz_pd,vf_pd,concen_pd

      implicit none

      integer :: i,j,k,n

      if(IsPeriodic)then
        vx_pd(-1     ,:,:) = vx_pd(nxmax-1,:    ,:)
        vx_pd(0      ,:,:) = vx_pd(nxmax,:    ,:)
        vx_pd(nxmax+1,:,:) = vx_pd(1    ,:    ,:)
        vx_pd(nxmax+2,:,:) = vx_pd(2    ,:    ,:)
      else
        vx_pd(-1     ,:      ,:) = vx_pd(1    ,:    ,:)
        vx_pd(0      ,:      ,:) = vx_pd(1    ,:    ,:)
        vx_pd(nxmax+1,:      ,:) = vx_pd(nxmax,:    ,:)
        vx_pd(nxmax+2,:      ,:) = vx_pd(nxmax,:    ,:)
      endif
      vy_pd(:      ,-1     ,:) = vy_pd(:    ,1    ,:)
      vy_pd(:      ,0      ,:) = vy_pd(:    ,1    ,:)
      vy_pd(:      ,nymax+1,:) = vy_pd(:    ,nymax,:)
      vy_pd(:      ,nymax+2,:) = vy_pd(:    ,nymax,:)
      vz_pd(:,:,     -1)       = vz_pd(:    ,:    ,1) ! A zero vel is not necessary for this slice
      vz_pd(:,:,      0)       = vz_pd(:    ,:    ,1) ! A zero vel is not necessary for this slice
      vz_pd(:,:,nzmax+1)       = vz_pd(:    ,:,nzmax)       
      vz_pd(:,:,nzmax+2)       = vz_pd(:    ,:,nzmax)

      do n=1,nsmax
    
      !***  Top (Z)
        do j=1,nymax
          do i=1,nxmax
            if((vz_pd(i,j,1)+vf_pd(i,j,1,n)).gt.0.0_ip) then    ! Bottom
              concen_pd(i,j,-1,n,ts0) = 0.0_ip
              concen_pd(i,j, 0,n,ts0) = 0.0_ip
            else
              concen_pd(i,j,-1,n,ts0) = concen_pd(i,j,1,n,ts0)
              concen_pd(i,j, 0,n,ts0) = concen_pd(i,j,1,n,ts0)
            endif
            if((vz_pd(i,j,nzmax)+vf_pd(i,j,nzmax,n)).lt.0.0_ip) then    ! Top
              concen_pd(i,j,nzmax+1,n,ts0) = 0.0_ip
              concen_pd(i,j,nzmax+2,n,ts0) = 0.0_ip
            else
              concen_pd(i,j,nzmax+1,n,ts0) = concen_pd(i,j,nzmax,n,ts0)
              concen_pd(i,j,nzmax+2,n,ts0) = concen_pd(i,j,nzmax,n,ts0)
            endif
          enddo
        enddo
   
      !***  Left/Right (X)
        if(IsPeriodic)then
          concen_pd(-1     ,:,:,n,ts0) = concen_pd(nxmax-1,:,:,n,ts0)
          concen_pd(0      ,:,:,n,ts0) = concen_pd(nxmax  ,:,:,n,ts0)
          concen_pd(nxmax+1,:,:,n,ts0) = concen_pd(1      ,:,:,n,ts0)
          concen_pd(nxmax+2,:,:,n,ts0) = concen_pd(2      ,:,:,n,ts0)
        else
          do k=1,nzmax
            do j=1,nymax
              ! Left side
              ! (0 if it's blowing in, dc/dx=0 if it's blowing out)
              if(vx_pd(1,j,k).gt.0.0_ip) then
                concen_pd(-1,j,k,n,ts0)    = 0.0_ip
                concen_pd( 0,j,k,n,ts0)    = 0.0_ip
              else
                concen_pd(-1,j,k,n,ts0)    = concen_pd(1,j,k,n,ts0)
                concen_pd( 0,j,k,n,ts0)    = concen_pd(1,j,k,n,ts0)
              endif
              ! Right side
              ! (0 if it's blowing in, dc/dx=0 if it's blowing out)
              if(vx_pd(nxmax,j,k).lt.0.0_ip) then
                concen_pd(nxmax+1,j,k,n,ts0) = 0.0_ip
                concen_pd(nxmax+2,j,k,n,ts0) = 0.0_ip
              else
                concen_pd(nxmax+1,j,k,n,ts0) = concen_pd(nxmax,j,k,n,ts0)
                concen_pd(nxmax+2,j,k,n,ts0) = concen_pd(nxmax,j,k,n,ts0)
              endif
            enddo
          enddo
        endif

      !***  Up/Down (Y)
        do k=1,nzmax
          do i=1,nxmax
             ! Up
             ! (0 if it's blowing in, dc/dy=0 if it's blowing out)
            if(vy_pd(i,1,k).gt.0.0_ip) then
              concen_pd(i,-1,k,n,ts0)    = 0.0_ip
              concen_pd(i, 0,k,n,ts0)    = 0.0_ip
            else
              concen_pd(i,-1,k,n,ts0)    = concen_pd(i,1,k,n,ts0)
              concen_pd(i, 0,k,n,ts0)    = concen_pd(i,1,k,n,ts0)
            endif
             ! Down
             ! (0 if it's blowing in, dc/dy=0 if it's blowing out)
            if(vy_pd(i,nymax,k).lt.0.0_ip) then
              concen_pd(i,nymax+1,k,n,ts0) = 0.0_ip
              concen_pd(i,nymax+2,k,n,ts0) = 0.0_ip
            else
              concen_pd(i,nymax+1,k,n,ts0) = concen_pd(i,nymax,k,n,ts0)
              concen_pd(i,nymax+2,k,n,ts0) = concen_pd(i,nymax,k,n,ts0)
            endif
          enddo
        enddo
    
        !***  Corners
        ! Each ghost of the eight corners are set to the value of the corner of the main grid
        concen_pd(     -1:0      ,     -1:0      ,     -1:0      ,n,ts0)=concen_pd(1    ,1    ,1    ,n,ts0)
        concen_pd(nxmax+1:nxmax+2,     -1:0      ,     -1:0      ,n,ts0)=concen_pd(nxmax,1    ,1    ,n,ts0)
        concen_pd(     -1:0      ,nymax+1:nymax+2,     -1:0      ,n,ts0)=concen_pd(1    ,nymax,1    ,n,ts0)
        concen_pd(nxmax+1:nxmax+2,nymax+1:nymax+2,     -1:0      ,n,ts0)=concen_pd(nxmax,nymax,1    ,n,ts0)
        concen_pd(     -1:0      ,     -1:0      ,nzmax+1:nzmax+2,n,ts0)=concen_pd(1    ,1    ,nzmax,n,ts0)
        concen_pd(nxmax+1:nxmax+2,     -1:0      ,nzmax+1:nzmax+2,n,ts0)=concen_pd(nxmax,1    ,nzmax,n,ts0)
        concen_pd(     -1:0      ,nymax+1:nymax+2,nzmax+1:nzmax+2,n,ts0)=concen_pd(1    ,nymax,nzmax,n,ts0)
        concen_pd(nxmax+1:nxmax+2,nymax+1:nymax+2,nzmax+1:nzmax+2,n,ts0)=concen_pd(nxmax,nymax,nzmax,n,ts0)

      enddo

      end subroutine Set_BC
