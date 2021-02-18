      subroutine Set_BC

      ! Set the Boundary conditions

      use precis_param

      use mesh,          only : &
         nxmax,nymax,nzmax,ts0,IsPeriodic

      use solution,      only : &
         vx_pd,vy_pd,vz_pd,vf_pd,concen_pd

      implicit none


      !------------------------------------------------------------------------
      !   VELOCITIES
      !------------------------------------------------------------------------
      ! Use zero-order extrapolation for all velocities (i.e. constant
      ! values extrapolation from edge cell to ghost cells).  This is
      ! important to supress any effects from the \Delta u terms of
      ! advection.  A linear extrapolation might seem better, but this
      ! could allow a non-physical inflow.
      if(IsPeriodic)then
        vx_pd(-1     ,:,:) = vx_pd(nxmax-1,:    ,:)
        vx_pd(0      ,:,:) = vx_pd(nxmax,:    ,:)
        vx_pd(nxmax+1,:,:) = vx_pd(1    ,:    ,:)
        vx_pd(nxmax+2,:,:) = vx_pd(2    ,:    ,:)
      else
        vx_pd(-1     ,:      ,:) = vx_pd(1    ,:    ,:)
        vx_pd( 0     ,:      ,:) = vx_pd(1    ,:    ,:)
        vx_pd(nxmax+1,:      ,:) = vx_pd(nxmax,:    ,:)
        vx_pd(nxmax+2,:      ,:) = vx_pd(nxmax,:    ,:)
      endif
      vy_pd(:      ,-1     ,:) = vy_pd(:    ,1    ,:)
      vy_pd(:      , 0     ,:) = vy_pd(:    ,1    ,:)
      vy_pd(:      ,nymax+1,:) = vy_pd(:    ,nymax,:)
      vy_pd(:      ,nymax+2,:) = vy_pd(:    ,nymax,:)

      vz_pd(:,:,     -1)       = vz_pd(:    ,:    ,1)
      vz_pd(:,:,      0)       = vz_pd(:    ,:    ,1)
      vz_pd(:,:,nzmax+1)       = vz_pd(:    ,:,nzmax)       
      vz_pd(:,:,nzmax+2)       = vz_pd(:    ,:,nzmax)

      ! Can't forget to apply the same extrapolation to fall
      ! velocities
      vf_pd(:,:,     -1,:)       = vf_pd(:    ,:    ,1,:)
      vf_pd(:,:,      0,:)       = vf_pd(:    ,:    ,1,:)
      vf_pd(:,:,nzmax+1,:)       = vf_pd(:    ,:,nzmax,:)
      vf_pd(:,:,nzmax+2,:)       = vf_pd(:    ,:,nzmax,:)

      !------------------------------------------------------------------------
      !   CONCENTRATIONS
      !------------------------------------------------------------------------
      ! Zero the concentration on all ghost cells.  
      !***  Left/Right (X)
      concen_pd(     -1:0      ,:,:,:,ts0) = 0.0_ip
      concen_pd(nxmax+1:nxmax+2,:,:,:,ts0) = 0.0_ip
      !***  Up/Down (Y)
      concen_pd(:,     -1:0      ,:,:,ts0) = 0.0_ip
      concen_pd(:,nymax+1:nymax+2,:,:,ts0) = 0.0_ip
      !***  Bottom/Top (Z)
      concen_pd(:,:,     -1:0      ,:,ts0) = 0.0_ip
      concen_pd(:,:,nzmax+1:nzmax+2,:,ts0) = 0.0_ip

      if(IsPeriodic)then
        concen_pd(-1     ,:,:,:,ts0) = concen_pd(nxmax-1,:,:,:,ts0)
        concen_pd(0      ,:,:,:,ts0) = concen_pd(nxmax  ,:,:,:,ts0)
        concen_pd(nxmax+1,:,:,:,ts0) = concen_pd(1      ,:,:,:,ts0)
        concen_pd(nxmax+2,:,:,:,ts0) = concen_pd(2      ,:,:,:,ts0)
      endif

      end subroutine Set_BC
