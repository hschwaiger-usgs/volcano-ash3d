 !****************************************************************************
 
 !This file contains the module AdvectionHorz and subroutine AdvectHorz.
 !  The subroutine calls other subroutines, advect_x and advect_y, that
 !  calculate horizontal advection in the x and y directions.  
 !  AdvectHorz is called during the time loop in the main program
 !  Ash3d.F90, after setting boundary conditions and before calculating
 !  vertical advection.
 
 !***************************************************************************
 
      module AdvectionHorz

      use precis_param

      implicit none

        ! Set everything to private by default
      private

        ! Publicly available subroutines/functions
      public AdvectHorz

      contains

!****************************************************************************

      subroutine AdvectHorz(itoggle)

      use AdvectionHorz_DCU, only : &
         advect_x,advect_y

      use solution,      only : &
         imin,imax,jmin,jmax,kmin,kmax

      use mesh,          only : &
         nxmax,nymax,nzmax

      integer,intent(in) :: itoggle

#ifdef FAST_SUBGRID
      INTERFACE
        subroutine get_minmax_index
        end subroutine get_minmax_index
      END INTERFACE

      call get_minmax_index
#else
      imin = 1
      imax = nxmax
      jmin = 1
      jmax = nymax
      kmin = 1
      kmax = nzmax

#endif
        ! Use dimension splitting with donor cell upwind (DCU)
      if(mod(itoggle,2).eq.0) then
          ! for even time steps, advect in y first
        call advect_y
        call advect_x
      else
        call advect_x
        call advect_y
      endif

      end subroutine AdvectHorz


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module AdvectionHorz

