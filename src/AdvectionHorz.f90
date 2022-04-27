      module AdvectionHorz

      use precis_param

      contains

!******************************************************************************

      subroutine AdvectHorz(itoggle)

      use AdvectionHorz_DCU, only : &
         advect_x,advect_y

      use solution,      only : &
         imin,imax,jmin,jmax,kmin,kmax

      use mesh,          only : &
         nxmax,nymax,nzmax

      implicit none

      integer,intent(in) :: itoggle


#ifdef FAST_SUBGRID
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

