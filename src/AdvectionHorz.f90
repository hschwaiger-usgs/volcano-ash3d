      module AdvectionHorz

      use precis_param

      contains

!******************************************************************************

      subroutine AdvectHorz(itoggle)

      use AdvectionHorz_DCU, only : &
         advect_x,advect_y

      implicit none

      integer,intent(in) :: itoggle

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

