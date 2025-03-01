!##############################################################################
!
! AdvectionHorz
!
! This module manages the advection routines used for horizontal advection.
! Only one subroutine is in this module, AdvectHorz, which controls the
! horizontal advection using the 1-d donor-cell-upwind (DCU) method.  This is
! where alternate advection schemes could be invoked, such and the 2-d
! corner-transport-upwind (CTU) or semi-Lagrange.
!
!##############################################################################
 
      module AdvectionHorz

      use precis_param

      implicit none

        ! Set everything to private by default
      private

        ! Publicly available subroutines/functions
      public AdvectHorz

      contains

      !------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  AdvectHorz(itoggle)
!
!  Called from: Ash3d.F90
!  Arguments:
!    itoggle : integer time step that is evaluated for even/odd
!
!  This subroutine is called from the main time loop in Ash3d.F90 and invokes
!  the horizontal advection routines.  This routine is called after the
!  boundary conditions are set and before the vertical advection routine and
!  the diffusion routines.  The default scheme is DCU, provided by the
!  module AdvectionHorz_DCU with x and y advection treated separately.  A
!  toggle variable is provided to determine which advection direction to
!  apply first (even-> y then x; odd-> x then y). 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine AdvectHorz(itoggle)

#ifndef FAST_SUBGRID
      use mesh,          only : &
         nxmax,nymax,nzmax

      use solution,      only : &
         imin,imax,jmin,jmax,kmin,kmax
#endif

      use AdvectionHorz_DCU, only : &
         advect_x,advect_y

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
        ! Use dimension splitting with donor-cell-upwind (DCU)
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

      end module AdvectionHorz

!##############################################################################
