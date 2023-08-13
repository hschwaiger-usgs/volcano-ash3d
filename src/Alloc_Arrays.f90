!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  alloc_arrays()
!
!  Called from: Ash3d.F90
!  Arguments:
!    none
!
!  This subroutine calls the allocation routines from each of the modules in
!  Ash3d_VariableModules.f90
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine alloc_arrays

      use precis_param

      use io_units

      use mesh,          only : &
         nxmax,nymax,nzmax,&
           Allocate_mesh

      use solution,      only : &
           Allocate_solution

      use wind_grid,     only : &
           Allocate_wind_grid

      use Diffusion,     only : &
         Allocate_Diff

      use Output_Vars,   only : &
         Allocate_Output_Vars

      use Source,        only : &
         Allocate_Source_grid

      implicit none

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"--------------------------------------------------"
        write(outlog(io),*)"---------- ALLOC_ARRAYS --------------------------"
        write(outlog(io),*)"--------------------------------------------------"
      endif;enddo

      ! io_data is allocated in Read_Control_File

      call Allocate_mesh
      call Allocate_solution
      call Allocate_wind_grid
      call Allocate_Output_Vars
      call Allocate_Source_grid
      call Allocate_Diff

      end subroutine alloc_arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  dealloc_arrays()
!
!  Called from: Ash3d.F90
!  Arguments:
!    none
!
!  This subroutine calls the deallocation routines from each of the modules in
!  Ash3d_VariableModules.f90
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine dealloc_arrays

      use io_units

      use io_data,       only : &
           Deallocate_io_data

      use mesh,          only : &
           Deallocate_mesh

      use solution,      only : &
           Deallocate_solution

      use wind_grid,     only : &
           Deallocate_wind_grid

      use Tephra,        only : &
           Deallocate_Tephra, Deallocate_Tephra_Met

      use Source,        only : &
           Deallocate_Source

      use Diffusion,     only : &
           Deallocate_Diff

      use Output_Vars,   only : &
           Deallocate_Output_Vars, Deallocate_NTime, Deallocate_Profile, Deallocate_Output_UserVars

      use Airports,      only : &
           Deallocate_Airports

      use Atmosphere,    only : &
           Deallocate_Atmosphere_Met

      use MetReader,     only : &
           MR_Reset_Memory

      implicit none

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Cleaning up allocated memory"
      endif;enddo

      call Deallocate_io_data
      call Deallocate_mesh
      call Deallocate_solution
      call Deallocate_wind_grid
      call Deallocate_Tephra
      call Deallocate_Tephra_Met
      call Deallocate_Source
      call Deallocate_Diff
      call Deallocate_Output_Vars
      call Deallocate_NTime
      call Deallocate_Profile
      call Deallocate_Output_UserVars
      call Deallocate_Airports
      call Deallocate_Atmosphere_Met

      call MR_Reset_Memory

      end subroutine dealloc_arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
