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

      write(global_production,*)"--------------------------------------------------"
      write(global_production,*)"---------- ALLOC_ARRAYS --------------------------"
      write(global_production,*)"--------------------------------------------------"

      ! io_data is allocated in Read_Control_File

      call Allocate_mesh
      call Allocate_solution
      call Allocate_wind_grid
      call Allocate_Output_Vars(nxmax,nymax,nzmax)
      call Allocate_Source_grid(nxmax,nymax,nzmax)
      call Allocate_Diff(nxmax,nymax,nzmax)

      end subroutine alloc_arrays
!###############################################################################


!###############################################################################
      subroutine dealloc_arrays

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
           Deallocate_Output_Vars, Deallocate_Output_UserVars

      use Airports,      only : &
           Deallocate_Airports

      use MetReader,     only : &
           MR_Reset_Memory

      implicit none

      call Deallocate_io_data
      call Deallocate_mesh
      call Deallocate_solution
      call Deallocate_wind_grid
      call Deallocate_Tephra
      call Deallocate_Tephra_Met
      call Deallocate_Source
      call Deallocate_Diff
      call Deallocate_Output_Vars
      call Deallocate_Output_UserVars
      call Deallocate_Airports

      call MR_Reset_Memory

      end subroutine dealloc_arrays

!##############################################################################

