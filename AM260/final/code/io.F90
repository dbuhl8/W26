module io

#include "definition.h"
  
  use grid_data, only : gr_xCoord, gr_V, gr_ibeg, gr_iend, gr_xbeg, gr_xend, gr_nx
  use sim_data, only : sim_name
  implicit none

  
  integer, save :: nCounter
  
contains

  subroutine io_writeGridConfig(ioCounter)
    implicit none

    integer, intent(IN) :: ioCounter
    character(len=35) :: ofile

    ! file name for ascii output
    ofile = 'grid_'//trim(sim_name)//'.dat'

    open(unit=20,file=ofile,status='unknown')
    write(20,'(2i5,2F16.8)') ioCounter, gr_nx, gr_xbeg, gr_xend
    close(20)
  end subroutine io_writeGridConfig

  subroutine io_writeOutput(t,nstep,ioCounter)
    implicit none

    real, intent(IN) :: t
    integer, intent(IN) :: nstep,ioCounter
    
    integer :: i,nVar,nCell
    character(len=35) :: ofile
    character(len=5)  :: cCounter


    ! convert conter number to character
    write(cCounter,910) ioCounter + 10000

    ! file name for ascii output
    ofile = 'slug_'//trim(sim_name)//'_'//cCounter//'.dat'

    open(unit=20,file=ofile,status='unknown')
    do i=gr_ibeg,gr_iend
       write(20,920) t,gr_xCoord(i),(gr_V(nVar,i),nVar=1,NUMB_VAR)
    end do

    
910 format(i5)
920 format(1x,f16.9,f16.8,1x,NUMB_VAR f32.16)
    
    close(20)
  end subroutine io_writeOutput
  
end module io
