subroutine soln_update(dt)

#include "definition.h"
  
  use grid_data
  use primconsflux, only : cons2prim,prim2cons

  implicit none
  real, intent(IN) :: dt
  integer :: i
  real :: dtx

  dtx = dt/gr_dx

  !! update conservative vars
  do i = gr_ibeg, gr_iend
     !! first convert prim vars to cons vars
     call prim2cons(gr_V(:,i),gr_U(:,i))
     
     !! let's update
     gr_U(:,i) = gr_U(:,i) - &
          dtx*(gr_flux(:,i+1) - gr_flux(:,i))
  end do


  !! get updated primitive vars from the updated conservative vars
  do i = gr_ibeg, gr_iend
     ! Eos is automatically called inside cons2prim
     call cons2prim(gr_U(:,i),gr_V(:,i))
  end do
  

  return
end subroutine soln_update
