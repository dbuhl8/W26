subroutine soln_update(dt)

#include "definition.h"
  
  use grid_data
  use primconsflux, only : cons2prim,prim2cons

  implicit none
  real, intent(IN) :: dt
  integer :: i, pidx1, pidx2
  real :: dtx

  dtx = dt/gr_dx

  !! update conservative vars
  do i = gr_ibeg, gr_iend+1
    ! convert prim to vars
    call prim2cons(gr_V(:,i),gr_U(:,i))
  
    ! forward euler 
    gr_U(:,i) = gr_U(:,i) - dtx*(gr_flux(:,i) - gr_flux(:,i-1))
  end do

  !pidx1 = gr_ibeg + gr_nx/2 - 1
  !pidx2 = gr_ibeg + gr_nx/2
  !print *, '   x          Frho         Fm          FE   '
  !print *, '________    ________    ________    ________'
  !print '(4(F8.3, "    "))', gr_xCoord(pidx1), gr_flux(:,pidx1)
  !print '(4(F8.3, "    "))', gr_xCoord(pidx2), gr_flux(:,pidx2)

  do i = gr_ibeg, gr_iend
    ! Eos is automatically called inside cons2prim
    call cons2prim(gr_U(:,i),gr_V(:,i))
  end do


  return
end subroutine soln_update
