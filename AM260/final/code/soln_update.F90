subroutine soln_update(dt)

#include "definition.h"
  
  use grid_data
  use primconsflux, only : cons2prim,prim2cons

  implicit none
  real, intent(IN) :: dt
  integer :: i, pidx1, pidx2, pidx3, pidx4
  real :: dtx

  dtx = dt/gr_dx
  pidx1 = gr_ibeg + gr_nx/4 + 1
  pidx2 = gr_ibeg + gr_nx/4 + 2
  pidx3 = gr_ibeg + gr_nx/4 - 1
  pidx4 = gr_ibeg + gr_nx/4 - 2
  !print *, 'dtx: ', dtx
  !print *, '   x           rho          m           E   '
  !print *, '________    ________    ________    ________'
  !print '(4(F8.3, "    "))', gr_xCoord(pidx1), gr_U(:,pidx1)
  !print '(4(F8.3, "    "))', gr_xCoord(pidx2), gr_flux(:,pidx2)
  !print '(4(F8.3, "    "))', gr_xCoord(pidx3), gr_flux(:,pidx3)
  !print '(4(F8.3, "    "))', gr_xCoord(pidx4), gr_flux(:,pidx4)
  !print *, '   x          Frho        Fm           FE   '
  !print *, '________    ________    ________    ________'
  !print '(4(F8.3, "    "))', gr_xCoord(pidx1), gr_flux(:,pidx1)
  !print '(4(F8.3, "    "))', gr_xCoord(pidx2), gr_flux(:,pidx2)
  !print '(4(F8.3, "    "))', gr_xCoord(pidx3), gr_flux(:,pidx3)
  !print '(4(F8.3, "    "))', gr_xCoord(pidx4), gr_flux(:,pidx4)

  !! update conservative vars
  do i = gr_ibeg, gr_iend+1
    ! convert prim to vars
    call prim2cons(gr_V(:,i),gr_U(:,i))

    !if (abs(gr_ibeg + gr_nx/4 - i) .le. 2) then
      !print '(4(F8.3, "    "))', gr_xCoord(i), gr_U(:,i)
    !end if
  
    ! forward euler 
    gr_U(:,i) = gr_U(:,i) - dtx*(gr_flux(:,i) - gr_flux(:,i-1))
  end do
  !print *, 'after update'

  !print *, '   x          rho        rhoU       rhoE +p '
  !print *, '________    ________    ________    ________'
  !print '(4(F8.3, "    "))', gr_xCoord(pidx1), gr_U(:,pidx1)

  do i = gr_ibeg, gr_iend
    ! Eos is automatically called inside cons2prim
    !if (abs(gr_ibeg + gr_nx/4 - i) .le. 2) then
      !print '(4(F8.3, "    "))', gr_xCoord(i), gr_U(:,i)
    !end if
    call cons2prim(gr_U(:,i),gr_V(:,i))
  end do

  !print *, 'after cons2prim'
  !print *, '   x           rho          u           p   '
  !print *, '________    ________    ________    ________'

  !do i = gr_ibeg, gr_iend
    !if (abs(gr_ibeg + gr_nx/4 - i) .le. 2) then
      !print '(4(F8.3, "    "))', gr_xCoord(i), gr_V(1:3,i)
    !end if
  !end do



  return
end subroutine soln_update
