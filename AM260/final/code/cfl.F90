subroutine cfl(dt)

#include "definition.h"
  
  use grid_data
  use sim_data, only: sim_cfl

  implicit none
  real, intent(OUT) :: dt
  integer :: i
  real :: maxSpeed, lambda, cs, maxDt

  maxSpeed = 1e-8
  maxDt = 0.01
  !! update conservative vars
  do i = gr_ibeg, gr_iend
    ! this needs to be changes to include all wave speeds (pressure, dens,
    ! and velx)
     cs = sqrt(gr_V(GAMC_VAR,i)*gr_V(PRES_VAR,i)/gr_V(DENS_VAR,i))
     lambda=(maxval(abs(gr_V(1:3,i))) + cs)
     maxSpeed=max(maxSpeed,lambda)
  end do

  ! cfl
  dt = sim_cfl*gr_dx/maxSpeed
  dt = min(dt,maxDt)

  return

end subroutine cfl
