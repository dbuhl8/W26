subroutine hll(vR,vL,Flux)

#include "definition.h"  

  use grid_data
  use primconsflux, only : prim2flux,prim2cons


  implicit none
  real, dimension(NUMB_VAR), intent(IN) :: vL,vR !prim vars
  real, dimension(NSYS_VAR), intent(OUT):: Flux 

  real, dimension(NSYS_VAR) :: FL,FR,uL,uR
  real :: sL,sR,aL,aR

  call prim2flux(vL,FL)
  call prim2flux(vR,FR)
  call prim2cons(vL,uL)
  call prim2cons(vR,uR)

  
  ! left and right sound speed a
  aL = sqrt(vL(GAMC_VAR)*vL(PRES_VAR)/vL(DENS_VAR))
  aR = sqrt(vR(GAMC_VAR)*vR(PRES_VAR)/vR(DENS_VAR))

  ! fastest left and right going velocities
  sL = min(vL(VELX_VAR) - aL,vR(VELX_VAR) - aR)
  sR = max(vL(VELX_VAR) + aL,vR(VELX_VAR) + aR)

  ! numerical flux
  if (sL >= 0.) then
     Flux(1:3) = FL(1:3)
  elseif ( (sL < 0.) .and. (sR > 0.) ) then
     Flux(1:3) = (sR*FL(1:3)-sL*FR(1:3)+sR*sL*(uR(1:3)-uL(1:3)))/(sR-sL)
  else
     Flux(1:3) = FR(1:3)
  endif

  return
end subroutine hll
