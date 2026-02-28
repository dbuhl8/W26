subroutine hll(vL,vR,Flux)
!subroutine hll(vL,vR,Flux,pt)

#include "definition.h"  

  use grid_data
  use sim_data, only: sim_smallpres
  use primconsflux, only : prim2flux,prim2cons


  implicit none
  real, dimension(NUMB_VAR), intent(IN) :: vL,vR !prim vars
  real, dimension(NSYS_VAR), intent(OUT):: Flux 

  real, dimension(NSYS_VAR) :: FL,FR,uL,uR
  real :: sL,sR,aL,aR, pres
  !logical, intent(IN) :: pt

  call prim2flux(vL,fL)
  call prim2flux(vR,fR)
  call prim2cons(vL,uL)
  call prim2cons(vR,uR)

  
  ! left and right sound speed a
  aL = sqrt(vL(GAME_VAR)*max(vL(PRES_VAR),sim_smallpres)/vL(DENS_VAR))
  aR = sqrt(vR(GAME_VAR)*max(vR(PRES_VAR),sim_smallpres)/vR(DENS_VAR))

  ! fastest left and right going velocities
  sL = min(vL(VELX_VAR) - aL,vR(VELX_VAR) - aR)
  sR = max(vL(VELX_VAR) + aL,vR(VELX_VAR) + aR)

  !if (pt) then
    !print '(A, F8.3)', 'aL: ', aL
    !print '(A, F8.3)', 'aR: ', aR
    !print '(A, F8.3)', 'sL: ', sL
    !print '(A, F8.3)', 'sR: ', sR
    !print '(A, 3(F8.3, "    "))', 'fL: ', fL
    !print '(A, 3(F8.3, "    "))', 'fR: ', fR
  !end if

  ! numerical flux
  if (sL > 0.) then
     Flux = fL
  elseif ( (sL < 0.) .and. (sR >= 0.) ) then
     Flux = (sR*fL-sL*fR+sR*sL*(uR-uL))/(sR-sL)
  else
     Flux = fR
  endif

  return
end subroutine hll
