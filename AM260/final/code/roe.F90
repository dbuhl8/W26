!subroutine roe(vL,vR,Flux)
subroutine roe(vL,vR,Flux,pt)

#include "definition.h"  

  use grid_data
  use primconsflux, only : prim2flux,prim2cons
  use eigensystem

  implicit none
  real, dimension(NUMB_VAR), intent(IN) :: vL,vR !prim vars
  real, dimension(NSYS_VAR), intent(OUT):: Flux 

  real, dimension(NSYS_VAR)  :: FL,FR,uL,uR
  real, dimension(NUMB_VAR)  :: vAvg
  real, dimension(NUMB_WAVE) :: lambda
  real, dimension(NSYS_VAR,NUMB_WAVE) :: reig, leig
  logical :: conservative
  real, dimension(NSYS_VAR) :: vec, sigma
  integer :: k, pidx1, pidx2
  logical, intent(in) :: pt
  
  ! set the initial sum to be zero
  sigma = 0.
  vec = 0.
  
  ! we need conservative eigenvectors
  conservative = .true.


  call averageState(vL,vR,vAvg)
  call eigenvalues(vAvg,lambda)
  call left_eigenvectors (vAvg,conservative,leig)
  call right_eigenvectors(vAvg,conservative,reig)
  
  call prim2flux(vL,FL)
  call prim2flux(vR,FR)
  call prim2cons(vL,uL)
  call prim2cons(vR,uR)


  do k = 1, NUMB_WAVE
     ! STUDENTS: PLEASE FINISH THIS ROE SOLVER
     !stop
    sigma = sigma + dot_product(leig(:,k),uR-uL)*abs(lambda(k))*reig(:,k)
  end do

  ! numerical flux
  Flux = 0.5*(FL + FR) - 0.5*sigma

  if (pt) then
    print '(A, 3(F8.3, "    "))', 'VAVG   :', vAvg(1:3)
    print '(A, 3(F8.3, "    "))', 'VL     :', vL(1:3)
    print '(A, 3(F8.3, "    "))', 'VR     :', vR(1:3)
    print '(A, 3(F8.3, "    "))', 'UL     :', uL(1:3)
    print '(A, 3(F8.3, "    "))', 'UR     :', uR(1:3)
    print '(A, 3(F8.3, "    "))', 'FL     :', FL(1:3)
    print '(A, 3(F8.3, "    "))', 'FR     :', FR(1:3)
    print '(A, 3(F8.3, "    "))', 'Flux   :', Flux(1:3)
  end if



  return
end subroutine roe
