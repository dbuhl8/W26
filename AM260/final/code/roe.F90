subroutine roe(vR,vL,Flux)

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
  integer :: kWaveNum
  
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


  do kWaveNum = 1, NUMB_WAVE
     ! STUDENTS: PLEASE FINISH THIS ROE SOLVER
     !stop
    sigma = sigma + dot_product(leig(:,kWaveNum),uR-uL)*abs(lambda(kWaveNum))&
      *reig(:,kWaveNum)
  end do
  
  ! numerical flux
  Flux = 0.5*(FL + FR) - 0.5*sigma


  return
end subroutine roe
