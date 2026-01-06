subroutine averageState(vL,vR,vAvg)

#include "definition.h"  

  implicit none
  real, dimension(NUMB_VAR), intent(IN) :: vL,vR !prim vars
  real, dimension(NUMB_VAR), intent(OUT) :: vAvg  !average state

  ! STUDENTS: PLEASE FINISH THIS SIMPLE AVERAGING SCHEME
  vAvg = (vL + vR)/2
   
  return
end subroutine averageState
