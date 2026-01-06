module slopeLimiter
  
#include "definition.h"
  
  !use grid_data


contains

  subroutine minmod(a,b,delta)
    implicit none
    real, intent(IN) :: a, b
    real, intent(OUT) :: delta

    delta = 0.5 * (sign(1.0,a) + sign(1.0,b))*min(abs(a),abs(b))

    return
  end subroutine minmod

  
  subroutine mc(a,b,delta)
    implicit none
    real, intent(IN) :: a, b
    real, intent(OUT) :: delta

    ! STUDENTS: PLEASE FINISH THIS MC LIMITER
    !stop
    delta = (sign(1.0,a)+sign(1.0,b))*min(min(abs(a),abs(b)),abs(a+b)/4)
    
    return
  end subroutine mc

  
  subroutine vanLeer(a,b,delta)
    implicit none
    real, intent(IN) :: a, b
    real, intent(OUT) :: delta

    ! STUDENTS: PLEASE FINISH THIS VAN LEER'S LIMITER
    !stop
    delta = a*b*abs(sign(1.0,a)+sign(1.0,b))/(a+b)
    ! note that the sign argument (intended to eliminate the if statement
    ! introduces a factor of 2 in the problem, hence the missing 2* in front of
    ! a*b
    
    return
  end subroutine vanLeer
  
end module slopeLimiter
