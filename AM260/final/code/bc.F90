module bc

#include "definition.h"

  use grid_data
  use sim_data, only : sim_bcType
  implicit none

contains

  subroutine bc_apply()
    implicit none
    if (sim_bcType == 'outflow') then
       call bc_outflow(gr_V)
    elseif (sim_bcType == 'reflect') then
       call bc_reflect(gr_V)
    elseif (sim_bcType == 'periodic') then
    elseif (sim_bcType == 'user') then

       ! STUDENTS: PLEASE IMPLEMENT THIS FOR THE SHU-OSHER PROBLEM
       stop
    endif

  end subroutine bc_apply


  subroutine bc_outflow(V)
    implicit none
    real, dimension(NUMB_VAR,gr_imax), intent(INOUT) :: V
    integer :: i

    do i = 1, gr_ngc
       ! on the left GC
       V(1:NUMB_VAR,i) = V(1:NUMB_VAR,gr_ibeg) 

       ! on the right GC
       V(1:NUMB_VAR,gr_iend+i) = V(1:NUMB_VAR,gr_iend)
    end do

    return
  end subroutine bc_outflow



  subroutine bc_reflect(V)
    implicit none
    real, dimension(NUMB_VAR,gr_imax), intent(INOUT) :: V
    integer :: i,k0,k1

    do i = 1, gr_ngc
       !k0 = 2*gr_ngc+1
       !k1 = gr_iend-gr_ngc
      k0 = gr_iend

       ! on the left GC
       V(       :,i) = V(       :,k0-i)

       ! on the right GC
       V(       :,k1+k0-i) = V(         :,k1+i)
    end do

    return
  end subroutine bc_reflect

  subroutine bc_periodic(V)
    implicit none
    real, dimension(NUMB_VAR,gr_imax), intent(INOUT) :: V
    real :: k0, k1

    k0 = 2*gr_ngc+1
    k1 = gr_iend-gr_ngc

    V(:,1:gr_ngc) = V(:,k1:gr_iend)
    V(:,k1:) = V(:,gr_ngc+1:2*gr_ngc)

    return
  end subroutine bc_periodic

  subroutine bc_user(V)
    implicit none
    real, dimension(NUMB_VAR,gr_imax), intent(INOUT) :: V
    ! STUDENTS: PLEASE IMPLEMENT THIS FOR THE SHU-OSHER PROBLEM
    V(DENS_VAR,:gr_ngc) = 3.857143
    V(VELX_VAR,:gr_ngc) = 2.629369
    V(PRES_VAR,:gr_ngc) = 10.33333

    V(DENS_VAR,gr_iend:) = 1 + 0.2*sin(5.0*gr_xCoord(gr_iend:))
    V(VELX_VAR,gr_iend:) = 0.0
    V(PRES_VAR,gr_iend:) = 1.0
    return
  end subroutine bc_user

end module bc
