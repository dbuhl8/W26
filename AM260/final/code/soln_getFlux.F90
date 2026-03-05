subroutine soln_getFlux()

#include "definition.h"  

  use grid_data
  use sim_data

  implicit none
  integer :: i

  if (sim_riemann == 'hll') then
     do i = gr_ibeg-1, gr_iend
        call hll(gr_vR(:,i),&
               gr_vL(:,i+1),&
               gr_flux(:,i))
        !if (i .eq. gr_ibeg + gr_nx/2 -1) then
          !call hll(gr_vR(:,i),&
                 !gr_vL(:,i+1),&
                 !gr_flux(:,i),.true.)
          !print *, '   v?         rho          u           p    '
          !print *, '________    ________    ________    ________'
          !print '(A8,"    ", 3(F8.3, "    "))', ' vR  i  ', gr_vR(1:3,i)
          !print '(A8,"    ", 3(F8.3, "    "))', ' vL i+1 ', gr_vL(1:3,i+1)
        !else 
          !call hll(gr_vR(:,i),&
                 !gr_vL(:,i+1),&
                 !gr_flux(:,i),.false.)
        !end if
     enddo

  elseif (sim_riemann == 'roe') then
     do i = gr_ibeg-1, gr_iend
        !call roe(gr_vR(:,i),&
                 !gr_vL(:,i+1),&
                 !gr_flux(:,i))
        if (i .eq. gr_ibeg + gr_nx/2 -1) then
          print *, "X = ", gr_xCoord(i)
          call roe(gr_vR(:,i),&
                 gr_vL(:,i+1),&
                 gr_flux(:,i),.true.)
          !print *, '   v?         rho          u           p    '
          !print *, '________    ________    ________    ________'
          !print '(A8,"    ", 3(F8.3, "    "))', ' vR  i  ', gr_vR(1:3,i)
          !print '(A8,"    ", 3(F8.3, "    "))', ' vL i+1 ', gr_vL(1:3,i+1)
          !print '(A8,"    ", 3(F8.3, "    "))', ' F i    ', gr_flux(1:3,i)
        else 
          call roe(gr_vR(:,i),&
                 gr_vL(:,i+1),&
                 gr_flux(:,i),.false.)
        end if
     enddo
  endif

  return
end subroutine soln_getFlux
