subroutine soln_PLM(dt)
#include "definition.h"  

  use grid_data
  use sim_data
  use slopeLimiter
  use eigensystem

  implicit none
  real, intent(IN) :: dt
  integer :: i

  real, dimension(NUMB_WAVE) :: lambda
  real, dimension(NSYS_VAR,NUMB_WAVE) :: reig, leig
  logical :: conservative
  real, dimension(NSYS_VAR) :: sigL,sigR
  integer :: kWaveNum, k, pidx1, pidx2
  real :: lambdaDtDx, ldL, ldR
  real, dimension(NUMB_VAR)  :: delV,delL,delR
  real, dimension(NUMB_WAVE) :: delW
  integer :: nVar

  pidx1 = gr_ibeg + gr_nx/2 -1
  pidx2 = gr_ibeg + gr_nx/2

  ! we need conservative eigenvectors
  conservative = .true.

  do i = gr_ibeg-1, gr_iend+1



    call eigenvalues(gr_V(DENS_VAR:GAME_VAR,i),lambda)
    call left_eigenvectors (gr_V(DENS_VAR:GAME_VAR,i),conservative,leig)
    call right_eigenvectors(gr_V(DENS_VAR:GAME_VAR,i),conservative,reig)

    ! primitive limiting
    if (.not. sim_charLimiting) then
      delL(1:3) = gr_V(1:3,i  )-gr_V(1:3,i-1)
      delR(1:3) = gr_V(1:3,i+1)-gr_V(1:3,i  )
      do nVar = DENS_VAR,PRES_VAR
        if (sim_limiter == 'minmod') then
          call minmod(delL(nVar),delR(nVar),delV(nVar))
        elseif (sim_limiter == 'vanLeer') then
          call vanLeer(delL(nVar),delR(nVar),delV(nVar))
        elseif (sim_limiter == 'mc') then
          call mc(delL(nVar),delR(nVar),delV(nVar))
        endif
      enddo
      do k = 1, NUMB_WAVE
        delW(k) = dot_product(leig(1:3,k),delV(1:3))
      enddo
    elseif (sim_charLimiting) then
      !stop
      !STUDENTS: PLEASE FINISH THIS CHARACTERISTIC LIMITING
      !(THE IMPLEMENTATION SHOULD NOT BE LONGER THAN THE PRIMITIVE LIMITING CASE)
      delL(1:3) = gr_V(1:3,i  )-gr_V(1:3,i-1)
      delR(1:3) = gr_V(1:3,i+1)-gr_V(1:3,i  )
      do k = 1, NUMB_WAVE
        ! slope limiting
        ldL = dot_product(leig(1:3,k),delL(1:3))
        ldR = dot_product(leig(1:3,k),delR(1:3))
        if (sim_limiter == 'minmod') then
          call minmod(ldL, ldR, delW(k))
        elseif (sim_limiter == 'vanLeer') then
          call vanLeer(ldL, ldR, delW(k))
        elseif (sim_limiter == 'mc') then
          call mc(ldL,ldR,delW(k))
        endif
      enddo
    endif

    ! set the initial sum to be zero
    sigL(1:3) = 0.
    sigR(1:3) = 0.

    do k = 1, NUMB_WAVE
      lambdaDtDx = lambda(k)*dt/gr_dx
      !if (sim_riemann == 'roe') then
      if (lambdaDtDx .gt. 0) then
        sigR(1:3) = sigR(1:3) + 0.5*(1.0 - lambdaDtDx)*reig(1:3,k)*delW(k)
      else
        sigL(1:3) = sigL(1:3) + 0.5*(-1.0 - lambdaDtDx)*reig(1:3,k)*delW(k)
      end if
      !elseif (sim_riemann == 'hll') then
        !sigR(1:3) = sigR(1:3) + 0.5*(1.0 - lambdaDtDx)*reig(1:3,k)*delW(k)
        !sigL(1:3) = sigL(1:3) + 0.5*(-1.0 - lambdaDtDx)*reig(1:3,k)*delW(k)
      !endif
    end do
    ! Now PLM reconstruction for dens, velx, and pres
    gr_vL(1:3,i) = gr_V(1:3,i) + sigL(1:3)
    gr_vR(1:3,i) = gr_V(1:3,i) + sigR(1:3)
    gr_vL(4:6,i) = gr_V(4:6,i)
    gr_vR(4:6,i) = gr_V(4:6,i)
    gr_vL(3,i) = max(gr_vL(3,i), sim_smallpres)
    gr_vR(3,i) = max(gr_vR(3,i), sim_smallpres)
    !gr_vL(1,i) = max(gr_vL(1,i), sim_smallpres)
    !gr_vR(1,i) = max(gr_vR(1,i), sim_smallpres)
    !if (i .eq. pidx1) then
      !print '(A, F8.3)', 'x: ', gr_xCoord(i)
      !print '(A, 3(F8.3, "    "))', 'dw     :', delW(:)
      !print '(A, 3(F8.3, "    "))', 'lambda :', lambda(:)
      !print '(A, 3(F8.3, "    "))', 'vL     :', gr_vL(1:3,i)
      !print '(A, 3(F8.3, "    "))', 'vR     :', gr_vR(1:3,i)
      !print '(A, 3(F8.3, "    "))', 'v      :', gr_V(1:3,i)
      !print '(A, 3(F8.3, "    "))', 'sigL   :', sigL(1:3)
      !print '(A, 3(F8.3, "    "))', 'sigR   :', sigR(1:3)
    !elseif (i .eq. pidx2) then
      !print '(A, F8.3)', 'x: ', gr_xCoord(i)
      !print '(A, 3(F8.3, "    "))', 'dw     :', delW(:)
      !print '(A, 3(F8.3, "    "))', 'lambda :', lambda(:)
      !print '(A, 3(F8.3, "    "))', 'vL     :', gr_vL(1:3,i)
      !print '(A, 3(F8.3, "    "))', 'vR     :', gr_vR(1:3,i)
      !print '(A, 3(F8.3, "    "))', 'v      :', gr_V(1:3,i)
      !print '(A, 3(F8.3, "    "))', 'sigL   :', sigL(1:3)
      !print '(A, 3(F8.3, "    "))', 'sigR   :', sigR(1:3)
    !end if

  end do

  return
end subroutine soln_PLM
