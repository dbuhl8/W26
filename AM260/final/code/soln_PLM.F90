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
  real, dimension(NSYS_VAR) :: vecL,vecR,sigL,sigR
  integer :: kWaveNum
  real :: lambdaDtDx, ldL, ldR
  real, dimension(NUMB_VAR)  :: delV,delL,delR
  real, dimension(NUMB_WAVE) :: delW
  integer :: nVar
  

  
  ! we need conservative eigenvectors
  conservative = .false.

  do i = gr_ibeg-1, gr_iend+1

     call eigenvalues(gr_V(DENS_VAR:GAME_VAR,i),lambda)
     call left_eigenvectors (gr_V(DENS_VAR:GAME_VAR,i),conservative,leig)
     call right_eigenvectors(gr_V(DENS_VAR:GAME_VAR,i),conservative,reig)

     ! primitive limiting
     if (.not. sim_charLimiting) then
        do kWaveNum = 1, NUMB_WAVE
           ! slope limiting
           ! deltas in primitive vars
           delL(DENS_VAR:PRES_VAR) = gr_V(DENS_VAR:PRES_VAR,i  )-gr_V(DENS_VAR:PRES_VAR,i-1)
           delR(DENS_VAR:PRES_VAR) = gr_V(DENS_VAR:PRES_VAR,i+1)-gr_V(DENS_VAR:PRES_VAR,i  )
           do nVar = DENS_VAR,PRES_VAR
              if (sim_limiter == 'minmod') then
                 call minmod(delL(nVar),delR(nVar),delV(nVar))
              elseif (sim_limiter == 'vanLeer') then
                 call vanLeer(delL(nVar),delR(nVar),delV(nVar))
              elseif (sim_limiter == 'mc') then
                 call mc(delL(nVar),delR(nVar),delV(nVar))
              endif
           enddo
           ! project primitive delta to characteristic vars
           delW(kWaveNum) = dot_product(leig(DENS_VAR:PRES_VAR,kWaveNum),delV(DENS_VAR:PRES_VAR))
        enddo
     elseif (sim_charLimiting) then
        !stop
        !STUDENTS: PLEASE FINISH THIS CHARACTERISTIC LIMITING
        !(THE IMPLEMENTATION SHOULD NOT BE LONGER THAN THE PRIMITIVE LIMITING CASE)
        do kWaveNum = 1, NUMB_WAVE
           ! slope limiting
           ! deltas in primitive vars
           delL(DENS_VAR:PRES_VAR) = gr_V(DENS_VAR:PRES_VAR,i  )-gr_V(DENS_VAR:PRES_VAR,i-1)
           delR(DENS_VAR:PRES_VAR) = gr_V(DENS_VAR:PRES_VAR,i+1)-gr_V(DENS_VAR:PRES_VAR,i  )
           ldL = dot_product(leig(DENS_VAR:PRES_VAR,kWaveNum),&
                  delL(DENS_VAR:PRES_VAR))
           ldR = dot_product(leig(DENS_VAR:PRES_VAR,kWaveNum),&
                  delR(DENS_VAR:PRES_VAR))
           !do nVar = DENS_VAR,PRES_VAR
           if (sim_limiter == 'minmod') then
              call minmod(ldL, ldR, delW(kWaveNum))
              !call minmod(delL(nVar),delR(nVar),delV(nVar))
              !call minmod(&
                !dot_product(leig(DENS_VAR:PRES_VAR,kWaveNum),&
                  !delL(DENS_VAR:PRES_VAR)),&
                !dot_product(leig(DENS_VAR:PRES_VAR,kWaveNum),&
                  !delR(DENS_VAR:PRES_VAR)),&
                !delW(kWaveNum))
           elseif (sim_limiter == 'vanLeer') then
              !call vanLeer(delL(nVar),delR(nVar),delV(nVar))
              call vanLeer(ldL, ldR, delW(kWaveNum))
              !call vanLeer(&
                !dot_product(leig(DENS_VAR:PRES_VAR,kWaveNum),&
                  !delL(DENS_VAR:PRES_VAR)),&
                !dot_product(leig(DENS_VAR:PRES_VAR,kWaveNum),&
                  !delR(DENS_VAR:PRES_VAR)),&
                !delW(kWaveNum))
           elseif (sim_limiter == 'mc') then
              !call mc(delL(nVar),delR(nVar),delV(nVar))
              call mc(ldL,ldR,delW(kWaveNum))
              !call mc(&
                !dot_product(leig(DENS_VAR:PRES_VAR,kWaveNum),&
                  !delL(DENS_VAR:PRES_VAR)),&
                !dot_product(leig(DENS_VAR:PRES_VAR,kWaveNum),&
                  !delR(DENS_VAR:PRES_VAR)),&
                !delW(kWaveNum))
           endif
           !enddo
           !delW(kWaveNum) = dot_product(leig(DENS_VAR:PRES_VAR,kWaveNum),delV(DENS_VAR:PRES_VAR))
           ! project primitive delta to characteristic vars
        enddo

     endif



     ! set the initial sum to be zero
     sigL(DENS_VAR:ENER_VAR) = 0.
     sigR(DENS_VAR:ENER_VAR) = 0.
     vecL(DENS_VAR:ENER_VAR) = 0.
     vecR(DENS_VAR:ENER_VAR) = 0.
     
     do kWaveNum = 1, NUMB_WAVE
        lambdaDtDx = lambda(kWaveNum)*dt/gr_dx

        
        if (sim_riemann == 'roe') then
           !stop
           ! STUDENTS: PLEASE FINISH THIS ROE SOLVER CASE
           ! THIS SHOULDN'T BE LONGER THAN THE HLL CASE
           if (lambdaDtDx .gt. 0) then
             vecR(1:3) = 0.5*(1.0 - lambdaDtDx)*reig(1:3,kWaveNum)*delW(kWaveNum)
             sigR(1:3) = sigR(1:3) + vecR(1:3)
           else
             vecL(1:3) = 0.5*(-1.0 - lambdaDtDx)*reig(1:3,kWaveNum)*delW(kWaveNum)
             sigL(1:3) = sigL(1:3) + vecL(1:3)
           end if

        elseif (sim_riemann == 'hll') then
              vecR(1:3) = 0.5*(1.0 - lambdaDtDx)*reig(1:3,kWaveNum)*delW(kWaveNum)
              sigR(1:3) = sigR(1:3) + vecR(1:3)

              vecL(1:3) = 0.5*(-1.0 - lambdaDtDx)*reig(1:3,kWaveNum)*delW(kWaveNum)
              sigL(1:3) = sigL(1:3) + vecL(1:3)
        endif

        ! Let's make sure we copy all the cell-centered values to left and right states
        ! this will be just FOG
        !gr_vL(DENS_VAR:NUMB_VAR,i) = gr_V(DENS_VAR:NUMB_VAR,i)
        !gr_vR(DENS_VAR:NUMB_VAR,i) = gr_V(DENS_VAR:NUMB_VAR,i)
        
        ! Now PLM reconstruction for dens, velx, and pres
        gr_vL(1:3,i) = gr_V(1:3,i) + sigL(1:3)
        gr_vR(1:3,i) = gr_V(1:3,i) + sigR(1:3)
     end do
  end do
  
  return
end subroutine soln_PLM
