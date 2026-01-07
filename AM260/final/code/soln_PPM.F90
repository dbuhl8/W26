subroutine soln_PPM(dt)

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
  real, dimension(NSYS_VAR) :: vecL,vecR,sigL,sigR,sig2L,sig2R,vec2L,vec2R
  integer :: kWaveNum, pt_idx = 64
  real :: lambdaDtDx, delC1, delC2, alpha, beta, csL, csR
  real, dimension(3)  :: delV,delL,delR,aLR,delLL,delRR,dVL,dVC,dVR
  real, dimension(3,3) :: C
  real, dimension(NUMB_WAVE) :: delW
  integer :: nVar
  

  ! we need conservative eigenvectors
  conservative = .false.

  do i = gr_ibeg-1, gr_iend+1
    ! step 1: parabolic profile
    C = 0.
    delLL(1:3) = gr_V(1:3,i-1) &
      -gr_V(1:3,i-2)
    delL(1:3)  = gr_V(1:3,i  ) &
      -gr_V(1:3,i-1)
    delR(1:3)  = gr_V(1:3,i+1) &
      -gr_V(1:3,i  )
    delRR(1:3) = gr_V(1:3,i+2) &
      -gr_V(1:3,i+1)

    if (i .eq. pt_idx) then
      print *,  " @ i = pt_idx, x = ", gr_xCoord(i)
      print *,  ""

      print *,  " DEL L @ Line 36"
      print *,  delL
      print *,  ""

      print *,  " DEL LL @ Line 36"
      print *,  delLL
      print *,  ""

      print *,  " DEL RR @ Line 36"
      print *,  delRR
      print *,  ""

      print *,  " DEL R @ Line 36"
      print *,  delR
      print *,  ""
    endif

    do nVar = DENS_VAR, PRES_VAR
      if (sim_limiter == 'minmod') then
         call minmod(delR(nVar),delRR(nVar),dvR(nVar))
         call minmod(delL(nVar),delR(nVar),dvC(nVar))
         call minmod(delLL(nVar),delL(nVar),dvL(nVar))
      elseif (sim_limiter == 'vanLeer') then
         !call vanLeer(delL(nVar),delR(nVar),delV(nVar))
         call vanLeer(delR(nVar),delRR(nVar),dvR(nVar))
         call vanLeer(delL(nVar),delR(nVar),dvC(nVar))
         call vanLeer(delLL(nVar),delL(nVar),dvL(nVar))

      elseif (sim_limiter == 'mc') then
         !call mc(delL(nVar),delR(nVar),delV(nVar))
         call mc(delR(nVar),delRR(nVar),dvR(nVar))
         call mc(delL(nVar),delR(nVar),dvC(nVar))
         call mc(delLL(nVar),delL(nVar),dvL(nVar))
      endif
    end do 
    
    ! this comes from equation 9.47 on page 158 of the typed lect note
    ! a+/- = (1/2)(v(i-1) + v(i)) - (1/6)(dv(i) - dv(i-1))
    ! vL = a-
    ! vR = a+
    vecL(1:3) = 0.5*(gr_V(1:3,i-1)+gr_V(1:3,i)) - (1./6)*&
      (dvC(1:3) - dvL(1:3))
    vecR(1:3) = 0.5*(gr_V(1:3,i)+gr_V(1:3,i+1)) - (1./6)*&
      (dvR(1:3) - dvC(1:3))

    if (i .eq. pt_idx) then
      print *,  " VEC L @ Line 78"
      print *,  vecL
      print *,  ""
      print *,  " VEC R @ Line 78"
      print *,  vecR
      print *,  ""
    endif

    C(:,3) = (6./(gr_dx**2))*(0.5*(vecL(1:3)+vecR(1:3))&
     - gr_V(1:3,i))
    C(:,2) = (vecR(1:3)-vecL(1:3))/gr_dx
    C(:,1) = gr_V(1:3,i) - C(:,3)*(gr_dx**2)/12.

    if (i .eq. pt_idx) then
      print *,  " C @ Line 92"
      print *, C(:,:)
      print *,  " "
    endif

    ! NEED TO IMPLEMENT CONDITION 1 on Page 159 of the lecture note. 
    ! it will probably be easier to just return the FOG method rather than shift
    ! the parameters vecL, vecR, C
    do nVar = DENS_VAR, PRES_VAR 
      if ((vecR(nvar) - gr_V(nvar,i))*(-vecL(nvar) + gr_V(nvar,i)) .le. 0) then
        gr_vL(:,i) = gr_V(:,i)
        gr_vR(:,i) = gr_V(:,i)
        return
      else
        alpha = (vecR(nVar)-vecL(nVar))**2
        beta = 6*(vecR(nVar)-vecL(nVar))*&
          (gr_V(nVar,i)-(vecR(nVar)+vecL(nVar))/2)
        if (-alpha .gt. beta) then
          vecR = 3*gr_V(:,i) - 2*vecL
          C(:,3) = (6./(gr_dx**2))*(0.5*(vecL(1:3)+&
            vecR(1:3))&
           - gr_V(1:3,i))
          C(:,2) = (vecR(1:3)-vecL(1:3))/gr_dx
          C(:,1) = gr_V(1:3,i) - C(:,3)*(gr_dx**2)/12.
          exit
        elseif (alpha .lt. beta) then
          vecL = 3*gr_V(:,i) - 2*vecR
          C(:,3) = (6./(gr_dx**2))*(0.5*(vecL(1:3)+&
            vecR(1:3)) - gr_V(1:3,i))
          C(:,2) = (vecR(1:3)-vecL(1:3))/gr_dx
          C(:,1) = gr_V(1:3,i) - C(:,3)*(gr_dx**2)/12.
          exit
        endif
      endif
    end do 


    if (i .eq. pt_idx) then
      print *,  " C @ Line 128"
      print *, C(:,:)
      print *,  " "
    endif

    ! step 2: characteristic tracing to compute the half step integration
    call eigenvalues(gr_V(:,i),lambda)
    call left_eigenvectors (gr_V(:,i),conservative,leig)
    call right_eigenvectors(gr_V(:,i),conservative,reig)

    ! set the initial sum to be zero
    sigL(DENS_VAR:ENER_VAR) = 0.
    sigR(DENS_VAR:ENER_VAR) = 0.
    vecL(DENS_VAR:ENER_VAR) = 0.
    vecR(DENS_VAR:ENER_VAR) = 0.
    sig2L(DENS_VAR:ENER_VAR) = 0.
    sig2R(DENS_VAR:ENER_VAR) = 0.
    vec2L(DENS_VAR:ENER_VAR) = 0.
    vec2R(DENS_VAR:ENER_VAR) = 0.
   
    do kWaveNum = 1, NUMB_WAVE
      ! lambdaDtDx = lambda*dt/dx
      lambdaDtDx = lambda(kWaveNum)*dt/gr_dx
      delC1 = gr_dx*dot_product(leig(:,kWaveNum),C(:,2))
      delC2 = (gr_dx**2)*dot_product(leig(:,kWaveNum),C(:,3))
      if (lambdaDtDx .gt. 0) then
        vecR = 0.5*(1.0 - lambdaDtDx)&
          *reig(:,kWaveNum)*delC1
        vec2R = 0.25*(1.0 - 2*lambdaDtDx+(4./3)&
          *lambdaDtDx**2)*&
          reig(:,kWaveNum)*delC2
        sigR = sigR + vecR
        sig2R = sig2R + vec2R
      else
        vecL = 0.5*(-1.0 - lambdaDtDx)*&
          reig(:,kWaveNum)*delC1
        vec2L(:) = 0.25*(1.0 + 2*lambdaDtDx+(4./3)&
          *lambdaDtDx**2)*reig(:,kWaveNum)*delC2
        sigL = sigL + vecL
        sig2L = sig2L + vec2L
      endif

      ! NOTE: I dont believe that PPM is concerned whether the fluxes are
      ! comptued using ROE or HLL so this section has been commented out and the
      ! core of the ROE section is used above.       
      !if (sim_riemann == 'roe') then
        !if (lambdaDtDx .gt. 0) then
          !vecR(1:3) = 0.5*(1.0 - lambdaDtDx)&
            !*reig(1:3,kWaveNum)*delC1
          !vec2R(1:3) = 0.25*(1.0 - 2*lambdaDtDx+(4./3)&
            !*lambdaDtDx**2)*&
            !reig(1:3,kWaveNum)*delC2
          !sigR(1:3) = sigR(1:3) &
            !+ vecR(1:3)
          !sig2R(1:3) = sig2R(1:3) &
            !+ vec2R(1:3)
        !else
          !vecL(1:3) = 0.5*(-1.0 - lambdaDtDx)*&
            !reig(1:3,kWaveNum)*delC1
          !vec2L(1:3) = 0.25*(1.0 + 2*lambdaDtDx+(4./3)&
            !*lambdaDtDx**2)*reig(1:3,kWaveNum)*delC2
          !sigL(1:3) = sigL(1:3) &
            !+ vecL(1:3)
          !sig2L(1:3) = sig2L(1:3) &
            !+ vec2L(1:3)
        !endif
      !elseif (sim_riemann == 'hll') then
        !vecR(1:3) = 0.5*(1.0 - lambdaDtDx)*&
          !reig(1:3,kWaveNum)*delW(kWaveNum)
        !sigR(1:3) = sigR(1:3) &
          !+ vecR(1:3)
        !vecL(1:3) = 0.5*(-1.0 - lambdaDtDx)&
          !*reig(1:3,kWaveNum)*delW(kWaveNum)
        !sigL(1:3) = sigL(1:3) &
          !+ vecL(1:3)
      !endif

      ! Now PPM reconstruction for dens, velx, and pres
      gr_vL(1:3,i) = C(:,1) + sigL + sig2L
      gr_vR(1:3,i) = C(:,1) + sigR + sig2R
      gr_vL(EINT_VAR,i) = sqrt(gr_vL(GAME_VAR,i)*gr_vL(PRES_VAR,i)/gr_vL(DENS_VAR,i))
      gr_vR(EINT_VAR,i) = sqrt(gr_vR(GAME_VAR,i)*gr_vR(PRES_VAR,i)/gr_vR(DENS_VAR,i))
   end do
 end do
 
  return
end subroutine soln_PPM
