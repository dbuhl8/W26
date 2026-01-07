module eigensystem

#include "definition.h"
  
  use grid_data

contains

  subroutine eigenvalues(V,lambda)
    implicit none

    real, dimension(NUMB_VAR), intent(IN)  :: V
    real, dimension(NUMB_WAVE),intent(OUT) :: lambda

    real :: a, u

    ! sound speed
    a = sqrt(V(GAMC_VAR)*V(PRES_VAR)/V(DENS_VAR))
    u = V(VELX_VAR)
    
    lambda(SHOCKLEFT) = u - a
    lambda(CTENTROPY) = u
    lambda(SHOCKRGHT) = u + a
    
    return
  end subroutine eigenvalues


  
  subroutine right_eigenvectors(V,conservative,reig)
    implicit none
    real, dimension(NUMB_VAR), intent(IN)  :: V
    logical :: conservative
    real, dimension(NSYS_VAR,NUMB_WAVE), intent(OUT) :: reig

    real :: cs, u, rho, pres, gam, ke
    
    ! sound speed, and others
    u = V(VELX_VAR)
    pres = V(PRES_VAR)
    rho = V(DENS_VAR)
    gam = V(GAMC_VAR)
    cs = sqrt(gam*pres/rho)
    ke = 0.5*u**2
    
    if (conservative) then
       !! Conservative eigenvector
       reig(DENS_VAR,SHOCKLEFT) = 1.
       reig(VELX_VAR,SHOCKLEFT) = u - cs
       reig(PRES_VAR,SHOCKLEFT) = ke + cs**2/(gam-1) - cs*u
       reig(:,SHOCKLEFT) = -0.5*rho/cs*reig(:,SHOCKLEFT)

       reig(DENS_VAR,CTENTROPY) = 1.
       reig(VELX_VAR,CTENTROPY) = u
       reig(PRES_VAR,CTENTROPY) = ke
       
       reig(DENS_VAR,SHOCKRGHT) = 1.
       reig(VELX_VAR,SHOCKRGHT) = u + a
       reig(PRES_VAR,SHOCKRGHT) = ke + cs**2/(gam-1) + cs*u
       reig(:,SHOCKRGHT) = 0.5*rho/cs*reig(:,SHOCKRGHT)

  else
  !! Primitive eigenvector
  !! STUDENTS: PLEASE FINISH THIS PRIMITIVE RIGHT EIGEN VECTORS
  !print*,'eigensystem.F90: right eigenvectors'
  !stop
       reig(DENS_VAR,SHOCKLEFT) = -0.5*rho/cs
       reig(VELX_VAR,SHOCKLEFT) = 0.5
       reig(PRES_VAR,SHOCKLEFT) = -0.5*rho*cs

       reig(DENS_VAR,CTENTROPY) = 1.
       reig(VELX_VAR,CTENTROPY) = 0
       reig(PRES_VAR,CTENTROPY) = 0
       
       reig(DENS_VAR,SHOCKRGHT) = 0.5*rho/cs
       reig(VELX_VAR,SHOCKRGHT) = 0.5
       reig(PRES_VAR,SHOCKRGHT) = 0.5*rho*cs
  endif

    
    return
  end subroutine right_eigenvectors


  subroutine left_eigenvectors(V,conservative,leig)
    implicit none
    real, dimension(NUMB_VAR), intent(IN)  :: V
    logical :: conservative
    real, dimension(NSYS_VAR,NUMB_WAVE), intent(OUT) :: leig

    real :: a, u, d, g, ke, hdai, hda
    
    ! sound speed, and others
    u = V(VELX_VAR)
    pres = V(PRES_VAR)
    rho = V(DENS_VAR)
    gam = V(GAMC_VAR)
    cs = sqrt(gam*pres/rho)
    ke = 0.5*u**2
    
    if (conservative) then
       !! Conservative eigenvector
       !! STUDENTS: PLEASE FINISH THIS CONSERVATIVE LEFT EIGEN VECTORS
       !print*,'eigeysystem.F90: left conservative eigenvectors'
       !stop
       leig(DENS_VAR,SHOCKLEFT) = -ke - cs*u/(gam-1)
       leig(VELX_VAR,SHOCKLEFT) = u + cs/(gam-1) 
       leig(PRES_VAR,SHOCKLEFT) = -1.

       leig(DENS_VAR,CTENTROPY) = (rho/cs)*(-ke + (cs**2)/(gam-1))
       leig(VELX_VAR,CTENTROPY) = rho*u/cs
       leig(PRES_VAR,CTENTROPY) = -rho/(gam-1)
       
       leig(DENS_VAR,SHOCKRGHT) = ke - cs*u/(gam-1)
       leig(VELX_VAR,SHOCKRGHT) = -u + cs/(gam-1) 
       leig(PRES_VAR,SHOCKRGHT) = 1.
       !leig(:,:) = ((gam-1)/(rho*cs))*leig(:,:)
       
    else
       !! Primitive eigenvector
       !! STUDENTS: PLEASE FINISH THIS PRIMITIVE LEFT EIGEN VECTORS
       !print*,'eigeysystem.F90: left prim eigenvectors'
       !stop
       leig(DENS_VAR,SHOCKLEFT) = 0
       leig(VELX_VAR,SHOCKLEFT) = 1.
       leig(PRES_VAR,SHOCKLEFT) = -1./(rho*cs)

       leig(DENS_VAR,CTENTROPY) = 1.
       leig(VELX_VAR,CTENTROPY) = 0
       leig(PRES_VAR,CTENTROPY) = -1./(cs**2)
       
       leig(DENS_VAR,SHOCKRGHT) = 0
       leig(VELX_VAR,SHOCKRGHT) = 1.
       leig(PRES_VAR,SHOCKRGHT) = 1./(rho*cs)

    endif
    
    return
  end subroutine left_eigenvectors


  
end module eigensystem
