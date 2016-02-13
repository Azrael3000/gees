!=====================================================
!
! Copyright 2012-2016 Arno Mayrhofer
! This program is licensed under the GNU General Public License.
!
! This file is part of Gees.
!
! Gees is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Gees is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Gees.  If not, see <http://www.gnu.org/licenses/>.
!
!=====================================================

module mod_fluxcalc

  use mod_helper

  ! physical parameters
  real(dp), parameter :: rho0 = 1d3
  real(dp), parameter :: xi = 7d0

  ! public variables
  public :: dp
  ! public functions
  public :: fluxcalc, bcs
  ! private functions
  private :: phi, ev, p_eos

contains

  subroutine fluxcalc(u, f, nx, c0)
    ! Task:
    !   Compute flux at midpoints
    ! Input:
    !   u  - state vector at center
    !   nx - number of gridpoints
    !   c0 - speed of sound
    ! Output:
    !   f  - fluxes at midpoints

    implicit none

    integer, intent(in) :: nx
    real(dp), dimension(-1:nx+1,2), intent(inout) :: u
    real(dp), dimension(nx,2), intent(inout) :: f
    real(dp), intent(in) :: c0

    integer :: i, j
    real(dp), dimension(2) :: ur, ul
    real(dp) :: a, pl, pr, nom, denom, ri, rim1
    real(dp), dimension(4) :: lam

    ! calculate flux at midpoints, f(i) = f_{i-1/2}
    do i=1,nx
      do j=1,2
        ! calculate r_{i} = \frac{u_i-u_{i-1}}{u_{i+1}-u_i}
        nom   = (u(i  ,j)-u(i-1,j))
        denom = (u(i+1,j)-u(i  ,j))
        ! make sure division by 0 does not happen
        if (abs(nom) < 1d-14) then ! nom = 0
          nom = 0d0
          denom = 1d0
        elseif (nom > 1d-14 .and. abs(denom) < 1d-14) then ! nom > 0 => r = \inf
          nom = 1d14
          denom = 1d0
        elseif (nom < -1d-14 .and. abs(denom) < 1d-14) then ! nom < 0 => r = 0
          nom = -1d14
          denom = 1d0
        end if
        ri = nom/denom
        ! calculate r_{i-1} = \frac{u_{i-1}-u_{i-2}}{u_i-u_{i-1}}
        nom   = (u(i-1,j)-u(i-2,j))
        denom = (u(i  ,j)-u(i-1,j))
        ! make sure division by 0 does not happen
        if (abs(nom) < 1d-14) then
          nom = 0d0
          denom = 1d0
        elseif (nom > 1d-14 .and. abs(denom) < 1d-14) then
          nom = 1d14
          denom = 1d0
        elseif (nom < -1d-14.and.abs(denom) < 1d-14) then
          nom = -1d14
          denom = 1d0
        end if
        rim1 = nom/denom
        ! u^l_{i-1/2} = u_{i-1} + 0.5*phi(r_{i-1})*(u_i-u_{i-1})
        ul(j) = u(i-1,j)+0.5d0*phi(rim1)*(u(i  ,j)-u(i-1,j))
        ! u^r_{i-1/2} = u_i + 0.5*phi(r_i)*(u_{i+1}-u_i)
        ur(j) = u(i,j)  -0.5d0*phi(ri  )*(u(i+1,j)-u(i  ,j))
      end do
      ! calculate eigenvalues of \frac{\partial F}{\parial u} at u^l_{i-1/2}
      lam(1) = ev(ul,c0, 1d0)
      lam(2) = ev(ul,c0,-1d0)
      ! calculate eigenvalues of \frac{\partial F}{\parial u} at u^r_{i-1/2}
      lam(3) = ev(ur,c0, 1d0)
      lam(4) = ev(ur,c0,-1d0)
      ! local propagation speed =
      ! max spectral radius (= max eigenvalue of dF/du) of flux Jacobians at (i-1/2)
      a = maxval(abs(lam),dim=1)
      ! calculate pressure via equation of state:
      pr = p_eos(ur(1),c0)
      pl = p_eos(ul(1),c0)
      ! calculate flux based on the Kurganov and Tadmor central scheme
      ! F_{i-1/2} = 0.5*(F(u^r_{i-1/2})+F(u^l_{i-1/2}) - a*(u^r_{i-1/2} - u^l_{i-1/2}))
      ! F_1 = rho * v = u_2
      f(i,1) = 0.5d0*(ur(2)+ul(2)-a*(ur(1)-ul(1)))
      ! F_2 = p + rho * v**2 = p + \frac{u_2**2}{u_1}
      f(i,2) = 0.5d0*(pr+ur(2)**2/ur(1)+pl+ul(2)**2/ul(1)-a*(ur(2)-ul(2)))
    end do

  end subroutine fluxcalc

  real(dp) function phi(r)
    ! Task:
    !   Compute limiter function for fluxes using the Ospre limiter
    ! Input:
    !   r    - ratio of successive gradients
    ! Output:
    !   phi  - limiter

    implicit none

    real(dp), intent(in) :: r

    phi = 0d0
    ! ospre flux limiter phi(r) = \frac{1.5*(r^2+r)}{r^2+r+1}
    if (r > 0d0) then
      phi = 1.5d0*(r**2+r)/(r**2+r+1d0)
    end if

    ! van leer
    !  phi = (r+abs(r))/(1d0+abs(r))
  end function phi

  real(dp) function ev(u,c0,sgn)
    ! Task:
    !   Compute the eigenvalues of the flux jacobian
    ! Input:
    !   u    - state vector at midpoint (left or right approximation)
    !   c0   - speed of sound
    !   sgn  - sign of the root
    ! Output:
    !   ev   - the eigenvalue

    implicit none

    real(dp), dimension(2), intent(inout) :: u
    real(dp), intent(in) :: sgn, c0

    ! calculate root of characteristic equation

    ! \lambda = v \pm c_0 (\frac{\rho}{\rho_0})^{(xi-1)/2}
    ev = u(2)/u(1) + sgn*c0*(u(1)/rho0)**((xi-1d0)/2d0)

    return
  end function ev

  real(dp) function p_eos(rho, c0)
    ! Task:
    !   Compute the pressure via the equation of state from the density
    !   Formula: p = \frac{\rho_0 c_0^2}{\xi} [(\frac{\rho}{\rho_0})^\xi - 1],
    !   where \rho_0 (=1000) is the reference density, c_0 the speed of sound and \xi the
    !   polytropic index (=7 for water).
    ! Input:
    !   rho   - density
    !   c0    - speed of sound
    ! Output:
    !   p     - pressure

    implicit none

    real(dp), intent(in) :: rho, c0

    p_eos = rho0*c0**2d0/xi*((rho/rho0)**xi-1d0)

    return
  end function p_eos

  subroutine bcs(u,p,v,nx,c0)
    ! Task:
    !   Update the velocity and pressure globally from the state vector
    !   Compute the boundary conditions on the ghost cells
    ! Input:
    !   u    - state vector
    !   nx   - number of gridpoints
    !   c0   - speed of sound
    ! Output:
    !   p   - pressure
    !   v   - velocity

    implicit none

    integer, intent(in) :: nx
    real(dp), intent(in) :: c0
    real(dp), dimension(-1:nx+1), intent(inout) :: p,v
    real(dp), dimension(-1:nx+1,2), intent(inout) :: u

    integer :: i

    ! calculate velocity and pressure
    do i=1,nx-1
      ! v = u_2 / u_1
      v(i) = u(i,2)/u(i,1)
      ! p = \frac{rho_0 c0^2}{xi}*((\frac{\rho}{\rho_0})^xi-1), (xi=7, rho_0=1d3)
      p(i) = p_eos(u(i,1),c0)
    end do

    ! calculate boundary conditions
    !  note: for periodic boundary conditions set
    !        f(0) = f(nx-1), f(-1) = f(nx-2), f(nx) = f(1), f(nx+1) = f(2) \forall f
    ! for rho: d\rho/dn = 0 using second order extrapolation
    u(0,1) = (4d0*u(1,1)-u(2,1))/3d0
    u(-1,1) = u(1,1)
    u(nx,1) = (4d0*u(nx-1,1)-u(nx-2,1))/3d0
    u(nx+1,1) = u(nx-1,1)
    ! for p: dp/dn = 0 (thus can use EOS)
    p(0) = p_eos(u(0,1),c0)
    p(-1) = p_eos(u(-1,1),c0)
    p(nx) = p_eos(u(nx,1),c0)
    p(nx+1) = p_eos(u(nx+1,1),c0)
    ! for v: v = 0 using second order extrapolation
    v(0) = 0d0
    v(-1) = v(2)-3d0*v(1)
    v(nx) = 0d0
    v(nx+1) = v(nx-2)-3d0*v(nx-1)
    ! for rho*v
    u(0,2) = u(0,1)*v(0)
    u(-1,2) = u(-1,1)*v(-1)
    u(nx,2) = u(nx,1)*v(nx)
    u(nx+1,2) = u(nx+1,1)*v(nx+1)

  end subroutine bcs

end module mod_fluxcalc
