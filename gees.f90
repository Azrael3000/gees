!=====================================================
!
! GEES - GPL Euler Equation Solver
!
!=====================================================
!
! Written 2012-2016 by Arno Mayrhofer
! Homepage: www.amconception.de
! Github: https://github.com/Azrael3000/gees
!
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
!
! Implementation details:
! Time-stepping: 4th order explicit Runge Kutta
! Flux calculation: MUSCL scheme (Kurganov and Tadmor central scheme)
! Flux limiter: Ospre
! Boundary conditions: Second order polynomial extrapolation
! Equation of state: Tait
!
!=====================================================

program gees

  use mod_helper
  use mod_fluxcalc

  implicit none

  ! output buffers
  !   p - pressure
  !   v - velocity
  real(dp), dimension(:), allocatable :: p,v
  ! work buffers
  !   u     - state vector
  !   f     - flux vector
  !   utild - temporary state vector
  !   uold  - state vector from previous time step
  real(dp), dimension(:,:), allocatable :: u,f, utild, uold
  ! floating point variables
  !   time  - simulation time
  !   dt    - time step size
  !   tend  - maximum simulation time
  !   dx    - grid spacing
  !   odt   - output time step size
  !   c0    - speed of sound
  real(dp) :: time, dt, tend, dx, odt, c0
  ! constants
  !   pi    - 3.1415927...12409435094120938202... or so
  real(dp), parameter :: pi = acos(-1d0)
  ! integer variables
  !   i     - iterator variable
  !   nt    - number of time steps
  !   it    - current time step
  !   nx    - number of gridpoints
  !   io    - current output index
  !   lout  - current status output index
  integer :: i, nt, it, nx, io, lout

  call get_input(nx, tend, dt, odt, c0)

  call initialize_simulation(tend, dt, nx, nt, dx, lout, io, time, p, v, f, u, utild, uold)

  ! init
  do i=1,nx-1
    u(i,1) = 1d3+cos(2*pi*real(i,8)/real(nx,8))*10d0 !rho
    u(i,2) = 0d0 !u(i,1)*0.1d0*sin(2*pi*real(i,8)/real(nx,8)) !rho*v
  end do
  ! set p,v from u and calcuate boundary values at bdry and ghost cells
  call bcs(u,p,v,nx,c0)

  ! output
  call write_output(time, odt, io, dt, it, nt, nx, p, u, v)

  ! loop over all timesteps
  do it=1,nt

    call show_progress(nt, it, lout)

    time = time + dt
    uold = u

    ! First Runge-Kutta step
    ! calculate flux at mid points using u
    call fluxcalc(u,f,nx,c0)
    do i=1,nx-1
      ! calc k1 = -dt/dx(f(u,i+1/2) - f(u,i-1/2))
      utild(i,:) = -dt/dx*(f(i+1,:)-f(i,:))
      ! u^n = uold^n + 1/6 k1
      u(i,:) = uold(i,:) + utild(i,:)/6d0
      ! utild = uold^n + 1/2 k_1
      utild(i,:) = uold(i,:) + utild(i,:)*0.5d0
    end do
    ! calculate p,v + bcs for uold + 1/2 k1 = utild
    call bcs(utild,p,v,nx,c0)

    ! Second Runge-Kutta step
    ! calculate flux at mid points using utild
    call fluxcalc(utild,f,nx,c0)
    do i=1,nx-1
      ! calc k2 = -dt/dx(f(utild,i+1/2) - f(utild,i-1/2))
      utild(i,:) = -dt/dx*(f(i+1,:)-f(i,:))
      ! u^n = u^n + 1/3 k2
      u(i,:) = u(i,:) + utild(i,:)/3d0
      ! utild = uold^n + 1/2 k_2
      utild(i,:) = uold(i,:) + utild(i,:)*0.5d0
    end do
    ! calculate p,v + bcs for uold + 1/2 k2 = utild
    call bcs(utild,p,v,nx,c0)

    ! Third Runge-Kutta step
    ! calculate flux at mid points using utild
    call fluxcalc(utild,f,nx,c0)
    do i=1,nx-1
      ! calc k3 = -dt/dx(f(utild,i+1/2) - f(utild,i-1/2))
      utild(i,:) = -dt/dx*(f(i+1,:)-f(i,:))
      ! u^n = u^n + 1/3 k3
      u(i,:) = u(i,:) + utild(i,:)/3d0
      ! utild = uold^n + k_3
      utild(i,:) = uold(i,:) + utild(i,:)
    end do
    ! calculate p,v + bcs for uold + k3 = utild
    call bcs(utild,p,v,nx,c0)

    ! Fourth Runge-Kutta step
    ! calculate flux at mid points using utild
    call fluxcalc(utild,f,nx,c0)
    do i=1,nx-1
      ! calc k4 = -dt/dx(f(utild,i+1/2) - f(utild,i-1/2))
      utild(i,:) = -dt/dx*(f(i+1,:)-f(i,:))
      ! u^n = u^n + 1/6 k4
      u(i,:) = u(i,:) + utild(i,:)/6d0
    end do
    ! calculate p,v + bcs for uold + dt*(k1 + 2*k2 + 2*k3 + k4)/6 = u
    call bcs(u,p,v,nx,c0)

    ! output
    call write_output(time, odt, io, dt, it, nt, nx, p, u, v)
  end do

  write(*,*) 'Calculation finished. Goodbye.'

end program gees
