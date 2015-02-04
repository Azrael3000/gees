!=====================================================
!
! GEES - GPL Euler equation solver
!
!=====================================================
!
! Written 2012 by Arno Mayrhofer (www.amconception.de)
!
!=====================================================
!
! Copyright 2012 Arno Mayrhofer
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

  use mod_fluxcalc

  implicit none

  real(dp), dimension(:),allocatable :: p,v
  real(dp), dimension(:,:), allocatable :: u,f, utild, uold
  real(dp) :: temps, dt, tend, dx, odt, pi, c0
  integer :: i, nt, it, nx, io, lout
  character(len=13) :: fname

  write(*,*) 'Welcome to GEES'
  write(*,*) 'your friendly GPL Euler Equation solver'
  write(*,*) 'written 2012 by Arno Mayrhofer (www.amconception.de)'
  write(*,*)
  write(*,*) 'Number of grid points:'
  read(*,*) nx
  write(*,*) nx
  write(*,*) 'End time:'
  read(*,*) tend
  write(*,*) tend
  write(*,*) 'Time-step size'
  read(*,*) dt
  write(*,*) dt
  write(*,*) 'Output dt:'
  read(*,*) odt
  write(*,*) odt
  write(*,*) 'Speed of sound:'
  read(*,*) c0
  write(*,*) c0

  ! number of time-steps
  nt = int(tend/dt+1d-14)
  ! grid size
  dx = 1D0/real(nx,8)
  pi = acos(-1d0)
  ! list ouput index
  lout = 1

  allocate(p(-1:nx+1),v(-1:nx+1),f(nx,2),u(-1:nx+1,2), &
  utild(-1:nx+1,2), uold(-1:nx+1,2))

  ! init
  do i=1,nx-1
    u(i,1) = 1d3+cos(2*pi*real(i,8)/real(nx,8))*10d0 !rho
    u(i,2) = 0d0 !u(i,1)*0.1d0*sin(2*pi*real(i,8)/real(nx,8)) !rho*v
  end do
  ! set p,v from u and calcuate boundary values at bdry and ghost cells
  call bcs(u,p,v,nx,c0)
  ! file output index
  io = 0
  ! simulation time
  temps = 0d0

  ! output
  write(fname,'(i5.5,a8)') io,'.out.csv'
  open(file=fname,unit=800)
  write(800,*) 'xpos,p,rho,v'
  do i=-1,nx+1
    write(800,'(4(a1,e20.10))') ' ', real(i,8)/real(nx,8),',', p(i),',', u(i,1),',', v(i)
  end do
  close(800)
  io = io+1

  ! loop over all timesteps
  do it=1,nt
    ! list output
    if (real(nt*lout)*0.1d0 < real(it)) then
      write(*,*) 'Calculated ', int(real(lout)*10.), '%'
      lout = lout + 1
    end if
    temps = temps + dt
    uold = u

    ! First Runge-Kutta step
    ! calculate flux at mid points using u
    call fluxcalc(u,v,f,nx,c0)
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
    call fluxcalc(utild,v,f,nx,c0)
    do i=1,nx-1
      ! calc k2 = -dt/dx(f(utild,i+1/2) - f(utild,i-1/2))
      utild(i,:) = -dt/dx*(f(i+1,:)-f(i,:))
      ! u^n = u^n + 1/6 k2
      u(i,:) = u(i,:) + utild(i,:)/6d0
      ! utild = uold^n + 1/2 k_2
      utild(i,:) = uold(i,:) + utild(i,:)*0.5d0
    end do
    ! calculate p,v + bcs for uold + 1/2 k2 = utild
    call bcs(utild,p,v,nx,c0)

    ! Third Runge-Kutta step
    ! calculate flux at mid points using utild
    call fluxcalc(utild,v,f,nx,c0)
    do i=1,nx-1
      ! calc k3 = -dt/dx(f(utild,i+1/2) - f(utild,i-1/2))
      utild(i,:) = -dt/dx*(f(i+1,:)-f(i,:))
      ! u^n = u^n + 1/6 k3
      u(i,:) = u(i,:) + utild(i,:)/6d0
      ! utild = uold^n + k_3
      utild(i,:) = uold(i,:) + utild(i,:)
    end do
    ! calculate p,v + bcs for uold + k3 = utild
    call bcs(utild,p,v,nx,c0)

    ! Fourth Runge-Kutta step
    ! calculate flux at mid points using utild
    call fluxcalc(utild,v,f,nx,c0)
    do i=1,nx-1
      ! calc k4 = -dt/dx(f(utild,i+1/2) - f(utild,i-1/2))
      utild(i,:) = -dt/dx*(f(i+1,:)-f(i,:))
      ! u^n = u^n + 1/6 k4
      u(i,:) = u(i,:) + utild(i,:)/6d0
    end do
    ! calculate p,v + bcs for (uold + k1 + k2 + k3 + k4)/6 = u
    call bcs(u,p,v,nx,c0)

    ! output
    if (abs(temps-odt*real(io,8)) < abs(temps+dt-odt*real(io,8)) .or. it == nt) then
      write(fname,'(i5.5,a8)') io,'.out.csv'
      open(file=fname,unit=800)
      write(800,*) 'xpos,p,rho,v'
      do i=-1,nx+1
        write(800,'(4(a1,e20.10))') ' ', real(i,8)/real(nx,8),',', p(i),',', u(i,1),',', v(i)
      end do
      close(800)
      io = io + 1
    end if
  end do

  write(*,*) 'Calculation finished. Goodbye.'

end program gees
