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

module mod_helper

  ! dp is for the kind definition of the reals, ensuring double precision
  integer, parameter :: dp = selected_real_kind(15, 307)

contains

  subroutine get_input(nx, tend, dt, odt, c0)
    ! Task:
    !   Read input from command line
    ! Output:
    !   nx    - number of grid points
    !   tend  - end time of simulation
    !   dt    - time-step size
    !   odt   - output time-step size
    !   c0    - speed of sound

    implicit none

    integer, intent(out) :: nx
    real(dp), intent(out) :: tend, dt, odt, c0

    write(*,*) 'Welcome to GEES'
    write(*,*) 'your friendly GPL Euler Equation Solver'
    write(*,*) 'written 2012-2016 by Arno Mayrhofer (www.amconception.de)'
    write(*,*)
    write(*,*) 'For the newest version, help and to report bugs please visit the'
    write(*,*) 'Github repository: https://github.com/Azrael3000/gees'
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
    write(*,*)

  end subroutine get_input

  subroutine show_progress(nt, it, lout)
    ! Task:
    !   Write percentage of completion (in 10% steps) to the command line (stdout)
    ! Input:
    !   nt    - Total number of iterations
    !   it    - Current iteration
    !   lout  - Current status output index (from 0 to 10)
    ! Output:
    !   lout  - Next status output index

    implicit none

    integer, intent(in) :: nt, it
    integer, intent(inout) :: lout

    ! list output
    if (real(nt*lout)*0.1d0 < real(it)) then
      write(*,*) 'Calculated ', int(real(lout)*10.), '%'
      lout = lout + 1
    end if

  end subroutine show_progress

  subroutine write_output(time, odt, io, dt, it, nt, nx, p, u, v)
    ! Task:
    !   Write data to file
    ! Input:
    !   time  - Total number of iterations
    !   odt   - Output time-step size
    !   io    - Current output index
    !   dt    - Integrator time-step size
    !   it    - Current iteration
    !   nt    - Number of total iterations
    !   nx    - Number of grid points
    !   p     - Pressure
    !   u     - State vector
    !   v     - Velocity
    ! Output:
    !   io    - Next output index

    implicit none

    integer, intent(in) :: it, nt, nx
    real(dp), intent(in) :: time, odt, dt
    real(dp), dimension(-1:nx+1,2), intent(in) :: u
    real(dp), dimension(-1:nx+1), intent(in) :: p,v
    integer, intent(inout) :: io

    ! file name for output
    character(len=13) :: fname
    ! iterator
    integer :: i

    if (abs(time-odt*real(io,8)) < abs(time+dt-odt*real(io,8))  &
        .or. it == nt                                           &
        .or. io == 0                                            &
       ) then
      write(fname,'(i5.5,a8)') io,'.out.csv'
      open(file=fname,unit=800)
      write(800,*) 'xpos,p,rho,v'
      do i=-1,nx+1
        write(800,'(4(a1,e20.10))') ' ', real(i,8)/real(nx,8),',', p(i),',', u(i,1),',', v(i)
      end do
      close(800)
      io = io + 1
    end if

  end subroutine write_output

  subroutine initialize_simulation(tend, dt, nx, nt, dx, lout, io, time, p, v, f, u, utild, uold)
    ! Task:
    !   Write data to file
    ! Input:
    !   tend  - Total simulation duration
    !   dt    - Integrator time-step size
    !   nx    - Number of grid points
    ! Output:
    !   nt    - Number of total iterations
    !   dx    - Grid size
    !   lout  - Current status output index (from 0 to 10)
    !   io    - Current output index
    !   time  - Current simulation time
    !   p     - pressure
    !   v     - velocity
    !   f     - fluxes
    !   u     - state array
    !   utild - temporary state array
    !   uold  - state array of previous time-step

    implicit none

    integer, intent(in) :: nx
    real(dp), intent(in) :: tend, dt
    integer, intent(out) :: nt, lout, io
    real(dp), intent(out) :: dx, time
    real(dp), dimension(:,:), allocatable, intent(out) :: u, f, utild, uold
    real(dp), dimension(:), allocatable, intent(out) :: p, v

    ! number of time-steps
    nt = int(tend/dt+1d-14)
    ! grid size
    dx = 1d0/real(nx,8)
    ! list ouput index
    lout = 0
    ! file output index
    io = 0
    ! simulation time
    time = 0d0

    ! Allocate arrays for pressure, velocity, flux and state vectors
    allocate( p(-1:nx+1),       &
              v(-1:nx+1),       &
              f(nx,2),          &
              u(-1:nx+1,2),     &
              utild(-1:nx+1,2), &
              uold(-1:nx+1,2))

  end subroutine initialize_simulation

end module mod_helper
