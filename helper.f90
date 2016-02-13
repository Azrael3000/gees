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


end module mod_helper
