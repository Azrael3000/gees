#=====================================================
#
# Copyright 2012 Arno Mayrhofer
# This program is licensed under the GNU General Public License.
#
# This file is part of Gees.
#
# Gees is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Gees is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Gees.  If not, see <http://www.gnu.org/licenses/>.
#
#=====================================================

FFLAGS = -Wall -fbacktrace -g -pg -O3
TARGET = gees
all: $(TARGET)

$(TARGET): fluxcalc.o $(TARGET).o
	gfortran $(FFLAGS) -o $(TARGET) fluxcalc.o $(TARGET).o

%.o: %.f90
	  gfortran $(FFLAGS) -c $< -o $@

clean:
	rm -f *.o *.mod $(TARGET)
