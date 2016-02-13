#=====================================================
#
# Copyright 2012-2016 Arno Mayrhofer
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

# Compiler flags
FFLAGS = -Wall -fbacktrace -g -pg -O3
# Name of the executable
TARGET = gees
# Compiler
COMPILER = gfortran
# All *.f90 files
F90FILES = $(wildcard *.f90)
# All object files created by compiling the *.f90 files
OFILES = $(patsubst %.f90, %.o, $(F90FILES))

# Make sure that clean is run every time it's called
.PHONY: clean

# Main make target is the executable
all: $(TARGET)

# To build the executable we require all object files
$(TARGET): $(OFILES)
# Linking the object files to create the target executable
	$(COMPILER) $(FFLAGS) -o $(TARGET) $(OFILES)

# Recipe how to build object files, which require the *.f90 files
%.o: %.f90
# Compiling one f90 file ($<) to create the corresponding .o file ($@)
	$(COMPILER) $(FFLAGS) -c $< -o $@

# Cleaning the workspace
clean:
# Remove all object and *.mod files as well as the executable and the performance file
	rm -f $(OFILES) *.mod $(TARGET) gmon.out
