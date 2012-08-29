#    Makefile for Icecamp L1L2 Model
#    Copyright (C) Joshua Charles Campbell 2012

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

FC=gfortran
FCFLAGS=-Wall -Wextra -g -fbounds-check -fbacktrace -ffree-line-length-none
LDFLAGS=-lblas
PROGRAMS=l1l2 fem
SLAPS=dbcg.o dcgn.o dgmres.o dlaputil.o dmvops.o mach.o dcg.o dcgs.o dir.o dmset.o domn.o xersla.o mrgrnk.o

default: test

%.o: %.f
	$(FC) $(FCFLAGS) -c $<

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

$(PROGRAMS): %: %.F90 $(SLAPS)
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

test: l1l2
	time ./l1l2
	./plot.py <plot.txt

clean:
	rm -f *.o *.mod *.MOD $(PROGRAMS)
