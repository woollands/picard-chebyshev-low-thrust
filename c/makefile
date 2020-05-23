#
#  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
#  DATE WRITTEN:     Mar 2020
#  LAST MODIFIED:    Mar 2020
#  AFFILIATION:      Jet Propulsion Laboratory, California Institute of Technology, Pasadena, CA
#  DESCRIPTION:      File to compile Bulirsch-Stoer numerical integrator
#

CC = gcc
CPP = g++

CFLAGS = -I$(ODIR) -I$(IDIR)

ODIR = obj
CSPICE = ../../../cspice/lib
IDIR = include
SDIR = src

LIBS = -lm -lnlopt

_DEPS1_inc = const.hpp states_twobody_thrust.hpp spacecraft_params.hpp eclipse_model.hpp
DEPS1_inc = $(patsubst %,$(IDIR)/%,$(_DEPS1_inc))

_SPICE = cspice.a
SPICE_lib = $(patsubst %,$(CSPICE)/%,$(_SPICE))

_OBJ1 = main_simulation.o states_twobody_thrust.o eclipse_model.o
OBJ1 = $(patsubst %,$(ODIR)/%,$(_OBJ1))

$(ODIR)/%.o: $(SDIR)/%.cpp
	$(CPP) -c -o $@ $< $(CFLAGS)

$(ODIR)/%.o: %.c $(DEPS1_inc)
	$(CC) -c -o $@ $< $(CFLAGS)

main_simulation: $(OBJ1) $(SPICE_lib)
	$(CPP) -g -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~
