FC = gfortran
FFLAGS = -O2 -ffixed-form -Wall -mcmodel=large

OBJS = fmcp_data.o \
       fmcp_read.o \
       fmcp_interactions.o \
       fmcp_targets.o \
       fmcp_flexibility.o \
       fmcp_dijkstra.o \
       fmcp_vmd.o \
       main.o

.PHONY: all clean

all: fmcp

fmcp: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

fmcp_data.o: fmcp_data.f
	$(FC) $(FFLAGS) -c $<

fmcp_read.o: fmcp_read.f fmcp_data.o
	$(FC) $(FFLAGS) -c $<

fmcp_interactions.o: fmcp_interactions.f fmcp_data.o
	$(FC) $(FFLAGS) -c $<

fmcp_targets.o: fmcp_targets.f fmcp_data.o
	$(FC) $(FFLAGS) -c $<

fmcp_flexibility.o: fmcp_flexibility.f fmcp_data.o
	$(FC) $(FFLAGS) -c $<

fmcp_dijkstra.o: fmcp_dijkstra.f fmcp_data.o
	$(FC) $(FFLAGS) -c $<

fmcp_vmd.o: fmcp_vmd.f fmcp_data.o
	$(FC) $(FFLAGS) -c $<

main.o: main.f fmcp_data.o fmcp_read.o fmcp_interactions.o fmcp_targets.o fmcp_flexibility.o fmcp_dijkstra.o fmcp_vmd.o
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.mod fmcp
