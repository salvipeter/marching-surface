all: cell-test fixed-volume

TRANSFINITE=/home/salvi/project/transfinite/
INCLUDES=-I/usr/include/eigen3 -I$(TRANSFINITE)/src/geom -I$(TRANSFINITE)/src/transfinite

CXXFLAGS=-Wall -g -pedantic -std=c++17 $(INCLUDES)
LIBS=-L$(TRANSFINITE)/release/geom -L$(TRANSFINITE)/release/transfinite -lgeom -ltransfinite

cell-test: cell-test.o gb-io.o superd-io.o
	g++ -o $@ $^ $(LIBS)

fixed-volume: fixed-volume.o gb-io.o
	g++ -o $@ $^ $(LIBS)
