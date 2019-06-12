all: cell-test

TRANSFINITE=/home/salvi/project/transfinite/

CXXFLAGS=-Wall -g -pedantic -std=c++17 -I$(TRANSFINITE)/src/geom -I$(TRANSFINITE)/src/transfinite
LIBS=-L$(TRANSFINITE)/release/geom -L$(TRANSFINITE)/release/transfinite -lgeom -ltransfinite

cell-test: cell-test.o gb-io.o
	g++ -o $@ $^ $(LIBS)
