CXX = g++

ROOT_DIR=/home/zhijie/root/root
ROOT_CFLAGS=$(shell $(ROOT_DIR)/bin/root-config --cflags)
ROOT_LDFLAGS=$(shell $(ROOT_DIR)/bin/root-config --ldflags)
ROOT_LIBS=$(shell $(ROOT_DIR)/bin/root-config --libs)

HEPMC2_DIR=/home/zhijie/hepmc2/install
HEPMC2_INC=$(HEPMC2_DIR)/include
HEPMC2_LIB=$(HEPMC2_DIR)/lib

FASTJET_DIR=/home/zhijie/fastjet3/install
FASTJET_CFLAGS=$(shell $(FASTJET_DIR)/bin/fastjet-config --cxxflags)
FASTJET_LDFLAGS=$(shell $(FASTJET_DIR)/bin/fastjet-config --libs)

CFLAGS=-I$(HEPMC2_INC) $(FASTJET_CFLAGS) $(ROOT_CFLAGS) -Wall
LDFLAGS=-L$(HEPMC2_LIB) -Wl,-rpath $(HEPMC2_LIB) $(FASTJET_LDFLAGS) $(ROOT_LDFALGS)
LIBS=-lHepMC $(ROOT_LIBS)

SRC=main.cpp Isolation.cpp Pdg_check.cpp
OBJ=main.o Isolation.o Pdg_check.o

all: $(SRC) $(OBJ) main

main: $(OBJ)
	$(CXX) $^ -o $@ $(LDFLAGS) $(LIBS)

$(OBJ): $(SRC)
	$(CXX) $(CFLAGS) -c $^

clean:
	rm *.o main
