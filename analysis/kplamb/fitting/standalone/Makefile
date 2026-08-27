ROOT_CONFIG ?= /sciclone/home/jrstevens01/builds/root/root-6.32.08/bin/root-config
AMPTOOLS ?= /sciclone/home/jrstevens01/wm_gluex/third_party/AmpTools-0.11.0

CXX ?= g++
CXXFLAGS ?= -O3 $(shell $(ROOT_CONFIG) --cflags) -std=c++17
CXXFLAGS += -D_GLIBCXX_USE_CXX11_ABI=0

INC_DIR := -I. -I$(shell $(ROOT_CONFIG) --incdir) -I$(AMPTOOLS)
LIB_DIR := -L$(shell $(ROOT_CONFIG) --libdir) -L$(AMPTOOLS)/lib
ROOTLIBS := $(shell $(ROOT_CONFIG) --libs)
LIBS := -lAmpTools -lIUAmpTools -lMinuitInterface -lUpRootMinuit $(ROOTLIBS) -lstdc++

TARGETS := fit khyperon_plotter

.PHONY: all clean

all: $(TARGETS)

fit: fit.o FSRootDataReader.o KHyperon.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(INC_DIR) $(LIB_DIR) $(LIBS)

khyperon_plotter: khyperon_plotter.o FSRootDataReader.o KHyperon.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(INC_DIR) $(LIB_DIR) $(LIBS)

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< $(INC_DIR) -o $@

clean:
	rm -f *.o fit khyperon_plotter
