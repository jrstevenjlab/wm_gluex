# Build the standalone KHyperon fit, plotter and generator.
#
# Requirements: ROOT (root-config on PATH), a C++17 compiler, make and git.
#
#   make amptools   # one time: download and build AmpTools into ./AmpTools
#   make            # build fit, khyperon_plotter and gen_khyperon
#
# Variables you can override on the command line or in the environment:
#   ROOT_CONFIG  root-config executable             [root-config]
#   AMPTOOLS     AmpTools source tree (the directory with IUAmpTools/ and lib/)
#                                                    [./AmpTools/AmpTools]
#   CXXSTD       C++ standard                        [c++17]
#   OLD_ABI=1    add -D_GLIBCXX_USE_CXX11_ABI=0, only needed when linking against
#                an AmpTools build made with the old libstdc++ ABI

ROOT_CONFIG ?= root-config
AMPTOOLS ?= $(CURDIR)/AmpTools/AmpTools
AMPTOOLS_VERSION ?= v0.15.3
AMPTOOLS_REPO ?= https://github.com/mashephe/AmpTools.git
CXXSTD ?= c++17

ROOT_CFLAGS := $(shell $(ROOT_CONFIG) --cflags 2>/dev/null)

# -std is added after the ROOT flags so it takes precedence
CXXFLAGS ?= -O3 $(ROOT_CFLAGS) -std=$(CXXSTD)
ifdef OLD_ABI
CXXFLAGS += -D_GLIBCXX_USE_CXX11_ABI=0
endif

INC_DIR := -I. -I$(shell $(ROOT_CONFIG) --incdir 2>/dev/null) -I$(AMPTOOLS)
LIB_DIR := -L$(shell $(ROOT_CONFIG) --libdir 2>/dev/null) -L$(AMPTOOLS)/lib
ROOTLIBS := $(shell $(ROOT_CONFIG) --libs 2>/dev/null)

# AmpTools 0.11 builds separate component archives; newer releases only
# build the combined libAmpTools.a
ifneq ($(wildcard $(AMPTOOLS)/lib/libIUAmpTools.a),)
AMPTOOLSLIBS := -lAmpTools -lIUAmpTools -lMinuitInterface -lUpRootMinuit
else
AMPTOOLSLIBS := -lAmpTools
endif
LIBS := $(AMPTOOLSLIBS) $(ROOTLIBS) -lstdc++

TARGETS := fit khyperon_plotter gen_khyperon

.PHONY: all clean amptools clean-amptools check-deps

all: $(TARGETS)

fit: fit.o FSRootDataReader.o KHyperon.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(INC_DIR) $(LIB_DIR) $(LIBS)

khyperon_plotter: khyperon_plotter.o FSRootDataReader.o KHyperon.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(INC_DIR) $(LIB_DIR) $(LIBS)

gen_khyperon: gen_khyperon.o KHyperon.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(INC_DIR) $(LIB_DIR) $(LIBS)

%.o: %.cc | check-deps
	$(CXX) $(CXXFLAGS) -c $< $(INC_DIR) -o $@

check-deps:
	@command -v $(ROOT_CONFIG) >/dev/null 2>&1 || { \
	  echo "ERROR: '$(ROOT_CONFIG)' not found.  Install ROOT and put root-config on your PATH"; \
	  echo "       (for example: source <root install>/bin/thisroot.sh), see README.md."; exit 1; }
	@test -f "$(AMPTOOLS)/lib/libAmpTools.a" || { \
	  echo "ERROR: no AmpTools build found in $(AMPTOOLS)."; \
	  echo "       Run 'make amptools' first, or set AMPTOOLS=/path/to/AmpTools/AmpTools."; exit 1; }

# Download and build AmpTools.  The ROOT flags are passed explicitly with
# -std=$(CXXSTD) because recent ROOT releases report -std=c++20, which the
# AmpTools sources do not compile with on macOS.
amptools:
	@command -v $(ROOT_CONFIG) >/dev/null 2>&1 || { \
	  echo "ERROR: '$(ROOT_CONFIG)' not found.  Install ROOT first, see README.md."; exit 1; }
	@if [ ! -d AmpTools ]; then \
	  git clone --depth 1 --branch $(AMPTOOLS_VERSION) $(AMPTOOLS_REPO) AmpTools; \
	fi
	$(MAKE) -C AmpTools/AmpTools AMPTOOLS_HOME=$(CURDIR)/AmpTools CXX="$(CXX)" \
	  CXX_FLAGS='-O3 -DAMPTOOLS_GDOUBLE_FP64 -DATVERSION=\"$(AMPTOOLS_VERSION)\" $(ROOT_CFLAGS) -std=$(CXXSTD)'

clean:
	rm -f *.o $(TARGETS)

clean-amptools:
	rm -rf AmpTools
