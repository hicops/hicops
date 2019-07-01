BUILD_MODE = RELEASE
OPENMP = YES
JOBS=1
CXXFLAGS = -Wall -fmessage-length=0 -std=c++11
LNKFLAGS = 

# Needs gperftools to be installed and configured before use
PROFILE = NO

PROJECT_ROOT = $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
SRC = $(PROJECT_ROOT)/src
MSTOOLKIT = $(PROJECT_ROOT)/mstoolkit

DEPS = src/dslim_query.cpp src/dslim.cpp src/hyperscore.cpp src/lbe_internal.cpp src/utils.cpp src/mods.cpp src/msquery.cpp include/common.h include/config.h include/hyperscore.h include/keyval.h include/utils.h include/mods.h include/msquery.h include/slm_dsts.h include/dslim.h include/slmerr.h include/lbe.h

LIBS = -ldslim -lmstoolkitlite 
LIBPATHS = -L$(SRC) -L$(MSTOOLKIT)
EXECUTEABLE = cfir.exe

ifeq ($(PROFILE),YES)
	LIBS += -lprofiler
	CXXFLAGS += -D_PROFILE -g
endif

ifeq ($(OPENMP),YES)
	CXXFLAGS += -fopenmp
	LNKFLAGS += -fopenmp
endif

# Setup build mode
ifeq ($(BUILD_MODE),DEBUG)
	CXXFLAGS += -g3 -O0 
else ifeq ($(BUILD_MODE),RELEASE)
	CXXFLAGS += -O3
else
	$(error Build mode $(BUILD_MODE) not supported by this Makefile)
endif

# Setup parallel build jobs
ifeq ($(OS),Windows_NT)
	JOBS := $(NUMBER_OF_PROCESSORS)
else
	JOBS := $(shell grep -c ^processor /proc/cpuinfo)
endif

# Export Build Variables
export BUILD_MODE
export OPENMP
export JOBS
export CXXFLAGS

all: mstoolkit lbe
	$(CXX) $(LNKFLAGS) $(LIBPATHS) -Wl,--start-group  $(LIBS) ./slm.o -Wl,--end-group -o $(EXECUTEABLE) -Wl,-Map=$(EXECUTEABLE).map

mstoolkit:
	$(MAKE) -j$(JOBS) -C $(MSTOOLKIT)

lbe: $(DEPS)
	$(MAKE) -j$(JOBS) -C $(SRC)
	$(CXX) -c -I./include -I./mstoolkit/include $(CXXFLAGS) -D_GLIBCXX_PARALLEL -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DHAVE_EXPAT_CONFIG_H -DGCC -D_NOTHERMORAW -D_NOSQLITE ./slm.cpp -o slm.o

clean:
	$(MAKE) -C $(SRC) clean
	rm -rf $(EXECUTEABLE) $(EXECUTEABLE).map slm.o

allclean:
	$(MAKE) -C $(SRC) clean
	rm -rf $(EXECUTEABLE) $(EXECUTEABLE).map slm.o
	$(MAKE) -C $(MSTOOLKIT) clean

.PHONY: all mstoolkit lbe clean allclean
