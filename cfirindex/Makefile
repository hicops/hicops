BUILD_MODE = RELEASE
OPENMP = YES
MPI = YES
JOBS=22
CXXFLAGS = -Wall -fmessage-length=0 -std=c++14
#LNKFLAGS = -Wl,--wrap,memcpy

# Change CXX to MPIC++ when MPI enabled
ifeq ($(MPI),YES)
	ifneq ($(OS),Windows_NT)
		CXX:= mpicxx
		CXXFLAGS += -DDISTMEM
	endif
endif

# Needs gperftools to be installed and configured before use
PROFILE = NO

PROJECT_ROOT = $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
SRC = $(PROJECT_ROOT)/src

DEPS = src/dslim_query.cpp src/dslim_comm.cpp src/scheduler.cpp src/dslim.cpp src/dslim_fileout.cpp src/lbe_internal.cpp src/utils.cpp src/mods.cpp src/msquery.cpp src/dslim_score.cpp src/dslim_scproc.cpp src/expeRT.cpp src/sgsmooth.cpp include/common.h include/config.h include/dslim_fileout.h include/keyval.h include/utils.h include/mods.h include/msquery.h include/slm_dsts.h include/dslim.h include/slmerr.h include/lbe.h include/sgsmooth.h include/expeRT.h include/lwvector.h include/dslim_score.h include/dslim_comm.h include/scheduler.h include/lwbuff.h include/lwqueue.h include/minheap.h

LIBS = -ldslim
LIBPATHS = -L$(SRC)
EXECUTEABLE = hicops.exe

ifeq ($(PROFILE),YES)
	LIBS += -lprofiler
	CXXFLAGS += -D_PROFILE -g
endif

ifeq ($(OPENMP),YES)
	CXXFLAGS += -fopenmp -D_GLIBCXX_PARALLEL
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
export CXX

all: hicops
	$(CXX) $(LNKFLAGS) $(LIBPATHS) -Wl,--start-group  $(LIBS) ./hicops.o -Wl,--end-group -o $(EXECUTEABLE) -Wl,-Map=$(EXECUTEABLE).map

hicops: $(DEPS)
	$(MAKE) -j$(JOBS) -C $(SRC)
	$(CXX) -c -I./include $(CXXFLAGS) -DGCC ./hicops.cpp -o hicops.o

clean:
	$(MAKE) -C $(SRC) clean
	rm -rf $(EXECUTEABLE) $(EXECUTEABLE).map hicops.o
	
distclean: clean

allclean: clean

.PHONY: all hicops clean distclean allclean
