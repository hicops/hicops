BUILD_MODE = RELEASE
OPENMP = YES
JOBS=1
CXXFLAGS = -Wall -fmessage-length=0 -std=gnu++0x
LNKFLAGS = -Wl,--wrap,memcpy

# Needs gperftools to be installed and configured before use
PROFILE = NO

PROJECT_ROOT = $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
SRC = $(PROJECT_ROOT)/src

DEPS = src/lbe_internal.cpp src/utils.cpp src/mods.cpp  include/common.h include/config.h include/utils.h include/mods.h include/slm_dsts.h include/slmerr.h include/lbe.h

LIBS = -ldslim
LIBPATHS = -L$(SRC)
EXECUTEABLE = counter.exe

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

all: lbe
	$(CXX) $(LNKFLAGS) $(LIBPATHS) -Wl,--start-group  $(LIBS) ./slm.o -Wl,--end-group -o $(EXECUTEABLE) -Wl,-Map=$(EXECUTEABLE).map

lbe: $(DEPS)
	$(MAKE) -j$(JOBS) -C $(SRC)
	$(CXX) -c -I./include -I./mstoolkit/include $(CXXFLAGS) -DGCC ./slm.cpp -o slm.o

clean:
	$(MAKE) -C $(SRC) clean
	rm -rf $(EXECUTEABLE) $(EXECUTEABLE).map slm.o

allclean:
	$(MAKE) -C $(SRC) clean
	rm -rf $(EXECUTEABLE) $(EXECUTEABLE).map slm.o

.PHONY: all lbe clean allclean
