BUILD_MODE = RELEASE
OPENMP = YES
JOBS=1
CXXFLAGS = -Wall -fmessage-length=0 -std=c++11
LNKFLAGS = 

PROJECT_ROOT = $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
SRC = $(PROJECT_ROOT)/src

DEPS = src/slm.cpp src/lbe_internal.cpp src/utils.cpp src/mods.cpp include/common.h include/config.h include/keyval.h include/utils.h include/mods.h include/slm_dsts.h include/slmerr.h include/lbe.h

LIBS = -ldslim 
LIBPATHS = -L$(SRC)
EXECUTEABLE = seqgen.exe

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

all: lbe
	$(CXX) $(LNKFLAGS) $(LIBPATHS) -Wl,--start-group $(LIBS) -Wl,--end-group -o $(EXECUTEABLE) -Wl,-Map=$(EXECUTEABLE).map

lbe: $(DEPS)
	$(MAKE) -j$(JOBS) -C $(SRC)

clean:
	$(MAKE) -C $(SRC) clean
	rm -rf $(EXECUTEABLE) $(EXECUTEABLE).map

.PHONY: all lbe clean