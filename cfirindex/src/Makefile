SRC_ROOT = $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

CXXFLAGS +=-DGCC

INCPATHS = -I../include

OBJS = dslim_query.o dslim.o lbe_internal.o utils.o mods.o msquery.o dslim_fileout.o scheduler.o dslim_comm.o dslim_score.o dslim_scproc.o
LIB = libdslim.a
OBJDIR = objs

.PHONY: all lbe clean

all: lbe $(OBJS)
	ar rcs $(LIB) $(OBJDIR)/*.o

lbe:
	mkdir -p $(OBJDIR)

%.o:	$(SRC)%.cpp

	$(CXX) -c $(CFLAGS) $(CXXFLAGS) $(CPPFLAGS) $(INCPATHS) -o $(OBJDIR)/$@ $<

clean:
	rm -rf $(OBJDIR) $(LIB)
