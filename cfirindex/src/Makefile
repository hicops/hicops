SRC_ROOT = $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

CXXFLAGS += -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DHAVE_EXPAT_CONFIG_H -DGCC -D_NOTHERMORAW -D_NOSQLITE 

INCPATHS = -I../include -I../mstoolkit/include

OBJS = dslim_query.o dslim.o lbe_internal.o utils.o mods.o msquery.o dslim_fileout.o
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
