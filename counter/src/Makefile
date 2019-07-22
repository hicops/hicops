SRC_ROOT = $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

CXXFLAGS +=

INCPATHS = -I../include -I../mstoolkit/include

OBJS = lbe_internal.o utils.o mods.o
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
