IDIR=include
CXX=g++
CFLAGS=-I$(IDIR) -O0 -fopenmp --std=c++17 -pg

OBJDIR=obj
SRCDIR=src
LDIR=lib
LDFLAGS=-lm -fopenmp -pg

DEPS = $(wildcard $(IDIR)/*.h)

SRCS=$(wildcard $(SRCDIR)/*.cc)
OBJ = $(patsubst $(SRCDIR)/%.cc,$(OBJDIR)/%.o,$(SRCS))


$(OBJDIR)/%.o: $(SRCDIR)/%.cc $(DEPS)
	@mkdir -p $(@D)
	$(CXX) -c -o $@ $< $(CFLAGS)

gsim: $(OBJ)
	$(CXX) -o $@ $^ $(LDFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(OBJDIR)/*.o *~ core $(INCDIR)/*~
