CXX = g++
ROOTFLAGS := $(shell root-config --cflags --glibs)

INCDIR = include
SRCDIR = src
BINDIR = bin

SRCFILES := $(wildcard $(SRCDIR)/*.cxx)
INCFILES := $(wildcard $(INCDIR)/*.h)

MACROS := $(notdir $(basename $(SRCFILES)))

all: $(MACROS)

$(MACROS): %: $(SRCDIR)/%.cxx $(INCFILES)
	@echo "building" $(BINDIR)/$@
	@mkdir -p $(BINDIR)
	@$(CXX) -o $(BINDIR)/$@ $(SRCDIR)/$@.cxx $(ROOTFLAGS)

clean:
	@rm -rf $(BINDIR)
