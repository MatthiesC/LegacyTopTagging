CC = g++
ROOTFLAGS = `root-config --cflags --glibs`


INCDIR   = include
SRCDIR   = src

SOURCES  := $(wildcard $(SRCDIR)/*.cxx)
INCLUDES := $(wildcard $(INCDIR)/*.h)

MACROS_temp := $(basename $(SOURCES))
MACROS := $(notdir $(MACROS_temp) )

all: $(MACROS)


$(MACROS): % : $(SRCDIR)/%.cxx $(INCLUDES)
	@echo 'building' $@
	@$(CC) -o $@ $(SRCDIR)/$@.cxx $(ROOTFLAGS)


clean:
	@rm -f $(MACROS)
