#.PHONY: Make-depend

WITH_OPENMP = 0     # warning does not work
FASTANDFURIOUS = 0 #1

CXX = g++
CC = gcc

LD            := $(CXX)

SRC_DIR = $(PWD)/src
PROG_DIR = $(PWD)/prog
UTL_DIR = $(PWD)/ utl
OBJ_DIR = $(PWD)/obj
LIB_DIR = $(PWD)/lib
BIN_DIR = $(PWD)/bin

CINT = rootcint
CXXFLAGS += -std=c++11 -O3 -Wall -Wextra -fPIC
CXXFLAGS += $(shell root-config --cflags)
CXXFLAGS += -I$(SRC_DIR) -I.
LDFLAGS  += -fPIC -ggdb3 -Wall
LDFLAGS  += $(shell root-config --libs)  -Wl,--no-as-needed -fPIC
#LDFLAGS  += $(shell root-config --ldflags) -lMinuit -lGeom
LDFLAGS  += $(shell root-config --ldflags) -lTreePlayer -lMinuit 
LDFLAGS  += -lgsl -lgslcblas
LDFLAGS  += -Wl,-rpath,$(LIB_DIR)
SOFLAGS   = -Wl,--no-as-needed -fPIC -ggdb3 -Wall -shared

ifeq ($(WITH_OPENMP), 1)
   CXXFLAGS += -D_WITH_OPENMP_
   LDFLAGS  += -fopenmp
   CXXFLAGS  += -fopenmp
endif

ifeq ($(FASTANDFURIOUS), 1)
   CXXFLAGS += -D_FASTANDFURIOUS_
endif


OBJS  = $(patsubst $(SRC_DIR)/%LinkDef.h, $(OBJ_DIR)/%Dict.o, $(wildcard $(SRC_DIR)/*LinkDef.h))
OBJS += $(patsubst $(SRC_DIR)/%.cc, $(OBJ_DIR)/%.o, $(wildcard $(SRC_DIR)/*.cc))
SOURCES := $(shell find $(SRC_DIR) -name '*.cc' ! -path '.svn/*')
SOURCES += $(shell find $(PROG_DIR) -name '*.cc' ! -path '.svn/*')

EXE = $(patsubst $(PROG_DIR)/%.cc, $(BIN_DIR)/%, $(wildcard $(PROG_DIR)/*.cc))
CLIBS = $(LIB_DIR)/libProp.so


all: $(CLIBS) $(EXE)

$(BIN_DIR)/runFit: $(PROG_DIR)/runFit.cc macros/fit.C
	$(CXX) -o $@ $^ $(CXXFLAGS)  $(LDFLAGS) -Llib -lProp

$(BIN_DIR)/runFit_massgroups: $(PROG_DIR)/runFit_massgroups.cc macros/fit.C
	$(CXX) -o $@ $^ $(CXXFLAGS)  $(LDFLAGS) -Llib -lProp

$(BIN_DIR)/%: $(PROG_DIR)/%.cc
	$(CXX) -o $@ $^ $(CXXFLAGS)  $(LDFLAGS) -Llib -lProp

$(CLIBS): $(OBJS)
	@mkdir -p $(LIB_DIR)
	$(LD) $(SOFLAGS) $(LDFLAGS) $^ -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc
	@mkdir -p $(OBJ_DIR)
	$(CXX) -c $(CXXFLAGS) $^ -o $@

$(SRC_DIR)/%Dict.cc: $(SRC_DIR)/%.h $(SRC_DIR)/%LinkDef.h
	@(echo generating $@ dictionary)
	($(CINT) -f $@ -c $(CINTCXXFLAGS)   $^)
	@mkdir -p lib
	@mv src/*Dict_rdict.pcm lib || true

clean:
	@rm -f *~ $(OBJS) core \
	       $(SRC_DIR)/*Dict.* $(SRC_DIR)/*~
	@rm -rf $(OBJ_DIR)
	@rm -rf $(LIB_DIR)
	@rm -f $(EXE)
	@rm -f Make-depend

Make-depend:
	$(CPP) $(CXXFLAGS) $(ROOTCXXFLAGS) $(CPPFLAGS) -MM $(SOURCES) > $@

#ifneq ($(MAKECMDGOALS),clean)
-include Make-depend
#endif
