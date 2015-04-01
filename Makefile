.PHONY: Make-depend

LD            := $(CXX)

SRC_DIR = src
PROG_DIR = prog
UTL_DIR = utl
OBJ_DIR = obj
LIB_DIR = lib
BIN_DIR = bin

CINT = rootcint
CXXFLAGS += -std=c++0x -O2 -Wall -Wextra -fPIC
CXXFLAGS += $(shell root-config --cflags)
CXXFLAGS += -I$(SRC_DIR) -I.
LDFLAGS  += -fPIC -ggdb3 -Wall
LDFLAGS  += $(shell root-config --libs)  -Wl,--no-as-needed -fPIC
LDFLAGS  += $(shell root-config --ldflags) -lMinuit -lGeom
LDFLAGS  += -lgsl -lgslcblas
SOFLAGS   = -Wl,--no-as-needed -fPIC -ggdb3 -Wall -shared


OBJS  = $(patsubst $(SRC_DIR)/%LinkDef.h, $(OBJ_DIR)/%Dict.o, $(wildcard $(SRC_DIR)/*LinkDef.h))
OBJS += $(patsubst $(SRC_DIR)/%.cc, $(OBJ_DIR)/%.o, $(wildcard $(SRC_DIR)/*.cc))
SOURCES := $(shell find $(SRC_DIR) -name '*.cc' ! -path '.svn/*')
SOURCES += $(shell find $(PROG_DIR) -name '*.cc' ! -path '.svn/*')

EXE = $(patsubst $(PROG_DIR)/%.cc, $(BIN_DIR)/%, $(wildcard $(PROG_DIR)/*.cc))
CLIBS = $(LIB_DIR)/libProp.so


all: Make-depend $(CLIBS) $(EXE)

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
