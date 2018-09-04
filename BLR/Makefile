CXXFLAGS += -std=c++11 -O3 -Wall -Wextra -fPIC
CXXFLAGS += $(shell root-config --cflags)
CXXFLAGS += -I.. -I.
LDFLAGS  += -fPIC -ggdb3 -Wall
LDFLAGS  += $(shell root-config --libs)  -Wl,--no-as-needed -fPIC
LDFLAGS  += $(shell root-config --ldflags) -lMinuit -lGeom

EXE = $(patsubst %.cxx, %, $(wildcard *.cxx))
EXE += $(patsubst %.c, %, $(wildcard *.c))
OBJ = $(patsubst %.cc, %.o, $(wildcard *.cc))

all: $(EXE) 

%.o: %.cc
	$(CXX) -c -o $@ $< $(CXXFLAGS) $(CPPFLAGS)

%: %.cxx $(OBJ)
	$(CXX) -o $@ $< $(OBJ) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS)

%: %.c
	$(CC) -o $@ $< -lm


clean:
	@rm -f $(LIB) $(OBJS) $(EXE) $(XML) $(OBJ) *Dict.cc Make-depend

