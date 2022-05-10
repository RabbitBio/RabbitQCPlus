DIR_INC := ./inc
DIR_SRC := ./src
DIR_OBJ := ./obj

PREFIX ?= /usr/local
BINDIR ?= $(PREFIX)/bin
INCLUDE_DIRS ?= 
LIBRARY_DIRS ?=

SRC := $(wildcard ${DIR_SRC}/*.cpp)
OBJ := $(patsubst %.cpp,${DIR_OBJ}/%.o,$(notdir ${SRC}))


SRC2 := $(wildcard ${DIR_SRC}/*.c)
OBJ += $(patsubst %.c,${DIR_OBJ}/%.o,$(notdir ${SRC2}))


TARGET := RabbitQCPlus

BIN_TARGET := ${TARGET}


CXX = g++

# you can add -DUSE_IGZIP to CXXFLAGS and -lisal to LIBS, to use igzip by default

# -DVec512 means using avx512 instruction set
# -DVec256 means using avx2 instruction set
# otherwise, let the compiler choose

# you can add -DVerbose to print more log information

CXXFLAGS := -DVec256 -std=c++11 -I./ -I./common -march=native  -g -O3  -w -fopenmp


CXX2 = gcc
CXXFLAGS2 := -g -O3 -w -Wextra -Wno-unknown-pragmas -Wcast-qual

LIBS := -lz -lpthread  -fopenmp  -lrt 


LD_FLAGS := $(foreach librarydir,$(LIBRARY_DIRS),-L$(librarydir)) $(LIBS)




${BIN_TARGET}:${OBJ}
	$(CXX) $(OBJ) -o $@ $(LD_FLAGS)

${DIR_OBJ}/%.o:${DIR_SRC}/%.cpp
	$(CXX) -c $< -o $@ $(CXXFLAGS)


${DIR_OBJ}/%.o:${DIR_SRC}/%.c
	$(CXX2) $(CXXFLAGS2) -c $< -o $@

.PHONY:clean
clean:
	rm $(DIR_OBJ)/*.o
	rm $(TARGET)

install:
	install $(TARGET) $(BINDIR)/$(TARGET)
	@echo "Installed."

