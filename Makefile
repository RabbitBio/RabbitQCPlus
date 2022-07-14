GCCVERSION = $(shell gcc --version | grep ^gcc | sed 's/^.* //g')

GCCVERSIONGTEQ4 := $(shell expr `gcc -dumpversion | cut -f1 -d.` \>= 4)
GCCVERSIONGTEQ5 := $(shell expr `gcc -dumpversion | cut -f1 -d.` \>= 5)
GCCVERSIONGTEQ6 := $(shell expr `gcc -dumpversion | cut -f1 -d.` \>= 6)
GCCVERSIONGTEQ7 := $(shell expr `gcc -dumpversion | cut -f1 -d.` \>= 7)
GCCVERSIONGTEQ8 := $(shell expr `gcc -dumpversion | cut -f1 -d.` \>= 8)
GCCVERSIONGTEQ9 := $(shell expr `gcc -dumpversion | cut -f1 -d.` \>= 9)
GCCVERSIONGTEQ10 := $(shell expr `gcc -dumpversion | cut -f1 -d.` \>= 10)

InstructSet :=

ifneq ($(MAKECMDGOALS),clean)
$(info The detected gcc version is $(GCCVERSION) )
endif


ifneq ($(MAKECMDGOALS),clean)
ifeq ($(GCCVERSIONGTEQ5),1)
ifeq ($(GCCVERSIONGTEQ7),1)
$(info gcc version >= 7, now use avx512 instruction set to accelerate code)
InstructSet := -DVec512
else
$(info gcc version >= 5 but < 7, now use avx2 instruction set to accelerate code)
InstructSet := -DVec256
endif
else
$(info gcc version < 5, now let the compiler automatically select the instruction set)
endif
endif


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

CXXFLAGS := $(InstructSet)
CXXFLAGS += -DVerbose -std=c++11 -I./ -I./common -march=native  -g -O3  -w -fopenmp


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

