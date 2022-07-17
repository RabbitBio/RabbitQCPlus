PRINTDEBUG := 0

#query gcc version
GCCVERSION := $(shell gcc -dumpfullversion -dumpversion)
GCC_GTEQ_485 := $(shell expr `gcc -dumpfullversion -dumpversion | sed -e 's/\.\([0-9][0-9]\)/\1/g' -e 's/\.\([0-9]\)/0\1/g' -e 's/^[0-9]\{3,4\}$$/&00/'` \>= 40805)
GCC_GTEQ_700 := $(shell expr `gcc -dumpfullversion -dumpversion | sed -e 's/\.\([0-9][0-9]\)/\1/g' -e 's/\.\([0-9]\)/0\1/g' -e 's/^[0-9]\{3,4\}$$/&00/'` \>= 70000)
GCC_GTEQ_1000 := $(shell expr `gcc -dumpfullversion -dumpversion | sed -e 's/\.\([0-9][0-9]\)/\1/g' -e 's/\.\([0-9]\)/0\1/g' -e 's/^[0-9]\{3,4\}$$/&00/'` \>= 100000)

#select which instruction set to use
sse2F := 0
sse4_1F := 0
avxF := 0
avx2F := 0
popcntF := 0
avx512fF := 0
avx512vlF := 0
avx512bwF := 0
ifneq ("$(shell lscpu | grep ^Flags | grep -w sse2)","")
sse2F := 1
endif
ifneq ("$(shell lscpu | grep ^Flags | grep -w sse4_1)","")
sse4_1F := 1
endif
ifneq ("$(shell lscpu | grep ^Flags | grep -w avx)","")
avxF := 1
endif
ifneq ("$(shell lscpu | grep ^Flags | grep -w avx2)","")
avx2F := 1
endif
ifneq ("$(shell lscpu | grep ^Flags | grep -w popcnt)","")
popcntF := 1
endif
ifneq ("$(shell lscpu | grep ^Flags | grep -w avx512f)","")
avx512fF := 1
endif
ifneq ("$(shell lscpu | grep ^Flags | grep -w avx512bw)","")
avx512bwF := 1
endif
ifneq ("$(shell lscpu | grep ^Flags | grep -w avx512vl)","")
avx512vlF := 1
endif
avx2Set := 0
avx512Set := 0
ifeq ($(sse4_1F)_$(sse2F)_$(avxF)_$(avx2F)_$(popcntF),1_1_1_1_1)
avx2Set := 1
endif
ifeq ($(avxF)_$(avx2F)_$(popcntF)_$(avx512fF)_$(avx512bwF)_$(avx512vlF),1_1_1_1_1_1)
avx512Set := 1
endif

ifeq ($(PRINTDEBUG),1)
$(info The detected gcc version is $(GCCVERSION))
$(info The detected gcc version is >= 4.8.5 : $(GCC_GTEQ_485))
$(info The detected gcc version is >= 7.0.0 : $(GCC_GTEQ_700))
$(info The detected gcc version is >= 10.0.0 : $(GCC_GTEQ_1000))
$(info This cpu has sse2 flag or not : $(sse2F))
$(info This cpu has sse4_1 flag or not : $(sse4_1F))
$(info This cpu has avx flag or not : $(avxF))
$(info This cpu has avx2 flag or not : $(avx2F))
$(info This cpu has avx512f flag or not : $(avx512fF))
$(info This cpu has avx512bw flag or not : $(avx512bwF))
$(info This cpu has avx512vl flag or not : $(avx512vlF))
$(info This cpu has popcnt flag or not : $(popcntF))
endif

InstructSet :=

ifneq ($(MAKECMDGOALS),clean)
$(info The detected gcc version is $(GCCVERSION) )
$(info         )
endif



ifneq ($(MAKECMDGOALS),clean)
ifeq ($(GCC_GTEQ_700)_$(avx512Set),1_1)
$(info Based on the detected gcc version and cpuflags, it was decided to use the avx512 instruction set to speed up the program)
InstructSet := -DVec512
else ifeq ($(GCC_GTEQ_485)_$(avx2Set),1_1)
$(info Based on the detected gcc version and cpuflags, it was decided to use the avx2 instruction set to speed up the program)
InstructSet := -DVec256
else
$(info Based on the detected gcc version and cpuflags, it was decided to let the compiler do automatic vectorization)
endif
endif

#select instruction set according to gcc version and cpuflags
#gcc >=7 && cpuflags include avx avx2 avx512f avx512vl avx512bw popcnt ----> use avx512
#gcc >=4.8.5 && cpuflags include sse4_1 sse2 avx avx2 popcnt ----> use avx2
#avx512 is preferred
#users can also manually modify 'InstructSet' to use either AVx2 or AVx512 or neither

$(info            )

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
CXXFLAGS += -std=c++11 -I./ -I./common -march=native  -g -O3  -w -fopenmp


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

