PRINTDEBUG := 0

GCCVERSION = $(shell gcc -dumpfullversion -dumpversion)

GCC_GTEQ_485 := $(shell expr `gcc -dumpfullversion -dumpversion | sed -e 's/\.\([0-9][0-9]\)/\1/g' -e 's/\.\([0-9]\)/0\1/g' -e 's/^[0-9]\{3,4\}$$/&00/'` \>= 40805)
GCC_GTEQ_700 := $(shell expr `gcc -dumpfullversion -dumpversion | sed -e 's/\.\([0-9][0-9]\)/\1/g' -e 's/\.\([0-9]\)/0\1/g' -e 's/^[0-9]\{3,4\}$$/&00/'` \>= 70000)
GCC_GTEQ_1000 := $(shell expr `gcc -dumpfullversion -dumpversion | sed -e 's/\.\([0-9][0-9]\)/\1/g' -e 's/\.\([0-9]\)/0\1/g' -e 's/^[0-9]\{3,4\}$$/&00/'` \>= 100000)
ifeq ($(PRINTDEBUG),1)
$(info The detected gcc version is $(GCCVERSION))
$(info The detected gcc version is >= 4.8.5 : $(GCC_GTEQ_485))
$(info The detected gcc version is >= 7.0.0 : $(GCC_GTEQ_700))
$(info The detected gcc version is >= 10.0.0 : $(GCC_GTEQ_1000))
endif

InstructSet :=

ifneq ($(MAKECMDGOALS),clean)
$(info The detected gcc version is $(GCCVERSION) )
endif


ifneq ($(MAKECMDGOALS),clean)
ifeq ($(GCC_GTEQ_485),1)
ifeq ($(GCC_GTEQ_700),1)
$(info gcc version >= 7.0.0, now use avx512 instruction set to accelerate code)
InstructSet := -DVec512
else
$(info gcc version >= 4.8.5 but < 7.0.0, now use avx2 instruction set to accelerate code)
InstructSet := -DVec256
endif
else
$(info gcc version < 4.8.5, now let the compiler automatically select the instruction set. However, this may affect performance and we recommend that you install gcc with >= 4.8.5)
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

