PRINTDEBUG := 0

$(info            )

DIR_INC := ./inc
DIR_SRC := ./src
SLAVE_DIR_SRC := ./slave
DIR_OBJ := ./obj


PREFIX ?= /usr/local
BINDIR ?= $(PREFIX)/bin
INCLUDE_DIRS ?= 
LIBRARY_DIRS ?=

SRC := $(wildcard ${DIR_SRC}/*.cpp)
OBJ := $(patsubst %.cpp,${DIR_OBJ}/%.o,$(notdir ${SRC}))


SRC2 := $(wildcard ${DIR_SRC}/*.c)
OBJ += $(patsubst %.c,${DIR_OBJ}/%.o,$(notdir ${SRC2}))


SRC3 := $(wildcard ${SLAVE_DIR_SRC}/*.cpp)
OBJ += $(patsubst %.cpp,${DIR_OBJ}/%.o,$(notdir ${SRC3}))


TARGET := RabbitQCPlus

BIN_TARGET := ${TARGET}


CXX = mpicxx
#CXX = swg++


CXXFLAGS := $(InstructSet)
CXXFLAGS +=  -DVerbose -std=c++11 -I./ -I./common -g -O3 -w


CXX2 = mpicc
#CXX2 = swgcc

CXXFLAGS2 :=  -g -O3 -w -Wextra -Wno-unknown-pragmas -Wcast-qual

LIBS := -lz -lpthread -lrt


LD_FLAGS := $(foreach librarydir,$(LIBRARY_DIRS),-L$(librarydir)) $(LIBS)

${BIN_TARGET}:${OBJ}
		$(CXX) -mdynamic $^ -o $@ $(LD_FLAGS)

${DIR_OBJ}/%.o:${DIR_SRC}/%.cpp
		$(CXX) -mhost -c $< -o $@ $(CXXFLAGS)
${DIR_OBJ}/%.o:${SLAVE_DIR_SRC}/%.cpp
		$(CXX) -mslave -msimd -c $< -o $@ $(CXXFLAGS)
${DIR_OBJ}/%.o:${DIR_SRC}/%.c
		$(CXX2) -mhost -c $< -o $@ $(CXXFLAGS2) 



.PHONY:clean
clean:
	rm $(DIR_OBJ)/*.o
	rm $(TARGET)

install:
	install $(TARGET) $(BINDIR)/$(TARGET)
	@echo "Installed."

