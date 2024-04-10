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


SRC4 := $(wildcard ${SLAVE_DIR_SRC}/lib/*.c)
OBJ += $(patsubst %.c,${DIR_OBJ}/%.o,$(notdir ${SRC4}))



TARGET := RabbitQCPlus
TARGET_MPI := RabbitQCPlus_mpi

BIN_TARGET := ${TARGET}
BIN_TARGET_MPI := ${TARGET_MPI}


CXX = mpicxx
#CXX = swg++


CXXFLAGS := $(InstructSet)
CXXFLAGS += -DVerbose -std=c++11 -I./ -I./common -g -O3 -w


CXX2 = mpicc
#CXX2 = swgcc

CXXFLAGS2 := -g -O3 -w

LIBS := -lz -lpthread -lrt


LD_FLAGS := $(foreach librarydir,$(LIBRARY_DIRS),-L$(librarydir)) $(LIBS)

#all: ${BIN_TARGET} ${BIN_TARGET_MPI}
all: ${BIN_TARGET_MPI}

${BIN_TARGET}:${OBJ}
		$(CXX) -mhybrid $^ -o $@ $(LD_FLAGS)

${BIN_TARGET_MPI}: ${OBJ}
		$(CXX) -mdynamic $^ -o $@ $(LD_FLAGS)

${DIR_OBJ}/%.o:${DIR_SRC}/%.cpp
		$(CXX) -mhost -c $< -o $@ $(CXXFLAGS)
${DIR_OBJ}/%.o:${SLAVE_DIR_SRC}/%.cpp
		$(CXX) -mslave -msimd -c $< -o $@ $(CXXFLAGS)
${DIR_OBJ}/%.o:${DIR_SRC}/%.c
		$(CXX2) -mhost -c $< -o $@ $(CXXFLAGS2) 
${DIR_OBJ}/%.o:${SLAVE_DIR_SRC}/lib/%.c
		$(CXX2) -mslave -msimd -c $< -o $@ $(CXXFLAGS2) 





.PHONY:clean
clean:
	rm -rf $(DIR_OBJ)/*.o
	#rm -rf $(TARGET)
	rm -rf $(TARGET_MPI)

install:
	install $(TARGET) $(BINDIR)/$(TARGET)
	@echo "Installed."

