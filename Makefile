PRINTDEBUG := 0
InstructSet :=

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


CXX = sw5g++


CXXFLAGS := $(InstructSet)
CXXFLAGS += -DVerbose -std=c++11 -I./ -I./common -g -O3 -w -fopenmp


CXX2 = sw5gcc
CXXFLAGS2 := -g -O3 -w -Wextra -Wno-unknown-pragmas -Wcast-qual

LIBS := -static -lz -Wl,-z,muldefs -lpthread -fopenmp -lrt


LD_FLAGS := $(foreach librarydir,$(LIBRARY_DIRS),-L$(librarydir)) $(LIBS)

${BIN_TARGET}:${OBJ}
		python3 xlink.py $(CXX) -mhybrid $^ -o $@ $(LD_FLAGS)

${DIR_OBJ}/%.o:${DIR_SRC}/%.cpp
		$(CXX) -mhost -msimd -c $< -o $@ $(CXXFLAGS)
${DIR_OBJ}/%.o:${SLAVE_DIR_SRC}/%.cpp
		$(CXX) -mslave -msimd -c $< -o $@ $(CXXFLAGS)
${DIR_OBJ}/%.o:${DIR_SRC}/%.c
		$(CXX2) -mhost -msimd -c $< -o $@ $(CXXFLAGS2) 



#${BIN_TARGET}:${OBJ}
#		python3 xlink.py $(CXX) -mhybrid $^ -o $@ $(LD_FLAGS)
#
#${DIR_OBJ}/%.o:${DIR_SRC}/%.cpp
#		$(CXX) -mhost -msimd -c $< -o $@ $(CXXFLAGS)
#
#${DIR_OBJ}/%.o:${SLAVE_DIR_SRC}/%.cpp
#		$(CXX) -mslave -msimd -c $< -o $@ $(CXXFLAGS)
#
#${DIR_OBJ}/%.o:${DIR_SRC}/%.c
#		$(CXX2) -mhost -msimd -c $< -o $@ $(CXXFLAGS2)
#


.PHONY:clean
clean:
	rm $(DIR_OBJ)/*.o
	rm $(TARGET)

install:
	install $(TARGET) $(BINDIR)/$(TARGET)
	@echo "Installed."

