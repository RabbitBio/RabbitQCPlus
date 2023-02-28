NVCC := nvcc
STD := c++14

OPT := -O3
DDEBUG := -O0 -g -Xptxas -v -UNDEBUG -DDEBUG
PROFILE := -O3 -lineinfo -g -Xptxas -v -UNDEBUG

CXXFLAGS = -std=$(STD) -Wall -Wextra -fopenmp
NVCCGENCODE = -arch=sm_35
NVCCFLAGS = $(OPT) -std=$(STD) $(NVCCGENCODE) -ccbin $(CXX) $(addprefix -Xcompiler ,$(CXXFLAGS)) --expt-extended-lambda

INCDIRS := include
INCPARAMS = $(addprefix -I, $(INCDIRS))

HELPERS := \
    ../include/cuda_helpers.cuh \
    ../include/hpc_helpers.h \
    ../include/io_helpers.h \
    ../include/timers.cuh \
    ../include/packed_types.cuh

.PHONY: all debug profile clean

all: helpers_example packed_types_example

debug: OPT := $(DEBUG)
debug: all

profile: OPT := $(PROFILE)
profile: all

helpers_example: helpers_example.cu $(HELPERS)
	$(NVCC) $(NVCCFLAGS) helpers_example.cu -o helpers_example

packed_types_example: packed_types_example.cu $(HELPERS)
	$(NVCC) $(NVCCFLAGS) packed_types_example.cu -o packed_types_example

clean:
	$(RM) -r helpers_example packed_types_example
