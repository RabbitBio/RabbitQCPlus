INCDIRS := include ../include

CC := g++
NVCC := nvcc
STD := c++14
OPT := 3
CCFLAGS := -O$(OPT) -std=$(STD) -Wall -Wextra -fopenmp
XCCFLAGS := $(addprefix -Xcompiler ,$(CCFLAGS))
NVCCGENCODE = -arch=sm_35
NVCCFLAGS := -O$(OPT) -std=$(STD) $(NVCCGENCODE) -ccbin $(CC) $(addprefix -Xcompiler ,$(CCFLAGS)) --expt-extended-lambda

INCPARAMS := $(addprefix -I, $(INCDIRS))

all: test

test: build/packed_types_test.o build/catch_main.o | bin
	$(NVCC) $(NVCCFLAGS) build/packed_types_test.o build/catch_main.o -o bin/$@

build/packed_types_test.o: src/packed_types_test.cu | build
	$(NVCC) $(NVCCFLAGS) $(INCPARAMS) -c $< -o $@

build/catch_main.o: src/catch_main.cu | build
	$(NVCC) $(NVCCFLAGS) $(INCPARAMS) -c $< -o $@

debug: OPT := 0
debug: CCFLAGS := -O$(OPT) -std=$(STD) -Wall -Wextra -fopenmp
debug: XCCFLAGS := $(addprefix -Xcompiler ,$(CCFLAGS))
debug: NVCCFLAGS := -O$(OPT) -std=$(STD) -ccbin $(CC) $(XCCFLAGS) $(NVCCGENCODE) --expt-extended-lambda -g -Xptxas -v -UNDEBUG -DDEBUG
debug: all

profile: NVCCFLAGS += -lineinfo -g -Xptxas -v -UNDEBUG
profile: all

clean:
	$(RM) -r build bin

build bin:
	mkdir -p $@

.PHONY: clean all build bin

