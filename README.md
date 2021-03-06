## ![image-20211220153709911](./RabbitQCPlus.png)

# RabbitQCPlus Features

A modern quality control tool for sequencing data.
- Single-threaded performance is improved by a factor of 2 at least.

- Solving the performance issues when processing gz files (more than 5x speedups compared with SOTA).

- Improving the efficiency of the time-consuming over-representation module by a factor of 5.

# Build

RabbitQCPlus can only support 64-bit Linux Systems.

### Dependancy

- gcc7.3.0 or newer
- [zlib](https://zlib.net/)

### Compilation

```bash
git clone https://github.com/RabbitBio/RabbitQCPlus.git
cd RabbitQCPlus
make
```

You can modify the compilation parameters in the Makefile to select the instruction set used for vectorization. 

``-DVec512`` means using the avx512 instruction set, and ``-DVec256`` means using the avx2 instruction set; otherwise, let the compiler choose.

You can also specify to use igzip in isal by default when processing compressed files (make sure it has been installed on the machine), which is much faster than the default zlib, but not as fast as pugz (>= 4 threads).

Specifically, add ``-DUSE_IGZIP`` to ``CXXFLAGS`` and ``-lisal`` to ``LIBS`` in the Makefile.

# Simple usage

## For next generation sequencing data

- For SE (not compressed)

```bash
./RabbitQCPlus -w 8 -i in1.fastq -o p1.fastq
```

- For SE (gzip compressed)

```bash
./RabbitQCPlus -w 4 -i in1.fastq.gz -o p1.fastq.gz --usePugz --pugzThread 2 --usePigz --pigzThread 4
```

- For PE (not compressed)

```bash
./RabbitQCPlus -w 8 -i in1.fastq -I in2.fastq -o p1.fastq -O p2.fastq
```

- For PE (gzip compressed)

```bash
./RabbitQCPlus -w 4 -i in1.fastq.gz  -I in2.fastq.gz -o p1.fastq.gz -O p2.fastq.gz --usePugz --pugzThread 2 --usePigz --pigzThread 2
```

## For third generation sequencing data

- not compressed

```bash
./RabbitQCPlus -w 4 -i in.fastq --TGS
```

- gzip compressed

```bash
./RabbitQCPlus -w 4 -i in.fastq.gz --TGS --usePugz --pugzThread 2
```

# Options

For more help information, please refer to `./RabbitQCPlus -h`.



