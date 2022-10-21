## ![logo](./pics/RabbitQCPlus.png)

# RabbitQCPlus Features

A modern quality control tool for sequencing data.
- Single-threaded performance is improved by a factor of 2 at least.

- Solving the performance issues when processing gz files (more than 4x speedups compared with SOAPnuke).

- Improving the efficiency of the time-consuming over-representation module by a factor of 5.

# Build

RabbitQCPlus can only support 64-bit Linux Systems.

### Dependancy

- gcc 4.8.5 or newer
- [zlib](https://zlib.net/)

### Compilation

```bash
git clone https://github.com/RabbitBio/RabbitQCPlus.git
cd RabbitQCPlus
make -j4
```
To improve the robustness of the software, we have implemented different software versions for different vectorized instruction sets. RabbitQCPlus can automatically detect the system CPU instruction set and compiler version at compile time to select the appropriate software version.

You can also specify the instruction set you want to use by manually modifying the ``InstructSet`` in the ``Makefile``. ``-DVec512`` means using the avx512 instruction set, and ``-DVec256`` means using the avx2 instruction set; otherwise, let the compiler choose.

# Simple usage

## For next generation sequencing data

- For SE (not compressed)

```bash
./RabbitQCPlus -w 8 -i in1.fastq -o p1.fastq
```

- For SE (gzip compressed)

```bash
./RabbitQCPlus -w 8 -i in1.fastq.gz -o p1.fastq.gz
```

- For PE (not compressed)

```bash
./RabbitQCPlus -w 8 -i in1.fastq -I in2.fastq -o p1.fastq -O p2.fastq
```

- For PE (gzip compressed)

```bash
./RabbitQCPlus -w 16 -i in1.fastq.gz  -I in2.fastq.gz -o p1.fastq.gz -O p2.fastq.gz
```

## For third generation sequencing data

- not compressed

```bash
./RabbitQCPlus -w 4 -i in.fastq --TGS
```

- gzip compressed

```bash
./RabbitQCPlus -w 6 -i in.fastq.gz --TGS
```

# Options

For more help information, please refer to `./RabbitQCPlus -h`.



# Performance results

Experiments have been conducted on a Linux server with 2 Intel Xeon Platinum 8260 CPUs, 1.5T RAM, 2T SSD, running Ubuntu 20.04 and GCC9. We have compared the performance of RabbitQCPlus with RabbitQC (v0.0.1), fastp (v0.23.2), SOAPnuke (v2.1.7), Trimmomatic (v0.40), and FASTQC (v0.11.9) using 370 million Illumina sequencing reads ([SRR7963242](https://www.ncbi.nlm.nih.gov/sra/?term=SRR7963242)).

We use default parameters of each application except for the number of threads on both plain and gzip-compressed FASTQ files.
When using default parameters, RabbitQCPlus performs the same or more comprehensive quality control operations compared to other applications.

## Plain FASTQ files

<img src="pics/plain.png" alt="plain" style="zoom:50%;" />

## Gzip-compressed FASTQ files

<img src="pics/gzip.png" alt="gzip" style="zoom:50%;" />

# Visual output

We visualized the information before and after data filtering, and [here](https://yanlifeng.github.io/someTest/example.html) is an example.
