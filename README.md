## ![image-20211220153709911](./RabbitQCPlus.png)

# RabbitQCPlus

A modern quality control tool for sequencing data.

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
./RabbitQCPlus -w 4 -i in.fastq --TGS --usePugz --pugzThread 2
```

# Options

For more help information, please refer to `./RabbitQCPlus -h`.

# Build

**For Linux and OSX:**

```
cd RabbitQCPlus && make clean && make
```

**For Windows:**

//TODO

