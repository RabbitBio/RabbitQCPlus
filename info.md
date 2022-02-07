- [x] Filter
- [x] Stats
- [x] Trimer->front tail trim
- [x] Trimer->auto adapter trim
- [x] Trimer->input adapter trim
- [x] Umi
- [x] PolyX
- [x] Overrepresented
- [x] Duplicate
- [ ] new Duplicate
- [ ] Deduplicate
- [x] Draw
- [x] 内存泄漏
- [x] check 正确性
- [ ] option里面的error detect
- [x] support long reads
- [x] support reading from STDIN and writing to STDOUT
- [x] support interleaved input
- [ ] ~~split the output to multiple files~~
- [x] Insert size estimation
- [x] support zip input output
- [x] optimize gzip
- [ ] optimize memory use when out is zip
- [x] Phred64
- [x] add pac gz in QC
- [ ] add 3ed function
- [ ] parallel in producer ??
- [ ] optimize write part --- reduce approx
- [ ] add comparison between muti datasets
- [ ] add N50 ?

### QC的几个问题

👇 

## AfterQC

- Different from most tools, AfterQC analyses the overlapping of paired sequences for pair-end sequencing data. Based on overlapping analysis, AfterQC can detect and cut adapters, and furthermore it gives a novel function to correct wrong bases in the overlapping regions. 

- filtersoutbadreads,detectsandeliminatessequencer’sbubble effects, trims reads at front and tail, detects the sequencing errors and corrects part of them, and finally outputs clean data and generates HTML reports with interactive figures

  

- Bubble detection and visualisation：使用聚类的方式检测出气泡在哪并且进行可视化展示

- Automatic trimming：即修建开头和结尾的低质量碱基，为了避免局部策略的弊端，提出了一种全局策略：先做一遍整体的统计，画出质量曲线；然后从中间开始往两边扫直到遇到不正常的部分(  

  

  1), too high or too low of mean base content percentages (i.e higher than 40%, or lower than 15%); 

  2), too significant change of mean base content percentages (i.e, ±10% change com- paring to neighbour cycle); 

  3), too high or too low of mean GC percentages (i.e higher than 70%, or lower than 30%); 

  4), too low of mean quality (i.e. less then Q20). Figure 4 gives an example how automatic trimming works.

  

     )->Before trimming happens, AfterQC will do pre-filtering quality control to calculate the base content and quality curves. Our algorithm initialises the central cycle as a good cycle, and then expands the good region by scanning the base content and quality curves cycle by cycle, until it meets the front or end, or meet a cycle con- sidered as abnormal. Then the cycles in the good region will be kept, and the rest cycles in the front and tail will be trimmed. 

- Filtering：quality filters->count the number of low quality bases or N, calculate the mean quality of each read, and then determine whether to keep or discard this read.   polyX filters://TODO

- Overlapping analysis and error correction：枚举offest找到双端数据的重叠部分（计算最小的编辑距离//TODO），然后进行修正->如果海明距离==编辑距离，就直接修正没对不一样的碱基（质量分低的改成高的）//TODO

- Sequencing error profiling：即根据👆提到的低质量改为高质量做下统计，然后分析一下，发现这个依赖于测序仪。

- Automatic adapter cutting：这里只讲到了双端数据去除adapter，方式就是基于overloaping analyze如果找到的best offest < 0，就断定出现了下图的状况：![image-20210701102351548](/Users/ylf9811/Library/Application Support/typora-user-images/image-20210701102351548.png)

  （不过前提必须是adapter添加在3‘端）

  这里提一下RabbitQC的adapter修剪：对于单端数据，开始先拿前100w条read检测adapter，然后每次直接修剪；对于双端的数据，过程和afterQC类似，先进行OverlapAnalysis（如果打开了correct参数会接着进行修正，包含了Sequencing error profiling部分，默认是不开的），如果offest<0，按照👆的方法去除adapter，否则，和单端一样（但是双端数据默认一开始不做adapter的检测）

- Quality profiling：strand bias profiling to reflect amplification bias, and per-cycle dis- continuity profiling to reflect sequencing quality instability. 前者是统计短kmer正反向计数是不是基本相同，后者统计不连续性是不是相对稳定。

- VS

  ![image-20210701110116765](/Users/ylf9811/Library/Application Support/typora-user-images/image-20210701110116765.png)

## fastp

- Adapter trimming：SE-adapter sequences are detected by assembling the high-frequency read tails//TODO；PE-adapt- er sequences are detected by finding the overlap of each pair.

  The adapter-sequence detection algorithm is based on two assumptions: the first is that only one adapter exists in the data; the second is that adapter sequences exist only in the read tails.

- Base correction：👌

- Sliding window quality pruning：和👆Automatic trimming类似，加上了滑动窗口的思想

- polyG and polyX tail trimming：//TODO具体咋做

- UMI preprocessing：//TODO a sample index or inserted DNA？论文里大概说了fastp把已有的umi工具的功能合成进来了，但是并没有详细说明功能是什么。

- Output splitting：splitting by file lines and splitting by file numbers//TODO

- Duplication evaluation：

- Overrepresented sequence analysis：



## MultiQC

a tool to create a single report visualising output from multiple tools across many samples, enabling global trends and biases to be quickly identified.

//TODO confounding batch effects ? 

- Leek,J.T. et al. (2010) Tackling the widespread and critical impact of batch ef- fects in high-throughput data. Nat. Rev. Genet., 11, 733739.
- Meyer,C.A. and Liu,X.S. (2014) Identifying and mitigating bias in next- generation sequencing methods for chromatin biology. Nat. Rev. Genet., 15, 709721.
- Taub,M.A. et al. (2010) Overcoming bias and systematic errors in next gener- ation sequencing data. Genome Med., 2, 87

## FQC

👆distinguishes FQC from similar tools (e.g. MultiQC; Ewels et al., 2016) designed to summarize FASTQ data sets as individual or groups of samples, but that lack the ability to display multiple, single-sample reports in a unified dashboard.

## fastqc

https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/

- Basic Statistics

- Per Base Sequence Quality

- Per Sequence Quality Scores

- Per Base Sequence Content

- Per Sequence GC Content

- Per Base N Content

- Sequence Length Distribution

- Duplicate Sequences

- Overrepresented Sequences

- Adapter Content

- Kmer Content

- Per Tile Sequence Quality



Parser 

io -> zip  | 流式处理 ｜

thread ->  threadPool ｜ 

functions

Vi 

## 0706

now cmd parser use CLI11, IO use RabbitIO, 



#### 1.文件命名规则

文件名全部小写，可以含下划线或连字符，按项目约定命名,且尽量保证文件名明确。比如：
`cmd_save_player_info_class.cc ，my_use_full_class.cc`

定义类的文件名一般是成对出现，如：`foo_bar.h foo_bar.cc`

若是类中含大量内联函数，我们还可使用-ini.h文件，使之文件内容更加清晰，于是又如：
`url_table.h url_table.cc url-table-ini.h`

#### 2.类命名规则

类型命名每个单词首字母大写，不含下划线，以名词形式。比如：`MyPalyerManager`，这对于所有类型命名一样，类结构体，枚举，类定义都是如此，比如：`MyExcitingEnum`

#### 3.变量命名规则

变量名一律小写，单词用下划线相连，例如：`int player_id; string table_name;`

特殊的是类成员变量，后跟下划线区别普通变量，比如：`player_name_ player_id_`

全局变量则以 g_ 开头，比如 ：`g_system_time`

当然结构体成员变量还是和普通变量一样,比如：`string name; int num_entries;`

#### 4.常量命名规则

k后面跟大写字母开头的单词，比如：
`const int kDaysInAWeek=7; const string kCompanyName=”Tecent”;`

#### 5.函数命名规则

常规函数每个单词首字母大写，使用命令式语气，比如：`OpenFile() CheckFileName()`，

而存取函数或短小的内联函数使用小写加下划线，且与访问变量相吻合，比如`set_num_errors();`

```cpp
class Player{ 
public: 
void set_player_id(const int player_id){return player_id_=player_id;} 
int get_player_id() const{return player-id_;} 
private: 
int palyer_id_; 
};
```

#### 7.枚举命名规则

枚举类名属于类型名，按类命名，枚举值全大写加下划线，比如：`ENUM_NAME` 。

#### 8.宏变量命名规则

如果你一定要用到宏，全大写加下划线，比如：`define PI_ROUND 3.0`。

#### 9.格式美化

可以借助工具进行美化。方便快捷。比如说我用的Qt里面的Beautifier，就可以进行一键格式化代码。

https://www.jianshu.com/p/f56383486520

#### 8.include规范

c库、c++库、其他库、本地文件





## 0707



根据那天和大师兄的讨论结果，基本确定每个功能一个类，暂时这几这几个类：

- Filter
- Trimer：包含front tail trim和adapter trim
- Umi
- PolyX
- Overrepresented
- Duplicate

| Version              | Se    | Pe   |
| -------------------- | ----- | ---- |
| count lines          | 5.60  |      |
| count bases thread 1 | 12.66 |      |
| count bases thread 4 | 3.42  |      |



👆是简单的se数据的简单信息统计，👇加一点fliter，方便做输出。

afterQC中是这样做的：

```
quality filters->count the number of low quality bases or N, calculate the mean quality of each read, and then determine whether to keep or discard this read.
```

fastp中是：

```
static const int PASS_FILTER = 0;
static const int FAIL_POLY_X = 4;
static const int FAIL_OVERLAP = 8;
static const int FAIL_N_BASE = 12;
static const int FAIL_LENGTH = 16;
static const int FAIL_TOO_int64_t = 17;
static const int FAIL_QUALITY = 20;
static const int FAIL_COMPLEXITY = 24;
```

这里先实现其中的 0 12 16 17 20

|                            | Se    |      |
| -------------------------- | ----- | ---- |
| add sample filter thread 1 | 13.69 |      |
| add sample filter thread 4 | 3.67  |      |
|                            |       |      |

有了简单的过滤之后就有输出过滤后read的必要的，下面先简单实现一版output。

才想neoReference的空间使用不会太多，所以直接存应该内存也是够的，不用分批处理。

这一版写的似乎并不巧妙，在统计信息和过滤的同时对pass_data进行拷贝，拷贝到连续的内存中，每64M做成一个string，然后用无锁队列维护，与此同时开一个写线程检测队列是否为空并进行输出。

## 0708

|                                                              | Se    |      |
| ------------------------------------------------------------ | ----- | ---- |
| add simple output thread 1（concurrentqueue.h）              | 38.42 |      |
| add simple output thread 4（concurrentqueue.h）              | 11.94 |      |
| adjust output block size and optimize queue(reserve) thread 1 | 36.74 |      |
| adjust output block size and optimize queue(reserve) thread 4 | 10.84 |      |
| adjust output block size and optimize queue(reserve) thread 1 just no write | 36.79 |      |
| adjust output block size and optimize queue(reserve) thread 4 just no write | 9.95  |      |

现在单线程慢是一次多余的拷贝，多线程加速比一般大概率是因为无锁队列，可以考虑换成原子操作。

👆考虑到队列操作并不多(fileSize/4M)，问题不大，重点还是优化那一次拷贝。

|                                             | Se    |      |
| ------------------------------------------- | ----- | ---- |
| One less memory copy thread 1               | 15.13 |      |
| One less memory copy thread 4               | 15.33 |      |
| One less memory copy thread 1 just no write | 15.03 |      |
| One less memory copy thread 4 just no write | 3.85  |      |
|                                             |       |      |

基本符合预期，减少拷贝之后快了一倍左右，但是多线程的时候卡在写的过程，把write注释就加速比很好了。

下面实现其他部分的信息统计功能，暂时先写一个类似fastqc的版本：

- Basic Statistics
  - filename
  - file type
  - reads number
  - read length
  - GC%
- Per Base Sequence Quality
  - 每个位置的平均质量分：位置-平均质量
- Per Sequence Quality Scores
  - 平均值质量分个数：read平均质量分-read条数
- Per Base Sequence Content
  - 碱基类型占比随位置分布图：位置-碱基占比
- Per Base GC Content
  - GC占比随位置分布图：位置-GC占比
- Per Sequence GC Content
  - 
- Per Base N Content
  - 
- Sequence Length Distribution
- Duplicate Sequences
- Overrepresented Sequences
- Adapter Content
- Kmer Content
- Per Tile Sequence Quality



## 0709

把简单的stateInfo从seqc中拿了出来，重写了一个完整的State类作为信息统计的功能模块，然后添加了几个统计功能，现在基本的统计功能大概都有了。

|                                                     | Se    |      |
| --------------------------------------------------- | ----- | ---- |
| Add some statistics for draw pic thread 1 no output | 23.49 |      |
| Add some statistics for draw pic thread 4 no output | 6.05  |      |
|                                                     |       |      |

## 0710

今天把上面功能们的对应Pe版本写了。

|                                                       | Se    | Pe    |
| ----------------------------------------------------- | ----- | ----- |
| Add some statistics for draw pic thread 1 no output   | 23.49 | 29.45 |
| Add some statistics for draw pic thread 4 no output   | 6.05  | 8.06  |
| Add some statistics for draw pic thread 1 with output | 25.02 | 31.46 |
| Add some statistics for draw pic thread 4 with output | 8.26  | 11.81 |
|                                                       |       |       |

一开始的版本有个bug是双端数据分开考虑质量分，然后分开过滤，这样可能导致过滤之后的p1.fq和p2.fq条数不一样，这肯定不合理。解决方法就是当r1 r2都pass filter的时候才输出。

## 0711

今天写一下adater的检测和cut。

对于功能是默认开启还是关闭的，采取的策略是默认开启，cmd中的trim_adapter_是总的控制开关，默认是打开的，对于单端的数据，只要trim_adapter_是true，就是在开始自动检测adapter，然后处理的时候进行trim；对于双端的数据，默认是不做自动检测adapter的，只要trim_adapter打开了就会做AnalyzeOverlap，依据这个的结果进行接头的去除，如果这个过程失败了，会再进行类似于se的过程。

首先是双端的数据，按照afterQC和fastp中的假设，双端数据除了中间overlap的地方，其他都是adapter，暂时订下只写一个类Adapter，里面有计算重复部分的函数，与此同时依据计算结果，把adapter找出来并trim。

|                                                              |      | Pe    |
| ------------------------------------------------------------ | ---- | ----- |
| Add Pe adapter trim by overlap analyze information and correction of data thread 1 no output |      | 40.67 |
| Add Pe adapter trim by overlap analyze information and correction of data thread 4 no output |      | 11.55 |
| Add Pe adapter trim by overlap analyze information and correction of data thread 4 with output |      | 42.60 |
| Add Pe adapter trim by overlap analyze information and correction of data thread 4 with output |      | 14.82 |

现在check正确性就是简单的和fastp的结果做比较，包括trim adapter之后的read数目，q20bases q30bases，输出文件大小等。

//TOOD 更完备的check正确性，add filter result too some data struct that can do some report.

## 0712

上午先把se的auto-detect-adapter弄好，两种模式，一是自己指定adapter，二是自动检测adapter。前者比较容易实现，只需要做比较简单的寻找和修剪就行了（这里可以采用ktrim的思路进行加速）；后者暂时想到的思路只有使用fastp的字典树进行统计。

上午遇到了一些问题，在测试correct函数的时候发现之前优化版本的diff统计值有点问题，和fastp的输出结果有出入，暂时回退会之前没有优化的版本。

淦！一个地方p2写成了p1找了2小时！

|                                                              |      | Pe    |
| ------------------------------------------------------------ | ---- | ----- |
| Add Pe adapter trim by overlap analyze information and correction of data thread 1 no output |      | 45.33 |
| Add Pe adapter trim by overlap analyze information and correction of data thread 4 no output |      | 12.38 |
| Add Pe adapter trim by overlap analyze information and correction of data thread 1 no output -c |      | 45.39 |
| Add Pe adapter trim by overlap analyze information and correction of data thread 4 no output -c |      | 12.31 |
| Add Pe adapter trim by overlap analyze information and correction of data thread 1 with output -c |      | 48.46 |
| Add Pe adapter trim by overlap analyze information and correction of data thread 4 with output -c |      | 14.73 |

先实现简单的版本--输入adapter只负责修剪掉

## 0713

👆

对于se数据基本的选项：不检测不修剪adapter、自动检测adapter并修剪、输入adapter并修剪

对于pe数据基本的选项：不检测不修剪adapter、自动检测adapter并修剪、输入adapter并修剪（r1 r2都要输入，如果只输入了一个就赋值给另一个）、根据overlap寻找adapter并修剪

cmdInfo中的trim_adapter就相当于总的开关，只要这个是false，所有关于adapter的操作都不进行了

se_auto_detect_adapter_是se自动检测adapter的开关，默认值是true，如果输入了adapter，这个值就置成false

pe_auto_detect_adapter_是pe自动检测的开关，默认是false，因为pe数据默认是利用overlap的数据进行adapter的处理

main函数中基本的逻辑判断大体是这样的：

|                                                            | SeAdapter | Pe    |
| ---------------------------------------------------------- | --------- | ----- |
| Add adapter trim for giving sequence no output -c thread 1 | 22.43     | 62.27 |
| Add adapter trim for giving sequence no output -c thread 4 | 6.00      | 18.54 |

关于自动检测adapter的版本，感觉fastp的实现不是很好，有机会和凯子哥讨论一下再写这部分。

除了对adapter的trim还有质量trim，包括两个方面，一是直接输入front和tail，二是从5‘或者3’端进行滑动窗口质量trim。这个功能模块基本上参考了fastp的代码，做了一点点的改动，测试如下：

|                                                              | SeAdapter | Pe    |
| ------------------------------------------------------------ | --------- | ----- |
| -5 -3 --trimFront1 3 --trimTail1 5  --adapter_seq1 [--adapter_seq2 -c] -w 1 | 21.03     | 61.10 |
| -5 -3 --trimFront1 3 --trimTail1 5  --adapter_seq1 [--adapter_seq2 -c] -w 4 | 5.89      | 16.26 |



下面实现Deduplicate的模块。

fastqc中的处理过程是取前1000000的read来分析以代表整个文件的情况，并且把长度>75的砍成了50，并且每个read使用hash密钥进行统计，此外对于双端的数据它是分开进行的重复统计，这样会导致结果偏高。

fastp相对于fastqc进行了改进，它统计的整个文件的所有read，对于一条read，前12位（如果有N就不统计这一条）hash成key，后32位hash成kmer32（除去最后5位），两个数组A[1<<24] B[1<<24]，只有当B[key]==kmer32的时候才有A[key]++；对于双端的数据，把r1前12位做key，r2后32位做kmer32，这样可以结合pe数据的特性取得更准确的结果，之前在毕业论文里也统计过这种记录方式的碰撞问题，确实存在问题，但问题不大。此外，fastp中明确说明了se啥的数据重复率可能被高估，pe则把握比较大。//TODO ？ 

暂时还没有读到相关的其他方法的论文，先把fastp的模块加进来。

|                                                              | SeAdapter | Pe    |
| ------------------------------------------------------------ | --------- | ----- |
| -5 -3 --trimFront1 3 --trimTail1 5  --adapter_seq1 [--adapter_seq2 -c] -w 1 | 22.81     | 61.60 |
| -5 -3 --trimFront1 3 --trimTail1 5  --adapter_seq1 [--adapter_seq2 -c] -w 4 | 6.24      | 16.73 |
|                                                              |           |       |

## 0714

上午看了看fastp和fastqc关于Overrepresented Sequences的部分。

Overrepresented Sequences即过度代表的序列，说白了就是把出现频率特别高的序列找出来作报告，他和👆duplicate模块一样只是发现问题，暂时还不能解决问题。

fastqc中的做法是只统计前1000000条read，找到比例超过0.1%的序列然后和常见的污染物列表比对，fastp指出了这种方法存在的问题，并进行了改进：统计前1.5Mbase中出现频率较高的序列，记录到hotSeqs中，然后对整个文件统计hotSeqs中序列的出现次数，依次来统计过度表示。

考虑了一下暂时不写了，fastp实质上还是统计了前1.5Mbase中的序列，我觉得不够合理，虽然统计整个文件的序列信息非常耗时，暂时等一下，下次开会找学长商量一下，等下周回来找找论文。



简单加一下poly模块

淦 fastp有个地方写的可能数组越界，找了半天。

|                                                              | SeAdapter | Pe    |
| ------------------------------------------------------------ | --------- | ----- |
| ./RabbitQCPlus -i $data/SRR2496709_1.fastq -I $data/SRR2496709_2.fastq -5 -3 --trimFront1 3 --trimTail1 5 -w 1 -g -x -c --adapter_seq1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_seq2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT | 23.10     | 62.27 |
| ./RabbitQCPlus -i $data/SRR2496709_1.fastq -I $data/SRR2496709_2.fastq -5 -3 --trimFront1 3 --trimTail1 5 -w 4 -g -x -c --adapter_seq1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_seq2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT | 6.91      | 17.19 |
|                                                              |           |       |



## 0720

沈阳回来了，继续干活

今天简单加了一下umi的功能，//TODO 一部分的参数检测还没有加，比如umiLen不能超过100，可以参考option.cpp。

因为对ref的操作不再仅仅是简单的substr，所以新的name地址可能和原来的并不连续，这里注意及时释放。

windows

|                                                              |      | Pe    |
| ------------------------------------------------------------ | ---- | ----- |
| ./RabbitQCPlus  -i ../../data/SRR2496709_1.fastq -I ../../data/SRR2496709_2.fastq --addUmi --umiLoc per_read -o p1.fq -O p2.fq --umiLen 4  -w 1 |      | 77.27 |
| ./RabbitQCPlus  -i ../../data/SRR2496709_1.fastq -I ../../data/SRR2496709_2.fastq --addUmi --umiLoc per_read -o p1.fq -O p2.fq --umiLen 4  -w 4 |      | 64.52 |
| ./RabbitQCPlus  -i ../../data/SRR2496709_1.fastq -I ../../data/SRR2496709_2.fastq --addUmi --umiLoc per_read --umiLen 4  -w 1 |      | 57.87 |
| ./RabbitQCPlus  -i ../../data/SRR2496709_1.fastq -I ../../data/SRR2496709_2.fastq --addUmi --umiLoc per_read --umiLen 4  -w 4 |      | 26.93 |

//TODO 似乎是个热点，暂时先实现功能，以后去服务器重点测性能

参考fastp（RabbitQC）添加了auto detect adapter模块

|                                                              | Se        | Pe         |
| ------------------------------------------------------------ | --------- | ---------- |
| ./RabbitQCPlus -w 1 -i ../../data/SRR2496709_1.fastq -I ../../data/SRR2496709_2.fastq --decAdaForPe | 28.7(0.8) | 82.14(1.5) |
| ./RabbitQCPlus -w 4 -i ../../data/SRR2496709_1.fastq -I ../../data/SRR2496709_2.fastq --decAdaForPe | 10.5(0.8) | 33.24(1.6) |
|                                                              |           |            |

## 0805

好啊，pac差不多了，来看看QC了。

基本的功能除了Overrepresented都有了，fastp的Overrepresented功能感觉不太好用，有点暴力，运行起来太慢了，准备过几天读读论文写写新方法，功能部分的开发暂时先这样，参考RabbitQCPlus搞一下性能。

## 0807

整理了一下现在有的功能，将分别进行测试：

fat节点，经典数据三件套

cost ： STD/RabbitQCPlus/ReRabbitQC

| version    | detect adapter cost /s | total cost /s     |
| ---------- | ---------------------- | ----------------- |
| init se    | 1.68/1.45/2.44         | 73.58/43.58/27.86 |
| init -o se | 1.68/1.47/2.37         | 76.97/50.12/31.64 |
| init pe    | x                      | 85.53/42.77/47.17 |
| init -o pe | x                      | 92.64/49.84/51.47 |

关于重复度分析的部分加了向量化，效果还可以。

//TODO correct部分把判断条件改了，感觉是没啥问题。

为了保证测试时间的时候Re和Plus执行的功能是一样的，就把Plus的UseLong打开了，并且在Re中添加了kmer统计的功能。

好啊，现在基本上的优化都加上了，除了detect adapter和stateRead的向量化，来测测时间：

| version               | detect adapter cost /s | total cost /s （no detect） |
| --------------------- | ---------------------- | --------------------------- |
| init se thread 1      | 1.67/1.52/0.41         | 71.35/45.58/26.25           |
| init se -o thread 1   |                        | 74.15/53.12/30.96           |
| init se thread 2      |                        | 39.17/24.23/14.95           |
| init se -o thread 2   |                        | 41.95/28.05/16.42           |
| init se thread 4      |                        | 20.34/12.59/7.77            |
| init se -o thread 4   |                        | 21.38/15.03/8.58            |
| init se thread 8      |                        | 11.30/6.94/4.15             |
| init se -o thread 8   |                        | 11.88/8.33/6.04             |
| init se thread 16     |                        | 6.07/3.97/3.54(12)          |
| init se -o thread 16  |                        | 6.95/8.24(10)/6.04(8)       |
| init se thread >16    |                        | 3.75(30)/3.57(24)/3.54(12)  |
| init se -o thread >16 |                        | 6.88(20)/8.24(10)/6.04(8)   |
|                       |                        |                             |
|                       |                        |                             |
| init pe thread 1      |                        | 82.77/45.89/23.10           |
| init pe -o thread 1   |                        | 91.24/51.76/27.28           |
| init pe thread 2      |                        | 42.35/23.29/11.78           |
| init pe -o thread 2   |                        | 47.79/30.05/14.32           |
| init pe thread 4      |                        | 22.07/12.32/6.16            |
| init pe -o thread 4   |                        | 25.09/15.84/7.69            |
| init pe thread 8      |                        | 11.90/6.82/4.58             |
| init pe -o thread 8   |                        | 13.32/9.00/4.62             |
| init pe thread 16     |                        | 6.26/4.43(14)/4.48(10)      |
| init pe -o thread 16  |                        | 7.53/7.08(10)/4.62(8)       |
| init pe thread 32     |                        | 5.08(22)/4.43(14)/4.48(10)  |
| init pe -o thread 32  |                        | 7.08(18)/7.08(10)/4.62(8)   |
|                       |                        |                             |
|                       |                        |                             |

## 0808

可以发现，pe的时候，Re比Plus要慢一点，用vtune看了看，居然是stateInfo那里，可能向量化确实是有用的，下午给Re的那一块写一下，现在睡会觉。

写完向量化确实单线程快了很多，直接起飞，但是峰值性能还是略差，简单分析了一下，猜想是生产者读的慢，拖慢了整体的速度，把消费者的功能都关了，果然pe的数据一个读2s一个2.9s，和测出来的数据基本上是一样的。

用vtune跑一下发现：

![image-20210808233224806](/Users/ylf9811/Library/Application Support/typora-user-images/image-20210808233224806.png)

![image-20210808233615857](/Users/ylf9811/Library/Application Support/typora-user-images/image-20210808233615857.png)

首先countline比较好处理，看看汇编发现是自动向量化的程度不同，至于为什么几乎一样的代码自动向量化的程度不同就先不管了，手写一下之后这不部分的时间差不多了👇：

![image-20210808233231909](/Users/ylf9811/Library/Application Support/typora-user-images/image-20210808233231909.png)

但是很奇怪的多了好几个copy函数，在总的显示的是这个👇：

![image-20210808233920868](/Users/ylf9811/Library/Application Support/typora-user-images/image-20210808233920868.png)

这块时间是哪里来的就比较迷。一开始以为是一个pool的原因，改成两个pool的版本还是一样。

经过一一的排查，最终确定是SwapBufferSize的原因，fastp设置的是1<<13，RabbitIO设置的是1<<20，这个参数可能要根据数据集或者机器配置调整，明天再看详细的，先睡觉。

## 0809

重新测试时间：

更新table👆。

有个小小的小问题，auto detect adapter的时候，比较慢，都没开多线程，而且Re要比STD和Plus还慢一点。原因是开发的时候懒得写，用的是Reference，多了很多拷贝啥的，现在改成neo就好了，顺便加加简单的并行，现在0.4s（Se），没有啥突出的热点，懒得再改了。

更新table👆。

关于overrepresent sequence部分，fastp基本的思路是先拿前10000条read每10、20、50、100、150长度统计子串，再把频率高的子串（之间子串去重）存下来，后面每隔100（默认）条read统计一下5种长度的子串，但是fastp（目前版本）似乎并没有把预处理时的高频子串复制到用到的map里面。

## 0810

突然发现Plus好像用的都是int，现在换成long重新跑一下：

更新table👆。

以及发现内存有的地方没有及时delete，导致内存使用量并不是定值，这里稍微改了一下。

对于overrepresent部分fastp的有点离谱，不好用，先不管了，先看看补充一下功能，然后赶紧画画图。

## 0811

基于zz毕设的画图结构，简单画了几个。

然后把TGS写了，顺便画好。

## 0812

今天找了找overrepresent的论文，根据论文里提到的数据找到一个会出现过度表示的数据。

是什么454巴拉巴拉，直接运行好像画的图很丑，不如fastp好看，以为是代码写的有bug，但找了找发现是fastp在画图的时候做了一点点美化，比如下标很密集的话就把长度100区间的val取平均值之类的。

此外，ReRabbitQC还加了自动的buffer extend等功能。

## 0813

添加了过度表示的功能以及绘图的功能，但是有点bug。仔细看了看也不是bug，就是打印的热点sequence，sampling=1的时候就是对的，但是=20就不对了，是因为多线程采样到的sequence不一样。

## 0814

简单优化了一下过度表示部分，把map换成了unordered_map，这样会出现一点点问题，就是在“去除子串”的过程里由于循环顺序变了导致结果会不一样，简单改了一下这块的逻辑，然后加了一下check函数。

好啊，现在过度表示部分的功能和check都做好了，开始着手优化。

6148(fat)

|                       | Over pre         | Over and other     |
| --------------------- | ---------------- | ------------------ |
| STD                   | 11.9781(10.9106) | 47.1976(43.9154)   |
| Re                    | 10.21341         | 55.28075           |
| Unordered_map         | 11.26912         | 30.11468           |
| All unordered_map     | 6.47571(6.73135) | 30.02323(21.20933) |
| 👇 fat                 | 👇 fat            | 👇 fat              |
|                       |                  |                    |
| STD new data          | 5.89946          | 91.3892            |
| Re unordered new data | 7.02323          | 27.63064           |
|                       |                  |                    |
|                       |                  |                    |

## 0816

换数据了，换机器了。

pre over似乎不是很好改，先不管了，主要是这个state里面的部分，花的时间很多很多，简单查了一下热点，本来以为是dist的统计很慢，想着可以用查分，但实际上是count，也就是说，绝大部分的seq都不在hot_seq里面，所以策略改为优化map。

先简单写个链式hash试试，在此之前先把内层的i和s调换一下顺序试试：

|                        | Over pre | Over and other |
| ---------------------- | -------- | -------------- |
| STD                    | 5.89946  | 91.3892        |
| Re change i s          | 6.76883  | 28.59623       |
| use hash by 👋          | 7.19042  | 18.57580       |
| STD thread 1           | x        | 2239.9         |
| Re thread 1            | x        | 350.865        |
| pre over use hash by 👋 | 0.39397  | 18.25773       |
| 07data STD             | 11.2057  | 43.6952        |
| 07data Re              | 1.48717  | 45.85724       |
| Pedata STD             | 4.20556  | 33.0239        |
| Pedata Re              | 0.56984  | 9.42415        |

现在常用的二代数据效果比较明显，但是下好的那个454的数据效果似乎不好。

好啊，今晚先简单写写论文。

## 0817

对pe数据添加了html report和over represented analyze。

然后实现一下InsertSize Analyze：

```c++
void PairEndProcessor::statInsertSize(Read *r1, Read *r2, OverlapResult &ov) {
    int isize = mOptions->insertSizeMax;
    if (ov.overlapped) {
        if (ov.offset > 0)
            isize = r1->length() + r2->length() - ov.overlap_len;
        else
            isize = ov.overlap_len;
    }

    if (isize > mOptions->insertSizeMax)
        isize = mOptions->insertSizeMax;

    mInsertSizeHist[isize]++;
}
```

这一块是fastp（RabbitQC）中用来统计insertSize的，也就是双端数据实际上的序列长度。

本来是想直接用fastp中的画图模块的，但是出现了bug，暂时写了个简单的折线图的版本，并且不再仅仅统计thread0的信息，是所有线程的都要统计。

## 0818

今天先测试一下RabbitIO的解压缩模块怎么样（only producer）：

|             | STD   | Re    |
| ----------- | ----- | ----- |
| Se data     | 23.08 | 22.72 |
| new Se data | 14.65 | 14.84 |
| Pe data     | 23.52 | 23.09 |

确实慢啊zlib，大师兄的意思好像libdeflate不支持流式处理，寒假搞得bam解压缩那一套又不能用（bam不是标准的gz格式，是分块压缩的），暂时先这样吧。

加一下gz输出部分。

输出的时候压缩也很慢，而且现在没有写 输出队列大小限制 ，也就是说可能内存需求很大很大。

## 0824

加了inteeleaved的input和output，输出没啥问题，输入的话，目前有点小bug，多线程的时候答案不太对，因为某个块是奇数个reads的时候会有问题，凯子哥之前写了这块的代码，下午问问他。

更新interleaved的部分。

## 0901

芜湖九月份了

加加stdin和stdout吧。

加了phred64。

发现了点子东西，最新的fastp里面更换了新的重复率判定方式（新的比旧的判定的重复度更低），并且添加了简单的去重机制。

## 0902

fastp中新的重复度检测模块是

？？

## 0911

怎么上面写了半句就没了。。

今天做pac突然找到了MinIONQC这个软件，他好像不是处理的fastq文件，但是是真对三代测序数据的，可以参考里面的指标丰富一下三代模块。

## 1031

抓住10月份的小已巴！

今天把pugz&pigz整理一下弄到ReRabbitQC里面去。

淦 怎么测的STD好久之前pe过滤之后的reads数目就不对了啊。

```
wa：
Filtering result:
reads passed filter: 23493144
reads failed due to low quality: 1714548
reads failed due to too many N: 7130
reads failed due to too short: 0
reads with adapter trimmed: 1143406
bases trimmed due to adapters: 15890010

ac：
Filtering result:
reads passed filter: 23490880
reads failed due to low quality: 1716812
reads failed due to too many N: 7130
reads failed due to too short: 0
reads with adapter trimmed: 1120462
bases trimmed due to adapters: 15528554
```

哦哦哦 想起来了

哎 还是怪自己这种记录文件写的不够清晰，很久之前为了方便向量化参考了最近版的fastp中的overlap方法，有些边界情况改了一下，具体的就是AdapterTrimmer::trimByOverlapAnalysis中和OverlapAnalysis::analyze中，这两个地方改一下就和RabbitQC一模一样了。

## 1101

qiao!

原来的rabbitio代码，把1<<13改回20就好了，但是和凯子哥分析了一波，感觉diff！=0 的时候往后找不大行，要往前，就把大师兄前段时间弄的rabbitQC里面的操作整一下，结果一直有bug，精准定位了一下，发现是swapbuff的大小也是SwapBufferSize，往前的话就不够大存不下了，修改方案是加了个新的变量GetNxtBuffSize。

## 1102

好啊，现在改了改rabbitQC现在关于se的bug，然后

然后没了

## 1103

今天加一下pxgz。

首先pugz吧，挺好弄的。

## 1104

今天加一下pigz。淦，delete位置不对。。。找了半天。只加了se的。

## 1120
pigz pe的时候出了很多问题，因为两个thread的全局变量都是共享的，只有一份会产生冲突，试了各种方法也不能解决这个问题，就准备先放一放，往后大改的时候顺面把动态create线程加上去。

然后华大那边说测新数据的时候有问题，数据发过来看了看，是结尾的地方有问题，虽然r1 r2的reads数目一样，但是如果大小差的很多，就可能倒数第二次read的时候，r1一下子全部读完了，但是r2还没有，就会导致直接退出，统计的reads数目就会少一点点。

解决方案就是分别开了两个eof，根据自己的eof判断要做什么。

## 1121
now still attemp to add pigz
let's define tnum = pigz thread numbers, for every global variable, we new tnum times, and define tid = {0,1,2...}, when call a function A(), we now call A(tid), when use a globa variable x, now use x[tid]. The only problem is when create pthread, we should padd tid to it.
It seems that this way it's ok, try try try.

## 1123
好啊，终于把pigz pe数据弄好了，但是可能有问题，因为只是把全局变量开了副本，但是两个线程似乎还是共用一套thread的，就比如joinall好像是等待所有的两个线程以及子线程完成才行。
//TODO replace all threads

## 1208
捣鼓了好几天的课设和比对算法，该干干正事了，着手看forg和写paper的书吧，给QC加点三代的玩意，改改画图，读读paper，弄弄fastp里面那个新的duplicate部分。

## 1209
干活干活，先把minion里面关于三代数据的质控研究一下看看能不能加进去。

## 1210 
把TGS的pugz部分加上了

## 1212
时光飞逝啊，功能都差不多了，看看之前大师兄发的综述里面提到的QC的paper们吧，找找思路。

## 1213
昨晚读了读那个综述，这篇paper讨论了ONT三代数据在长度、准确率和吞度量上的提升，以及现有的基于该类数据的各种方法和软件，最后讨论了存在的限制和解决方案。
今天上午简单过了一下NanoPack这个软件，它是python开发的，里面大致有几个脚本，能够对三代数据进行QC，有几个图还是值得参考的，具体的可以看补充材料。
然后SQANTI这个软件，paper太长了，没咋读懂，感觉是和转录组有关的//TODO
LongQC


## 1217
因为要写专利啥的，需要全面测一测性能的提升。

| data type & function & thread num | RabbitQC      | RabbitQCPlus   | fastp | FASTQC |
| --------------------------------- | ------------- | -------------- | ----- | ------ |
| se 7.5G & all & thread 1          | 1.7+83.1      | 0.4+30.7       | 113   | 105    |
| se 7.5G & all & thread 2          | 1.7+43.3      | 0.4+16.7       |       |        |
| se 7.5G & all & thread 4          | 1.7+22.6      | 0.4+9.3        |       |        |
| se 7.5G & all & thread 8          | 1.7+12.1      | 0.4+5.0        |       |        |
| se 7.5G & all & thread 16         | 1.7+6.9       | 0.4+4.5        |       |        |
| se 7.5G & all & thread 32         | 1.7+6.4       |                |       |        |
|                                   |               |                |       |        |
| pe 3.3+3.3G & all & thread 1      | 100.8         | 29.1           | 107   |        |
| pe 3.3+3.3G & all & thread 2      | 51            | 15.2           |       |        |
| pe 3.3+3.3G & all & thread 4      | 25.7          | 8.5            |       |        |
| pe 3.3+3.3G & all & thread 8      | 13.8          | 6.0            |       |        |
| pe 3.3+3.3G & all & thread 16     | 7.9           | 6              |       |        |
| pe 3.3+3.3G & all & thread 32     | 7.6(20)       | 6              |       |        |
|                                   |               |                |       |        |
| se 7.5G & overrep & thread 1      | 2100+80=2180  | 270+30=300     |       |        |
| se 7.5G & overrep & thread 2      | 1100+45=1145  | 138+17=155     |       |        |
| se 7.5G & overrep & thread 4      | 580+24=604    | 68+10=78       |       |        |
| se 7.5G & overrep & thread 8      | 305+14=319    | 37+5=42        |       |        |
| se 7.5G & overrep & thread 16     | 164+8=172     | 19+5=24        |       |        |
| se 7.5G & overrep & thread 32     | 90+8=98       | 10+4.5=14.5    |       |        |
|                                   |               |                |       |        |
| pe 3.3+3.3G & overrep & thread 1  | 1724+100=1824 | 230+29=259     |       |        |
| pe 3.3+3.3G  & overrep & thread 2 | 874+51=925    | 118.5+15=133.5 |       |        |
| pe 3.3+3.3G & overrep & thread 4  | 444+26=470    | 60+9=69        |       |        |
| pe 3.3+3.3G & overrep & thread 8  | 232+14=246    | 30.9+6=36.9    |       |        |
| pe 3.3+3.3G & overrep & thread 16 | 132+8=140     | 14.6+6=20.6    |       |        |
| pe 3.3+3.3G & overrep & thread 32 | 76+8=84       | 7.3+6=13.3     |       |        |
|                                   |               |                |       |        |
|                                   |               |                |       |        |
| se 7.5G & gz all & thread 1       | 100           | 33             | 128   | 126    |
| se 7.5G & gz all & thread 2       | 100           | 18             | 106   | 132    |
| se 7.5G & gz all & thread 4       | 100           | 11             | 108   | 128    |
| se 7.5G & gz all & thread 8       | 100           | 10             | 107   | 129    |
| se 7.5G & gz all & thread 16      | 100           | 10             | 105   | 128    |
| se 7.5G & gz all & thread 32      | 100           | 10             | 110   | 130    |
|                                   |               |                |       |        |
| pe 3.3+3.3G & gz all & thread 1   | 100           | 34             | 115   |        |
| pe 3.3+3.3G & gz all & thread 2   | 50            | 17             |       |        |
| pe 3.3+3.3G & gz all & thread 4   | 45            | 9              |       |        |
| pe 3.3+3.3G & gz all & thread 8   | 45            | 6              |       |        |
| pe 3.3+3.3G & gz all & thread 16  | 45            | 6              |       |        |
| pe 3.3+3.3G & gz all & thread 32  | 45            | 6              |       |        |



## 0109
考完试了，生病也好的差不多了，这两天加一加三代数据的新指标
 - 长度-含量
 - 一些基本的统计量

## 0111
昨天有点摸啊，支楞起来啊。昨天看了longQC的论文，它指出有的三代数据中质量分是不能用的，所以采用的是比对的策略找到 无意义read
今天上午做一下明天小课堂的ppt吧，顺便梳理一下paper的东西。

## 0113
好像0112写了点什么东西，但是在实验室的电脑上，人已经回家了，先不管了。
昨天讲小课堂发现过度表示模块似乎可以用布隆过滤器，应该是有效果的。



## 0115

回家玩的差不多了，今晚写写布隆过滤器试试，由于vpn受不太了就在本地跑一跑。

先测测本地现在版本的速度，以及抄一抄原来fat节点的速度：

|                               | fat RQC      | fat RQCP    | mac RQCP | mac RQCPP |
| ----------------------------- | ------------ | ----------- | -------- | --------- |
| se 7.5G & overrep & thread 1  | 2100+80=2180 | 270+30=300  | 305      | 287       |
| se 7.5G & overrep & thread 2  | 1100+45=1145 | 138+17=155  |          |           |
| se 7.5G & overrep & thread 4  | 580+24=604   | 68+10=78    | 91       | 83        |
| se 7.5G & overrep & thread 8  | 305+14=319   | 37+5=42     |          |           |
| se 7.5G & overrep & thread 16 | 164+8=172    | 19+5=24     |          |           |
| se 7.5G & overrep & thread 32 | 90+8=98      | 10+4.5=14.5 |          |           |
|                               |              |             |          |           |
|                               |              |             |          |           |
|                               |              |             |          |           |

## 0117

今天开始写写paper了，看了看之前写的一点点，有点啰嗦，而且新加的功能体现的不明显。

RQCP的几个亮点工作是：单线程性能（普通PC也能达到系统性能峰值）、过度表示、压缩文件。其中单线程性能是重点，采用了向量化、全新的存储技术等，快了3倍左右；过度表示使用了全新的数据结构存储和布隆过滤器，快了接近一个数量级；压缩文件部分用了pigz和pugz，多线程能快一个数量级，但是pugz不稳定（不过两线程从来没出错过），pigz双端数据目前该的不够优雅，可能有bug，并且似乎有的时候pigz压缩的文件pugz解压会出错。

关于bioinformatics怎么写准备看几篇读过的paper，下面详细列一下每篇的结构：

#### RabbitQC：

- Abstract：强调质控的重要性，指出现有软件不够快--引出RQC软件，支持多种测序数据，快1-2个数量级

- Introduction：介绍了质控的重要性，指出了现有软件的问题——多核利用率低，且不支持三代数据，NanoQC又太慢。提出RQC，并实际列了一下加速效果。

- Methods：先介绍了一下整体的框架，然后详细介绍了框架中提升性能的每个关键点；在功能上，二代数据借用了fastp的功能，三代数据用c++实现了NanoQC中的功能。

- Results：一个表：二代se、二代pe、三代PacBio、三代MinION数据上和现有软件（带版本号）的性能比较，介绍了测试环境、测试数据（有一个大于SSD大小的）、软件参数、输出结果、加速比分析。

#### fastp：

- Abstract：介绍了现在有好多指控相关的软件提供了独立的功能，但是用每个软件都读写一遍数据太慢了，就把所有的功能合成到了一个总的软件上，并用c++重写了，有了性能提升。

- Introduction：介绍了质控的内容和重要性，指出了现有软件的问题，提出了fastp，简单列举了一下他的优点。

- Methods：针对每个模块的优化做了介绍。

- Results：从性能和功能两个方面展示了fastp的创新之处。


#### BGSA：

- Abstract：提出比对算法的重要用途，提出BGSA，列出它相比现有软件在哪些功能上取得了多少倍的加速。
- Introduction：介绍了比对算法的背景和重要性，列举了现有软件支持的各类比对算法；提出了BGSA软件，简单介绍他是个啥
- Materials and methods：介绍了BGSA在功能上和实现上的创新。
- Results：BGSA对比现有软件在各种比对算法上的加速效果。

#### ktrim：

- Abstract：提出去除适配器的意义，指出现有软件不够快；提出ktrim，简单列举了它的特点，性能提升。
- Introduction：介绍了适配器以及去除他的重要性，指出现有的软件不够快，提出了ktrim，它在功能上支持多种数据多种功能，性能上快了2-18倍。
- Implementation：介绍了krim的工作原理，在模拟数据和真实数据上都很准确。
- Results：krim和现有软件在功能、性能、准确率上的比较，做了总结。



## 0118

晚上动笔写了写，感觉还差点东西

- [ ] 过度表示模块的生物学意义
- [x] 怎么解释和RQC、fastp的关系（功能、代码、重构）
- [ ] 一些软件的现在版本性能情况，NanoQC、AfterQC、FASTQC、fastp等
- [x] 重写ORP模块
- [ ] 调查其他软件ORP模块的速度
- [x] RabbitQC和fastp中ORP的 i+=step 这个？？？？
- [ ] 软件名在论文中怎么表示，在句首要大写吗 fastp
- [x] 简单优化一下生产者，能不能多线程峰值性能快一倍
- [ ] 论文中写几倍几倍？
- [x] 功能上的优化怎么说？
- [x] 布隆过滤器要加到论文里吗？
- [x] 优化生产者？？
- [ ] 性能提升是在methods里面说一下，还是不提，全都放在results里面？
- [ ] results里面加具体的实验数据是啥
- [x] 软件名在句首写？
- [ ] gz部分要不要加上 用pugz解压到内存（硬盘上）再读写这种策略
- [ ] 加上⬆️

## 0119

新写的ORP模块还是去fat节点测一测吧

|                               | fat RQC      | fat RQCP    | mac RQCP | mac RQCPP | fat RQCPP   |
| ----------------------------- | ------------ | ----------- | -------- | --------- | ----------- |
| se 7.5G & overrep & thread 1  | 2100+80=2180 | 270+30=300  | 305      | 287       | 267+30      |
| se 7.5G & overrep & thread 2  | 1100+45=1145 | 138+17=155  |          |           |             |
| se 7.5G & overrep & thread 4  | 580+24=604   | 68+10=78    | 91       | 83        |             |
| se 7.5G & overrep & thread 8  | 305+14=319   | 37+5=42     |          |           |             |
| se 7.5G & overrep & thread 16 | 164+8=172    | 19+5=24     |          |           | 17.5+5=22.5 |
| se 7.5G & overrep & thread 32 | 90+8=98      | 10+4.5=14.5 |          |           | 9+4.5=13.5  |
|                               |              |             |          |           |             |
|                               |              |             |          |           |             |

上午测ORP的时候随手跑了一下vtune，发现慢在hash计算上，实际上的访存并不慢，所以加了布隆过滤器效果提升不大。优化hash计算的话，有两个路子，内层和外层。

这里先回忆了一下字符串hash的一些知识https://oi-wiki.org/string/hash/

原来的hash方式用的是5进制64位整数自然溢出，但是似乎这种未定义行为不同编译器结果不一样，规范一点禁止这种未定义行为，换成正儿八经的多项式字符串hash吧。但是但是%对M取模的话太慢了，还是自然溢出了。

但是有个问题就是，字符串hash的话错误率是有点高的，假设两个长度相等的串A B，h(A)=sigma(A[i] * x^i % M)，h(B)=sigma(B[i] * x^i % M)，有f(x)=sigma((B[i]-A[i]) * x^i % M)，碰撞就是fx有根，一元n次多项式不超过n个根，假设%M是个质数，在0～M-1里面均匀分布，那也有t=n/M的概率碰撞，这里假设M是1e18，n是reads长度算是1e2，t=1e-16。但如果有T次(7.5G 100长度的se文件大约有1e4个hotseqs和1e9个子串询问，T=1e13)，那碰撞的概率是1-(1-t)^T，这个数大约是0.0011，即千分之一的概率碰撞，感觉是不怎么靠谱的，这里测试了一下M取1e9，计数就错的很离谱了。但是测试的数据里面计数都没有出现问题的，即使出现问题，对统计高频字串也没啥影响。必要时可以采用双hash的方式。

测试的时候发现了hash上有些数据比较奇怪，debug之后发现是AGCTN对应0-4，AAA和A的hash值是一样的，所以出错了。但是原来的base2val确实是这样的，这里只改了hash部分用到的。

对于外层，相当于10、20、40、100大小的窗口滑动，能不能用类似kmer的方式只弄一遍拿到这些窗口的hash值呢？不好办，因为100长度的窗口用kmer的话相当于只看了最后几个碱基。可以试一试多项式hash询问子串的方式，能降低一点复杂度。因为自然溢出，所以多项式hash字串搞不大行，kmer的话，不行不行。

对于内层，思路就是向量化，但是有着明显依赖所以还是要借助外层实现向量化，这个方法听大师兄说过但是没有试过，应该也不难写，关键是写之前分析好究竟有没有用，别再呼啦啦写完没啥用。

## 0120

先改改论文吧，感觉外层向量化不那么好写哦。

methods部分想分三部分写：基本质控功能优化、过度表示模块优化、读写压缩文件优化。注意说明过度表示没有算在基本的质控功能里面。对于基本质控功能的优化，主要写以下两个部分：可拓展的向量化技术、内存拷贝优化。

今晚突然发现好像开启向量化的话多线程会降频，第一次遇到这种情况比较好奇，找个检测主频的api看看。开开检测所有read的ORP功能，AVX512和自动向量化在32线程下主频是2700和3400，再试试AVX2，但是好像AVX2有的部分功能没写完。。。今晚加加班写了他。算了太多了，暂时用自动向量化的版本放在AVX2上吧。

|                           | no vec    | vec256    | vec512    |      |
| ------------------------- | --------- | --------- | --------- | ---- |
| SRR25 thread 1 base func  | 0.6+36.2  | 0.6+37.6  | 0.6+28.1  |      |
| SRR25 thread 16 base func | 0.6+3.9   | 0.6+3.9   | 0.6+3.9   |      |
| SRR25 thread 32 base func | 0.6+3.9   | 0.6+3.9   | 0.6+3.9   |      |
| SRR25 thread 1 ORP func   |           |           |           |      |
| SRR25 thread 4 ORP func   |           |           |           |      |
| SRR25 thread 8 ORP func   |           |           |           |      |
| SRR25 thread 16 ORP func  |           |           |           |      |
| SRR25 thread 32 ORP func  | 1.4+13.8  | 1.3+14.1  | 1.4+16.1  |      |
|                           |           |           |           |      |
| SRR241 thead 1 base func  | 0.3+27.7  | 0.3+25.0  | 0.3+15.2  |      |
| SRR241 thead 16 base func | 0.3+2.1   | 0.3+1.95  | 0.3+1.8   |      |
| SRR241 thead 32 base func | 0.3+1.7   | 0.3+1.7   | 0.3+1.63  |      |
| SRR241 thead 1 ORP func   | 1.0+172.6 | 1.0+169.4 | 1.0+168.8 |      |
| SRR241 thead 4ORP func    | 1.0+44.3  | 1.0+43.5  | 1.0+43.4  |      |
| SRR241 thead 8 ORP func   | 1.0+23.5  | 1.0+23.2  | 1.0+22.9  |      |
| SRR241 thead 16 ORP func  | 1.0+12.2  | 1.0+12.1  | 1.0+12.3  |      |
| SRR241 thead 32 ORP func  | 1.0+6.7   | 1.0+6.8   | 0.8+7.4   |      |
|                           |           |           |           |      |

记录了一下具体的时间，发现简单的指控功能里面Vec512多线程并不慢。。。。但是ORP的时候就慢了。。。。但是ORP的部分我没写向量化啊？？？明天跑跑vtune看看问题在哪，睡觉睡觉。。 

## 0121 

改专利了，难。

## 0122

先解决前天遗留问题，为什么多线程ORP的时候Vec512反而慢了一些，通过跑vtune发现，vec512确实能加速stateinfo的过程，能快一倍多，但是stateORP函数和HashQueryAndAdd函数却慢了一些。有点子奇怪的是vtune里面stateORP函数都看不了热点分析，只有一个ifxxx，试着把ifxxx放在函数外面试试。

能看到时间了，热点语句的汇编都是一样的，应该就是因为降频了，大概是2600和3100的区别，简单的基本质控功能看不出来可能是因为时间太短了，没有降频或者影响不明显。

好啊，现在这个问题解决了，尝试一下那天大师兄说的两种方式，先来看看nthash。

简单看了看nthash的paper和github，感觉就是多项式hash的翻版，但是适用于基因数据，并且在hash分布和碰撞上都有很好的实验结果，此外还支持一个kmer产生多个hash值，这个暂时不知道咋用。

明天简单试试RQCP中换上这个借口看看效果。

## 0123

换了nthash，快是快了不少，但是答案不太对，数目对不上，查到的cnt也不一样，并且查询次数都不一样，那肯定是写的有问题咯。

一开始发现是原来fastp写的有点奇怪，可能是个bug，i+step<slen应该是s+step<=slen，导致nthash的查询次数比原来的多很多。这个问题修改之后nthash又比原来的版本少了一点点。

找到bug了，怪自己读论文的时候不够仔细，nthash只能处理包含AGCT的基因数据，当kmer遇到N的时候就会返回false，迭代器会忽略false，就会出现少一点点查询的情况。

针对N的处理是在insert或者枚举kmer的时候只要包含N就按照原来的策略算hash值。

并且nthash对s和他的反向互补链的到的hash值是相通的，这个应该可以调参数修改的按理说，大概估计约摸着把这个改了应该就答案正确了。

感觉NTF这个函数就是我们需要的，最起码它会区分s和它的反向互补链，感觉他对N好像也没啥要求，明天换成这个试试。

等不及了，刚才在test.cpp里面试了试，果然，和猜想的是一样的，而且这github上的例子有错误，参数列表和代码里面不一样，一开始用hash老是乱的。

马上12点了，睡觉睡觉，明天再弄。

## 0124

今天弄了弄，kmer长度%4==0的时候NTF64函数老是出错？？？？和NTC、NTS算的hash值不一样， 读了读源码，NTF算第一个kmer似乎是优化过的，但是好像不对？？？？还是我使用的姿势不对？？？

直接用了最简单的循环算第一个kmer（参考的NTS中的），然后问题就都解决了，单线程快了一倍多，多线程快了接近一倍，但是目前NNNNN这种子串的hash值似乎有点子奇怪，回家再研究研究，争取今晚收尾，明天开始写压缩部分的论文了。

似乎全N的子串nthash给出的hash值全是0，这个暂时的解决方案是在把hash值和序列长度都存在hash链表中。实际上没有加slen，用了seq.length()。

下面是新的数据测试

|                           | no vec（no bf）            | vec256   | vec512（no bf）            |
| ------------------------- | -------------------------- | -------- | -------------------------- |
| SRR25 thread 1 base func  | 0.6+36.2                   | 0.6+37.6 | 0.6+28.1                   |
| SRR25 thread 4 base func  | 0.6+10.2                   |          | 0.6+8.0                    |
| SRR25 thread 8 base func  | 0.6+5.5                    |          | 0.6+4.4                    |
| SRR25 thread 16 base func | 0.6+3.9                    | 0.6+3.9  | 0.6+3.9                    |
| SRR25 thread 32 base func | 0.6+3.9                    | 0.6+3.9  | 0.6+3.9                    |
|                           |                            |          |                            |
| SRR25 thread 1 ORP func   | 1.4+36.2+101.3             |          | 1.4+28+108                 |
| SRR25 thread 4 ORP func   | 1.4+10.2+25.8              |          | 1.4+8+27.4                 |
| SRR25 thread 8 ORP func   | 1.4+5.5+13.8               |          | 1.4+4.4+14.7               |
| SRR25 thread 16 ORP func  | 1.4+3.9+6.5                |          | 1.4+3.9+6.6                |
| SRR25 thread 32 ORP func  | 1.4+3.9+2.6（1.4+3.9+3.1） |          | 1.4+3.9+3.1（1.4+3.9+3.7） |

AVX512带来的降频还是有影响。

然后测了测，bf的效果就比较明显了，现在计算和访存的耗时基本差不多了，计算还是多一点，这里计算还是可以向量化应该，暂时就不做了，先写paper。

## 0125

好啊，先去paper那里改改ORP的部分。

简单改改生产者。

|                         | se SRR25（small） | pe SRR24（512）（512+small） | se SRR25 gz in | se SRR24 gz in（avx512） |      |
| ----------------------- | ----------------- | ---------------------------- | -------------- | ------------------------ | ---- |
| thread 1 only producer  | 2.0（1.52）       | 3.3（3.0）（2.6）            |                |                          |      |
| thread 2 only producer  | 2.0（1.52）       | 3.3（3.0）（2.6）            |                |                          |      |
| thread 4 only producer  | 2.0（1.52）       | 3.3（3.0）（2.6）            |                |                          |      |
| thread 8 only producer  | 2.0（1.52）       | 3.3（3.0）（2.6）            |                |                          |      |
| thread 16 only producer | 2.0（1.52）       | 3.3（3.0）（2.6）            |                |                          |      |
|                         |                   |                              |                |                          |      |
|                         |                   |                              |                |                          |      |

## 0126

似乎啊，这个之前的rabbitQC（大师兄更新过的）峰值性能还能再高一点，之前测的不太准，先重新测一下时间吧：

（）表示带输出的性能

| data type & function & thread num | RabbitQC（2022.1.26version） | RabbitQCPlus  | fastp | FASTQC |
| --------------------------------- | ---------------------------- | ------------- | ----- | ------ |
| se 7.5G & all & thread 1          | 0.6+74.7                     | 0.6+27.7      | 113   | 105    |
| se 7.5G & all & thread 2          | 0.6+39.2                     | 0.6+15.5      |       |        |
| se 7.5G & all & thread 4          | 0.6+20.7                     | 0.6+8.0       |       |        |
| se 7.5G & all & thread 8          | 0.6+11.3                     | 0.6+4.3       |       |        |
| se 7.5G & all & thread 16         | 0.6+6.2                      | 0.6+3.1       |       |        |
| se 7.5G & all & thread 32         | 0.6+4.0                      | 0.6+2.7（28） |       |        |
| se 7.5G & all & thread 48         |                              |               |       |        |
|                                   |                              |               |       |        |
| pe 3.3+3.3G & all & thread 1      | 88.7                         | 23.8          | 107   |        |
| pe 3.3+3.3G & all & thread 2      | 44.1                         | 12.7          |       |        |
| pe 3.3+3.3G & all & thread 4      | 22.8                         | 6.6           |       |        |
| pe 3.3+3.3G & all & thread 8      | 12.2                         | 4.0           |       |        |
| pe 3.3+3.3G & all & thread 16     | 7.0                          |               |       |        |
| pe 3.3+3.3G & all & thread 32     | 6.4（28）                    |               |       |        |
| pe 3.3+3.3G & all & thread 48     |                              |               |       |        |

这个版本是生产者1<<12偏移量，好像改变block的大小还能更快，但是要事先检测reads的长度是不是正常的二代reads，如果是的话，就可以把偏移量和块大小弄小点，这样应该是最快的版本。困了，明天写吧。

此外，最新的RQC版本中的init adapter过程似乎也挺快的，RQCP中的还需要优化。





## 0128

上午先把论文写写吧，ORP部分写完。

准备写gz的部分了。

## 1029

简单写了写最后的部分。

## 0205

过完年啦，开始干活咯，键盘好像不知道里面掉出啥东西来了，嘶。

## 0206

今天用grammarly检查一下语法，直接在线改了。

关于画图还是画表？

展示三块性能，单线程基础质控功能（se pe）、过度表示模块、读写gz文件

![image-20220206115647544](/Users/ylf9811/Library/Application Support/typora-user-images/image-20220206115647544.png)

图的话大概这样，分图1abc，表的话

好啊，简单和大师兄聊了聊，fastp就都用小写，2 times faster than xx 这种写法没问题，另外pugz解压在RQCP处理的方式可以写到paper里面，具体解压到内存和硬盘都写上。

下面测测pugz加RQCP协同使用的效果和直接RQCP的pugz模块

| data/pugz/RQCP/pigz/1:hdd 2:mem | pugz+RQCP | RQCP |
| ------------------------------- | --------- | ---- |
| se/4/8/16/1                     |           |      |
| se/4/8/16/2                     |           |      |
|                                 |           |      |
|                                 |           |      |
|                                 |           |      |
|                                 |           |      |
|                                 |           |      |

## 0207

测完了👆的数据。

上午改改论文中大师兄昨天说的

- [ ] 过度表示模块的生物学意义
- [ ] 一些软件的现在版本性能情况，NanoQC、AfterQC、FASTQC、fastp等
- [ ] 调查其他软件ORP模块的速度
- [x] 软件名在论文中怎么表示，在句首要大写吗 fastp
- [x] 论文中写几倍几倍？
- [ ] 把methods中的性能部分挪到results里面
- [ ] results里面加具体的实验数据是啥，补充results的内容
- [x] 软件名在句首写？
- [ ] gz部分要不要加上 用pugz解压到内存（硬盘上）再读写这种策略
- [ ] 加上⬆️
- [x] 加上参考文献
- [ ] 加表格
- [x] add bf filter in paper



好啊，论文简单改改之后把估计mxlen的功能加上去，期间遇到了缺省函数参数和重载函数的一点点问题，补充一下基础知识：

c++中的函数参数允许设置缺省值，但只能是最后几个连续的参数。

重载的定义都很熟了。

晚上改代码的时候对于一个含有缺省参数的函数进行了重载，并且把多加的参数放在了最后，显然调用的时候就会有歧义，导致编译出错。
