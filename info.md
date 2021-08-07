- [x] Filter

- [x] Stats

- [x] Trimer->front tail trim

- [x] Trimer->auto adapter trim

- [x] Trimer->input adapter trim

- [x] Umi

- [x] PolyX

- [ ] Overrepresented

- [x] Duplicate

- [ ] new Duplicate

- [ ] Deduplicate

- [ ] Draw

- [ ] 内存泄漏

- [ ] check 正确性

- [ ] option里面的error detect

  

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
