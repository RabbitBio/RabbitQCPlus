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

