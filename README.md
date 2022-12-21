
# BamDeal
<b>BamDeal: a comprehensive toolkit for bam manipulation</b>
</br> BamDeal is a full-featured toolkit for comprehensive analysis of bam file,
BamDeal is implemented in C/C++ language, available for Linux and Mac OS X operating system.
###  1) Download and Install
------------
The <b>new version</b> will be updated and maintained in <b>[hewm2008/BamDeal](https://github.com/hewm2008/BamDeal)</b>, please click below Link to download the latest version
</br><p align="center"><b>[hewm2008/BamDeal](https://github.com/hewm2008/BamDeal)</b></p>

<b> [Download](https://github.com/hewm2008/BamDeal/archive/v0.27.tar.gz) </b>
</br> </br>
For <b>linux/Unix </b> static
</br>you can use the statically compiled programs <i>directly</i>
<pre>
             chmod 755 ./bin/BamDeal_Linux
             ./bin/BamDeal_Linux       # mv BamDeal_Linux    BamDeal
</pre>
</br>

For <b>linux/Unix </b> and <b>macOS </b> compile
  </br> </br> Pre-installations of 4 libraries or softs are required before installing BamDeal
  </br> 1 htslib: [samtools-1.12/htslib-1.12](https://sourceforge.net/projects/samtools/files/samtools)   htslib >=1.12
  </br> 2 g++   : g++ with [--std=c++11](https://gcc.gnu.org/) > 4.8+ is recommended
  </br> 3 zlib  : [zlib](https://zlib.net/) > 1.2.3 is recommended
  </br> 4 R     : [R](https://www.r-project.org/) with [ggplot](http://ggplot.yhathq.com/) is recommended


- To compile BamDeal, do [ ./configure] first, then [make] 
- Final software can be found in the direcoty [bin/BamDeal]
<pre>
        tar -zxvf  BamDeal-XXX.tar.gz
        cd BamDeal-XXX;
        chmod 755  configure;      ./configure
        make;
        mv BamDeal bin/;   ./bin/BamDeal 
</pre>

**Note1:** If fail in [./configure], to find library **_htslib_**, you can add LDFLAGS to find the htslib library and <b>CXXFLAGS</b> and CFLAGS to find the htslib header [ ./configure LDFLAGS=-L/usr/lib/ <b>CXXFLAGS=-I/usr/include/</b> CFLAGS=-I/usr/include/];
</br>  export 'LIBRARY/CFLAGS/CXXFLAGS' environment variables to set the htslib path also can be a way to compile, such : export LD_LIBRARY_PATH=/usr/lib:$LD_LIBRARY_PATH;

**Note2:** you can use the statically compiled programs directly [chmod 755 ./bin/BamDeal_Linux; <b>./bin/BamDeal_Linux</b>]


### 2) Features 
------------

### Parameter description</b>
```php

Program: BamDeal
Version: 0.27   hewm2008@gmail.com      2021-12-02

        Usage:

                convert        convert tools
                modify         modify tools
                statistics     statistics analysis tools
                visualize      visualize tools for bam

                Help           Show help in detail

```

#### convert

|Module |    Function   |       Description                                                |
|:-----:|:--------------|:-----------------------------------------------------------------|
|Summary|               |                                                                  |
|       |soap2bam       |soap    -->  bam/sam Format                                       |
|       |bam2soap       |bam/sam -->  soap    Format                                       |
|       |bam2fq         |bam/sam -->  Fastq   Format                                       |
|       |bam2fa         |bam/sam -->  Fasta   Format                                       |


#### modify

|Module |    Function   |       Description                                                |
|:-----:|:--------------|:-----------------------------------------------------------------|
|Summary|               |                                                                  |
|       |bamFilter      |filter low quality read in bam                                    |
|       |bamSplit       |split single/muti-Bam by chr                                      |
|       |bamAssign      |split single/muti-Bam by assign chr                               |
|       |bamCat         |Merge/Cat diff header muti bam to one bam                         |
|       |bamRand        |random out partly of bam read                                     |
|       |bamSubChr      |get/remove some chr form bam                                      |
|       |bamShiftQ      |modify seq phred quality in bam                                   |
|       |bamLimit       |Limit big bam to muti subbam by fix line                          |

#### statistics

|Module |    Function   |       Description                                                |
|:-----:|:--------------|:-----------------------------------------------------------------|
|Summary|               |                                                                  |
|       |Coverage       |Calculate Genome Coverage/Depth/GC Dis based Bam (v1.32)          |
|       |BasesCount     |Calculate Genome every Site's four base Depth                     |
|       |DeteCNV        |Detect CNV/Deletion Region by merge Depth info based Bam          |
|       |DeteSV         |Detect SV by Pair End Read insert size in Bam                     |
|       |LowDepth       |GiveOut bed file of low Depth Region(may BigDeletion)             |
|       |SiteRead       |Extract the read (with base) covering the SNP site                |

#### visualize

|Module |    Function   |       Description                                                |
|:-----:|:--------------|:-----------------------------------------------------------------|
|Summary|               |                                                                  |
|       |StatQC         |Basic Stat and qualty control,Show gc-depth... result Fig         |
|       |DepthCov       |Show pdf Fig of Depth Dis & Depth~Coverage                        |
|       |DepthGC        |Show pdf Fig of Depth~RefGC                                       |
|       |DepthSlide     |Show Manhattan Fig of Depth sliding Windows along genome          |


### 3) Examples
------------
* 1) convert 
```
    #convert soap 2 bam 
   ./bin/BamDeal   convert   soap2bam   -InSoap <in.soap> -OutBam <out.bam>  -Dict Ref.fa
    #convert bam 2 soap 
    ./BamDeal  convert   bam2soap   -InFile <in.bam>  -OutPut <Out.soap>
```

* 2)  modify 
```
    #  cat diff header bam  to one bam 
    ./BamDeal    modify  bamCat  -InFile A.bam -InFile B.bam  -OutFile C.bam
    #  merge  muti sort bam  to an sort bam 
    ./BamDeal    modify  bamCat  -InList  <bam.Sort.list>  -OutFile <out.sort.bam>  -Merge
    #  split bam by chr 
     ./BamDeal   modify   bamSplit  -InList  <bam.list>  -ReSetHead   -OutDir  ./
```

* 3)  statistics 
```
    # detect CNV/Deletion Region based bam file  
      ./BamDeal   statistics  DeteCNV  -List  <bam.list>  -OutPut  <outPrefix>
    # stat coverage and depth each chr 
      ./BamDeal  statistics  Coverage   -List  <bam.list>  -OutPut  <out>  -Ref  Ref.fa  -Stat
    # GiveOut bed file of low Depth Region 
      ./BamDeal   statistics  LowDepth   -InList     <bam.list>  -OutPut  <out.bed>
```


* 4)  statistics 
```
    #  Show Manhattan Fig of Depth sliding Windows along genome
      ./BamDeal    visualize  DepthSlide   -InList  <bam.list>  -Ref  <Ref.fa> -OutPut  <outPrefix>
    # Show pdf Fig of Depth~RefGC
       ./BamDeal    visualize DepthGC  -InList  <bam.list>  -Ref  <Ref.fa> -OutPut  <outPrefix>
    # Show pdf Fig of Depth Dis & Depth~Coverage 
     ./BamDeal    visualize DepthCov    -InList     <bam.list> -OutPut  <outPrefix>
```

see more other Usage in the Documentation

### 4) Format
------------
Format Introduction

* [sam/bam format](https://samtools.github.io/hts-specs/SAMv1.pdf)
* [soap format](http://soap.genomics.org.cn/soapaligner.html)
     [chinese soap introduction](http://blog.sina.com.cn/s/blog_70b2b6020101b609.html)


### 5) Discussing
------------
- [:email:](https://github.com/hewm2008/BamDeal) hewm2008@gmail.com / hewm2008@qq.com
- join the<b><i> QQ Group : 125293663</b></i>


######################swimming in the sky and flying in the sea ########################### ##

