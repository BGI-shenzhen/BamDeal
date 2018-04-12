# BamDeal
<b>BamDeal: a comprehensive toolkit for bam manipulation</b>
</br></br> BamDeal is a full-featured toolkit for comprehensive analysis of bam file,

BamDeal is implemented in C/C++ language, available for Linux and Mac OS X operating system. 
###  1) Download and Install
------------
  <b> [Download](https://github.com/BGI-shenzhen/BamDeal/raw/master/BamDeal-0.19.tar.gz) </b>
  </br> </br> Pre-installations of 3 libraries are required before installing BamDeal
  </br> 1 htslib: [samtools-1.6/htslib-1.6](https://sourceforge.net/projects/samtools/files/samtools)）
  </br> 2 g++ with [--std=c++11](https://gcc.gnu.org/) > 4.8  is recommended
  </br> 3 zlib : [zlib](https://zlib.net/) > 1.2.3 is recommended
  </br> 4 R : [R](https://www.r-project.org/) is recommended


- To compile BamDeal, do [ ./configure] first, then [make] 
- Final software can be found in the direcoty [bin/BamDeal]

 For <b>linux /Unix </b> and <b>macOS </b>
<pre>
        tar -zxvf  BamDeal-XXX.tar.gz
        cd BamDeal-XXX;                             
        ./configure	   
        make ;                        
        ./bin/BamDeal                             
</pre>
**Note:**  If fail in [./configure], to find library ***_htslib_**, you can add  LDFLAGS  to find the htslib library  and  CPPFLAGS to find the htslib header [ ./configure LDFLAGS=-L/usr/lib/ CPPFLAGS=-I/usr/include/]


### 2) Features 
------------

### Parameter description</b>
```php

Program: BamDeal
Version: 0.19   hewm2008@gmail.com      Apr 12 2018

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
|       |bamCat         |Merge/Cat diff header muti bam to one bam                         |
|       |bamRand        |random out partly of bam read                                     |
|       |bamShiftQ      |modify seq phred quality in bam                                   |


#### statistics

|Module |    Function   |       Description                                                |
|:-----:|:--------------|:-----------------------------------------------------------------|
|Summary|               |                                                                  |
|       |Coverage       |Calculate Genome Coverage/Depth/GC Dis based Bam                  |
|       |DeteCNV        |detect CNV/Deletion Region by merge Depth info based Bam          |
|       |LowDepth       |GiveOut bed file of low Depth Region(may BigDeletion)             |
|       |               |                                                                  |
|       |               |                                                                  |

#### visualize

|Module |    Function   |       Description                                                |
|:-----:|:--------------|:-----------------------------------------------------------------|
|Summary|               |                                                                  |
|       |StatCheck      |Basic Statistics and qualty control,Show result Fig               |
|       |DepthCov       |Show pdf Fig of Depth Dis & Depth~Coverage                        |
|       |DepthGC        |Show pdf Fig of Depth~GC                                          |
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


### 5）Discussing
------------
- [:email:](https://github.com/BGI-shenzhen/BamDeal) hewm2008@gmail.com / hewm2008@qq.com
- join the<b><i> QQ Group : 125293663</b></i>


######################swimming in the sky and flying in the sea ########################### ##
