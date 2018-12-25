#RNA-seq

## PIPELINE

![pipeline.png](images/RNA-SEQ.png)

## 目前基于RNACocktail

目前版本：[RNACocktail_v0.2.2](https://github.com/bioinform/RNACocktail/archive/v0.2.2.tar.gz)

还有一种直接使用docker的，在其他机器上做了测试，用docker的确完全不用自己配置任何东西，就能使用，非常方便。但是docker还未掌握。也想从头看看这个流程，所以还是先自己配置一遍。

## 安装

python版本需求：`python2.7`

python requirements：[requirements](./requirements.txt)

软件需求：[softwares_require](./softwares_require.txt)

```bash
pip install -r requirements.txt
pip install https://github.com/bioinform/RNACocktail/archive/v0.2.2.tar.gz
```

## 注意

安装好后，要不将所有的软件都放到`PATH`目录下，或者修改default.py这个文件。

**查看了源码发现，原来整个流程只支持他的默认参数，也就是比对只能用一个软件。每种都只能用一个软件**

**推荐修改default.py这个文件**

一般情况下，这个文件的目录如下：

`$LIB/lib/python2.7/site-packages/src/default.py`

需要修改的为如下这些(从57行开始)

```python
JAVA_XMS = "-Xms1g" #根据需求进行修改
JAVA_XMG = "-Xmx5g" #根据需求进行修改
JAVA_OPT= "%s %s"%(JAVA_XMS,JAVA_XMG)


HISAT2 = "hisat2"   #修改为绝对路径
HISAT2_SPS = "hisat2_extract_splice_sites.py"   #.py结尾的可以不用修改，是该包中已经安装到$PATH中了
SAMTOOLS = "samtools"   #修改为绝对路径
STRINGTIE = "stringtie" #修改为绝对路径
SALMON = "salmon"   #修改为绝对路径
R_CMD = "R"     #修改为绝对路径
FEATURECOUNTS = "featureCounts" #修改为绝对路径
VELVETG = "velvetg" #修改为绝对路径
VELVETH = "velveth" #修改为绝对路径
OASES = "oases" #修改为绝对路径
LORDEC = "lordec-correct"   #修改为绝对路径
STARLONG = "STARlong"   #修改为绝对路径
SAM2PSL = "sam2psl.py" #修改为绝对路径（和FusionCatcher一个路径）
IDP = "runIDP.py" #修改为绝对路径
IDPFUSION = "runIDP.py" #修改为绝对路径
GMAP="gmap" #修改为绝对路径
STAR_DIR = "/us/local/bin"  #修改为绝对路径
BOWTIE2_DIR = "/us/local/bin"   #修改为绝对路径
PICARD = "picard.jar"   #修改为绝对路径
GATK = "GenomeAnalysisTK.jar"   #修改为绝对路径
JAVA = "java"   #修改为绝对路径
GIREMI = "giremi"   #修改为绝对路径
HTSLIB = ""  #修改为绝对路径
FUSIONCATCHER= "fusioncatcher"  #修改为绝对路径


