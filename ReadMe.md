# PNGseqR

## Rapid candidate gene Selection and Visualization through multiple Bulk Segregate Analysis with next-generation sequencing

PNGseqR provides a convenient tool to map genetic loci controlling both quantitative and qualitative trait by performing BSA analysis, which integrates multiple algorithms and methods to indicate the magnitude of genome wide signals and to define the candidate regions, and prioritizes the causal genes through DEG and GO analysis. It takes single nucleotide polymorphism (SNP) markers from Next-generation sequencing (NGS) data in Variant Call Format (VCF) format as the input file, and provides four published BSA algorithms to rapidly define the candidate region through permutation test and fractile quantile. Users can choose the appropriate methods according to their data and experimental design. In addition, it also supports performing DEG analysis and GO analysis in order to prioritize the target gene. Once the analysis is completed, the plots can conveniently be exported.

PNGseqR is still under development and is offered with out any guarantee, our article has also been published，user can download through  https://www.mdpi.com/2223-7747/11/14/1821/pdf for more detailed information

PNGseqR was developed in R version `3.5.0` and depends on a number of packages for various aspects of its implementation 

`c("vcfR", "dplyr", "tidyr", "readr", "locfit", "DESeq2", "ggplot2", "topGO", "rtracklayer", "grDevices", "graphics","methods","stats","utils")`


## Table of Contents
  - [Installation](#installation)
  - [Input files](#input-files)
      - [vcf_data.example](#vcf)
      - [DEG_reads.example](#DEG-reads)
      - [GO_Term.example](#GO-term)
  - [Function arguments](#function-arguments)
  - [Use Examples](#use-examples)
      - [Figure 1](#figure-1)
      - [Figure 2](#figure-2)
  - [Supplementary materials](#Pipeline to get genotype file in VCF format)
  
<p>&nbsp;</p>
<p>&nbsp;</p>


 ## Installation
 

PNGseqR can be install using `devtools`, either directly from GitHub,

`devtools::install_github("https://github.com/smilejeen/PNGseqR")`

or by downloading the repository to your computer, unzipping, and installing the `PNGseqR` folder.

`devtools::install("PNGseqR")`

or users can download the archive file from “https://github.com/smilejeen/PNGseqR/raw/main/PNGseqR_0.0.0.9.tar.gz”,

and installed the package by following command:

`install.packages("{Your PATH}/PNGseqR_0.0.0.9.tar.gz", repos = NULL, type = "source")`

**Note:** Apart from regular package dependencies, there are some Bioconductor tools that we use as well, as such you will be prompted to install support for Bioconductor, if you haven’t already. Because of this, in order to install PNGseqR from github you will be required to install some compiling tools (Rtools for Windowsand Mac, respectively).

**If you use PNGseqR in published research, please cite:**

> 

We also recommend citing the paper for the corresponding method you work with.

△SNP method:

> Takagi, H., Abe, A., Yoshida, K., Kosugi, S., Natsume, S., Mitsuoka, C., Uemura, A., Utsushi, H., Tamiru, M., Takuno, S., Innan, H., Cano, L. M., Kamoun, S. and Terauchi, R. (2013), QTL-seq: rapid mapping of quantitative trait loci in rice by whole genome resequencing of DNA from two bulked populations. *Plant J*, 74: 174–183. [doi:10.1111/tpj.12105](https://onlinelibrary.wiley.com/doi/full/10.1111/tpj.12105)

G prime method:

> Magwene PM, Willis JH, Kelly JK (2011) The Statistics of Bulk Segregant Analysis Using Next Generation Sequencing. *PLOS Computational Biology* 7(11): e1002255. [doi.org/10.1371/journal.pcbi.1002255](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002255)

Euclidean Distance method:

> Hill JT, Demarest BL, Bisgrove BW, Gorsi B, Yost HJ (2013) Mmappr: mutation mapping analysis pipeline for pooled rna-seq. Genome Research 23(4):687-697. [doi: 10.1101/gr.146936.112](https://genome.cshlp.org/content/23/4/687)

Empirical Bayesian method:

> Liu, S., Yeh, C. T., Tang, H. M., Dan, N., Schnable PS. (2012), Gene mapping via bulked segregant rna-seq (bsr-seq). Plos One, 7(5), e36406. [doi:10.1371/journal.pone.0036406](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0036406)

## Input files
At a minimum, PNGseqR requires one input data set: 

<ol>
<li>genotype information of SNP Calling in VCF format(produce from GATK https://github.com/broadinstitute/gatk/releases VCF format http://gatkforums.broadinstitute.org/discussion/1268/how-should-i-interpret-vcf-files-produced-by-the-gatk) </li>
</ol>


Three optional data frames may also be supplied:

<ol>
<li>a dataset containing reads count of replicate samples for bulks</li>
<li>go terms corresponding to gene id for analysis species</li> 
<li>a dataset in GTF format contain the position information of genes</li>
</ol>

The formatting parameters of required and all optional input files are summarized below.

#### vcf_data
`vcf_data` is a dataset containing SNP information in  Variant Call Format(VCF) one row per SNP, with eleven columns, containing the following columns:

**Column**|**Description**
-----|:-----
  `CHROM`|The chromosome which SNP fall in. Data type: character
  `POS`|The position on each chromosome of SNP, in base pairs. Data type: character
  `ID`|The name of each SNP. Data type: character  
  `REF`|The base type at this position of reference genome. Data type: character
  `ALT`|Compare to the reference genome,the of alternative base type at this position. Data type: character
  `QUAL`|Quality value of variation site detection. Data type: character 
  `FILTER`|Used to judge if the site passed the filter criteria. Data type: character
  `FORMAT`|Including the genotype information(GT),allele depth in each sample information(AD), total depth of this site information(DP), Genotype quality value(GQ), Likelihood value of genotype(PL). Data type: character  
  `BULK1`|The value cprresponding to FORMAT in BULK1 sample, i.e. "B73_MU". Data type: character 
  `BULK2`|The value cprresponding to FORMAT in BULK2 sample, i.e. "B73_WT". Data type: character 

 ```
> library(vcfR) 
> vcf_data.example <- read.vcfR("snp.vcf",verbose = TRUE) # snp.vcf is a example dataset
> head(vcf_data.example)

		CHROM POS      ID REF ALT QUAL      FILTER	FORMAT           B73_MU                           B73_WT                        
[1,]	"1"   "49527"  NA "A" "G" "428.44"  NA    "GT:AD:DP:GQ:PL" "0/1:8,7:15:99:215,0,222"        "0/1:4,6:10:99:221,0,130"
[2,]	"1"   "51053"  NA "A" "T" "6117.44" NA    "GT:AD:DP:GQ:PL" "0/1:87,161:248:99:4854,0,2383"  "0/1:80,51:131:99:1271,0,2314"
[3,]	"1"   "55373"  NA "A" "T" "7487.44" NA    "GT:AD:DP:GQ:PL" "0/1:108,137:246:99:4585,0,3353" "0/1:75,87:163:99:2910,0,2376"
[4,]	"1"   "139436" NA "G" "A" "637.44"  NA    "GT:AD:DP:GQ:PL" "0/1:11,12:23:99:391,0,362"      "0/1:4,8:12:99:254,0,134"    
[5,]	"1"   "199964" NA "C" "T" "1026.96" NA    "GT:AD:DP:GQ:PL" "1/1:0,15:15:45:625,45,0"        "0/1:1,11:12:0:406,0,0"       
[6,]	"1"   "199995" NA "T" "C" "1120.34" NA    "GT:AD:DP:GQ:PL" "1/1:1,17:18:12:670,12,0"        "1/1:0,11:11:33:456,33,0"     
```

<p>&nbsp;</p>

#### DEG_reads
`DEG_reads` is an optional dataframe containing multiple replicates reads count of each bulk, one row per gene, with at least five columns:

**Column**|**Description**
-----|:-----
  `Gene id`|Gene symbol/name used for differential expression analysis refers to. Data type: character
  `BULK1-1`|Reads number of each gene in the first sample of bulk1 replicate. Data type: integer
  `BULK1-2`|Reads number of each gene in the second sample of bulk1 replicate. Data type: integer
  `BULK1-n`|Reads number of each gene in the nth sample of bulk1*Note: DEG_reads can contain multiple samples such as BULK1-4, BULK1-5 but the number of sample must corresponding with BULK2-n*. Data type: integer 
  `BULK2-1`|Reads number of each gene in the first sample of bulk2 replicate. Data type: integer
  `BULK2-2`|Reads number of each gene in the second sample of bulk2 replicate. Data type: integer
  `BULK2-n`|Reads number of each gene in the nth sample of bulk2*Note: DEG_reads can contain multiple samples such as BULK2-4, BULK2-5 but the number of sample must corresponding with BULK1-n*. Data type: integer 

```
> data("DEG_reads.example")
> tail(DEG_reads.example)
             Gene id  MU1  MU2  WT1  WT2
46113 Zm00001d054106 1606 1810 1318 1431
46114 Zm00001d054107 4746 4490 5527 7359
46115 Zm00001d054109    0    0    0    0
46116 Zm00001d054110   55   46   16    4
46117 Zm00001d054111    6    4   27   10
46118 Zm00001d054112   26    5   29   23
```
  
<p>&nbsp;</p>

#### GO_Term
`GO_Term` is an optional data frame, one row per gene and a go term, which should contain the following columns: 

*Note: The example is a dataset containing maize GO Terms of B73 version4 from agrigo(http://bioinfo.cau.edu.cn/agriGO/index.php)*

**Column**|**Description**
-----|:-----
`Gene id`|Gene symbol/name for gene ontology analysis refers to. Data type: character
`GO Term`|GO term corresponding to the gene. Data type: character

```
> data("GO_Term.example")
> head(GO_Term.example,25)
   Gene id    GO Term
1  Zm00001d012464 GO:0044282
2  Zm00001d012464 GO:0005777
3  Zm00001d012464 GO:0055114
4  Zm00001d012464 GO:0010264
5  Zm00001d012464 GO:0071482
6  Zm00001d012464 GO:0003959
7  Zm00001d012464 GO:0003954
8  Zm00001d012464 GO:0046496
9  Zm00001d012464 GO:0006091
10 Zm00001d012464 GO:0031304
11 Zm00001d012464 GO:0050660
12 Zm00001d012464 GO:0005759
13 Zm00001d012464 GO:0005975
14 Zm00001d010658 GO:0016020
15 Zm00001d010658 GO:0009410
16 Zm00001d010658 GO:0009863
17 Zm00001d010658 GO:0043565
18 Zm00001d010658 GO:0006355
19 Zm00001d010658 GO:0046983
20 Zm00001d010658 GO:0005737
21 Zm00001d010658 GO:0005634
22 Zm00001d010658 GO:0009627
23 Zm00001d010658 GO:0003700
24 Zm00001d010659 GO:0080092
25 Zm00001d010659 GO:0061025
```

<p>&nbsp;</p>

#### gtf_data
`gtf_data` is an optional GRanges object with 6 ranges and 15 metadata columns which should contain the following columns: 

*Note: The example is a dataset containing part of maize genes positional information of B73 version4 from ensembl plant(ftp://ftp.ensemblgenomes.org/pub/plants/release-34/gtf/zea_mays/Zea_mays.AGPv4.34.gtf.gz)*

**Column**|**Description**
-----|:-----
`seqnames`|The ID of the sequence, which can be the ID of chromosome, scaffold or contig.. Data type: factor
`start`|The start position of the feature on the sequence. Data type: integer
`end`|The end position of the feature on the sequence. Data type: integer
`width`|The width of the feature on the sequence. Data type: integer
`source`|The source of this file. Data type: factor
`type`|The feature represented by the region between start and end can be gene, non coding RNA, TE, exon, trascript, etc. Data type: factor
`score`|A floating-point value. When there is a value, it represents the reliability of the above features. Data type: numeric
`phase`|Only valid for comment type "CDS", indicating the position of the starting code. The valid values are 0, 1 and 2. Data type: integer
`gene_id`|The name of this sequence. Data type: character
`gene_biotype`|optional column.The biological type of this sequence such as gene, lincRNA, snoRNA, etc. Data type: character
`transcript_id`|The name of this sequence, if the sequence is a transcript. Data type: character
`transcript_source`|optional column.The source of this transcript. Data type: character
`transcript_biotype`|optional column.The biological type of this sequence. Data type: character
`exon_number`|optional column.The exon number of the transcript. Data type: character
`exon_id`|optional column.The name of exon if above feature is an exon. Data type: character
`protein_id`|optional column.The name of corresponding protein. Data type: character
`protein_version`|optional column.The version of protein function annotation. Data type: character
`gene_name`|optional column.indicate the function of gene. Data type: character

```
> library(rtracklayer)
> gtf_data.example <- import("zma.v4.example.gtf")
> head(gtf_data.example)
GRanges object with 6 ranges and 15 metadata columns:
      seqnames      ranges strand |   source        type     score     phase        gene_id  gene_source   gene_biotype       transcript_id transcript_source transcript_biotype  exon_number                exon_id          protein_id protein_version   gene_name
         <Rle>   <IRanges>  <Rle> | <factor>    <factor> <numeric> <integer>    <character>   <character>    <character>         <character>       <character>        <character>      <character>          <character>         <character>     <character> <character>
  [1]        8 77622-80307      + |  gramene gene               NA      <NA> Zm00001d008171    gramene protein_coding                <NA>              <NA>               <NA>              <NA>                   <NA>                <NA>            <NA>        <NA>
  [2]        8 77622-80307      + |  gramene transcript         NA      <NA> Zm00001d008171    gramene protein_coding Zm00001d008171_T001           gramene     protein_coding              <NA>                   <NA>                <NA>            <NA>        <NA>
  [3]        8 77622-78130      + |  gramene exon               NA      <NA> Zm00001d008171    gramene protein_coding Zm00001d008171_T001           gramene     protein_coding              <NA>                   <NA>                <NA>            <NA>        <NA>
  [4]        8 77951-78130      + |  gramene CDS                NA         0 Zm00001d008171    gramene protein_coding Zm00001d008171_T001           gramene     protein_coding                 1 Zm00001d008171_T001...                <NA>            <NA>        <NA>
  [5]        8 77951-77953      + |  gramene start_codon        NA         0 Zm00001d008171    gramene protein_coding Zm00001d008171_T001           gramene     protein_coding                 1                   <NA> Zm00001d008171_P001               1        <NA>
  [6]        8 78839-79033      + |  gramene exon               NA      <NA> Zm00001d008171    gramene protein_coding Zm00001d008171_T001           gramene     protein_coding                 2 Zm00001d008171_T001...                <NA>            <NA>        <NA>
```

<p>&nbsp;</p>

## Function arguments

**Required Arguments**

Function | Parameter | Description
-----|-----|-----
`vcf2table()`| | The function converts file in VCF format into a data frame
``|file|The name of genotype file in VCF format from GATK SNP calling pipeline
``|BulkA|The sample name of one bulk
``|BulkB|The sample name of the other bulk
``|chromlist|A string vector of the chromosomes to be used in the analysis, used to filter out unneeded contigs etc.
`BSA_filter()`| |The function to filter SNPs based on read depth and quality of each bulk
``|file|The data frame exported by vcf2table()
``|Bulk.DP|The minimum total read depth for a SNP
``|Bulk2.Ref|The minimum read depth of a reference SNP allele 
``|Bulk2.Alt|The minimum read depth of a alternative SNP allele
``|Depth.diff|The maximum absolute differences between the read depths of the bulks
``|min.GQ|The minimum Genotype Quality output by GATK
``|verbose|Default value is “TRUE”. If the value is "TRUE", the number of SNPs filtered in each step will print
`DeltaSNP()`| |The function to perform BSA analysis with △SNP algorithm
``|file|The data frame exported by BSA_filter()
``|Windowsize|A numeric value indicates the window size (in base pairs) containing SNPs subjected to calculating △SNP	
``|threshold|From 0 to 1, default value is “0.95”. A numeric threshold value to extract candidate region by fractile quantile
``|pval|A numeric threshold value. Default value is “0.05”
`Bayesian()`| |The function to perform BSA analysis with empirical bayesian algorithm 
``|file|The data frame exported by BSA_filter()
``|Bulksize|The number of samples in each pool
``|Genomelength|The total length of reference genome
``|Windowsize|A numeric value indicates the window size (in base pairs) containing SNPs subjected to calculating empirical bayesian
``|threshold|From 0 to 1, default value is “0.95”. A numeric threshold value to extract candidate region by fractile quantile
``|pval|A numeric threshold value. Default value is “0.05”
`ED()`| |The function to perform BSA analysis with Euclidean distance algorithm
``|file|The data frame exported by BSA_filter() 
``|power|A numeric value indicates the power to raise ED. It was suggested to be 4~6
``|Windowsize|A numeric value indicates the window size (in base pairs) containing SNPs subjected to calculating Euclidean distance
``|threshold|From 0 to 1, default value is “0.95”. A numeric threshold value to extract candidate region by fractile quantile
``|pval|A numeric threshold value, default value is “0.05”
`Gprime()`| |The function to perform BSA analysis with G-test algorithm	
``|file|The data frame exported by BSA_filter()
``|Windowsize|A numeric value indicates the window size (in base pairs) containing SNPs subjected to calculating G statistic
``|threshold|From 0 to 1, default value is “0.95”. A numeric threshold value to extract candidate region by fractile quantile
``|pval|A numeric threshold value. Default value is “0.05”
`plot_BSA()`| |The function to output BSA analysis result into scatter and line plots
``|file|The data frame exported by the BSA analysis functions 
``|chromlist|A string vector of the chromosome names used to plot line plots specific chromosomes, default value is “NULL”
``|type|A character to indicate the selected algorithm	
``|file.name|A character to name the output plot
``|threshold|From 0 to 1, default value is “NULL”. The red threshold line displayed in the output plot
``|pval|Default value is “NULL”. A numeric threshold value drawn in blue line in the output plot
`DEG_analysis()`| |The function to perform differential expression analysis
``|file|A data frame containing read count of each pool
``|rep.num|A numeric value indicates the number of replications in differential expression analysis
``|sig.level|A numeric value indicates the FDR threshold in differential expression analysis 
``|exp.fold|A numeric value indicates the log2FoldChange threshold in differential expression analysis 
``|plot|Default value is “TRUE”. If "TRUE", the volcano plot will be created
`GO_analysis()`| |The function to perform GO analysis
``|GO_term|A data frame contains at least two columns, one column is genes' name and the other column is the GO term
``|DEG_id|A vector containing the names of the differential expression genes
``|type|A character indicates the mode of GO analysis. There are four options: "CC" indicates cellular component; "BP" indicates biological process, "MF" indicates molecular function, and "all" indicates all three modes. The default is “all”
``|algorithm|A character indicates the algorithm for GO analysis. Three options "classic","weight01", and "elim" indicate three different algorithms
``|statistic|A character indicates two different significance test methods. Two options "fisher" and "ks" indicates two different methods
``|top_num|A numeric value to indicate the number of significant GO terms
``|plot|Default value is “TRUE”. If "TRUE", the acyclic, bubble, and histogram plots will be created
`result_anno()`| |The function to extract genes fall in the candidate region and judge whether gene is differential expression
``|BSA_data|A data frame contains three columns, including the chromosome name, the start and the end positions of the candidate regions
``|gtf_data|The name of reference genome in GTF format
``|DEG_data|A vector containing the names of differential expression genes,default is "NULL"

<p>&nbsp;</p>
<p>&nbsp;</p>

## Use Examples:

This is a basic example which shows you how to import and analyze
NGS-BSA data.

``` r

#load the package
library("PNGseqR")

#Set file and sample names
fl <- "snp.vcf"
MutBulk <- "B73_MU"
WtBulk <- "B73_WT"



#Choose which chromosomes will be included in the analysis (i.e. exclude smaller contigs)
seqname <- c(1:10)


#Import SNP data from file
df <- vcf2table(
        file = fl,
        BulkA = MutBulk,
        BulkB = WtBulk,
        chromlist = seqname
     )

head(df) # The genotype data frame convert from vcf file

#Filter SNPs based on some criteria
data_filt <- BSA_filter(file=df,
                        Bulk.DP=6,
                        Bulk2.Ref=3,
                        Bulk2.Alt=3,
                        Depth.diff=100,
                        min.GQ=99,
                        verbose = TRUE)

head(data_filt) # The remaining SNPs after filter
					 
```

###Four algorithms for BSA analysis, for more detailed instructions please read the supplement materials
``` r		
			 
#BSA analysis

#Run Delta SNP analysis
df_deltaSNP <- DeltaSNP(file = data_filt,
                        Windowsize = 5e6,
                        threshold = 0.995,
                        pval = 0.001)

head(df_deltaSNP$total.result) # The total analysis result 
head(df_deltaSNP$sigthreshold) # The SNPs selected based on the fractile quntile threshold
head(df_deltaSNP$sigpval.result) #The SNPs selected based on the permutation test threshold

#Run Baysian analysis
df_bayes <- Bayesian(file = data_filt,
                     Bulksize = 30,     #The number of samples in the mixing bulks is 30
                     Genomelength = 2338,
                     Windowsize = 5e6,
                     threshold = 0.995,
                     pval = 0.001)

#Run Euclidean Distance analysis
df_ED <- ED(file = data_filt,
            power = 4,
            Windowsize = 5e6,
            threshold = 0.995,
            pval = 0.001)

#Run G-test analysis
df_G <- Gprime(file = data_filt,
               Windowsize = 5e6,
               threshold = 0.995,
               pval = 0.001)
			   
```

#### Figure 1
![](https://github.com/smilejeen/PNGseqR/blob/main/BSA_analysis_plots.jpg)
Figure 1 Candidate region identified by different algorithms in PNGseqR for maize small kernel mutant. Plots produced by the *plot_BSA()* function with a 5 Mb sliding window: (A-D).The scatter plots exported from the BSA results, the used algorithms are Δ(SNP-index), empirical bayes, ED, and G-test. (E-H) the tricube-smoothed values on the genome. Red line threshold shows the region containing SNPs whose tricube-smoothed values are higher than 99.5% of all values, and blue line threshold shows the SNPs whose P-values are lower than 0.001 based on permutation test.
 
``` r			  
#Plot

##deltaSNP
plot_BSA(file = df_deltaSNP$total.result,
         chromlist = seqname,
         type = "DeltaSNP",
         file.name = "DS",
         threshold = 0.995,
              pval = 0.001)
         
##Baysian
plot_BSA(file = df_bayes$total.result,
         chromlist = seqname,
         type = "Bayesian",
         file.name = "Bayesian",
         threshold = 0.995,
         pval = 0.001)

##Euclidean Distance
plot_BSA(file = df_ED$total.result,
         chromlist = seqname,
         type = "ED",
         file.name = "ED",
         threshold = 0.995,
         pval = 0.001)

##Gprime
plot_BSA(file = df_G$total.result,
         chromlist = seqname,
         type = "Gprime",
         file.name = "Gprime",
         threshold = 0.995,
         pval = 0.001)
		 
```

###Plots for DEG analysis and GO analysis,the plot can be found 
#### Figure 2
![](https://github.com/smilejeen/PNGseqR/blob/main/expression_analysis_plots.jpg)
Figure 2 Differential expression analysis and GO analysis for maize small kernel mutant. (A) The volcano plot shows the DEGs for maize small kernel mutant, the analysis was performed by using *DEG_analysis()*. (B) The histogram shows the result of GO analysis for maize small kernel mutant, the analysis was performed by using *GO_analysis()*. Top 20 terms of each mode (BP, CC, MF) are showed in this plot. (C)  A part of directed acyclic plot plotted by using the “MF” mode of *GO_analysis()*. (D)  Bubble plot shows the result of GO analysis, the GO terms are same as (B). 


``` r
#Differential expression analysis

## Set data and parameters
DEG_reads <- DEG_reads.example
rep <- 2
pvalue <- 0.05
log2fold <- 1.5
res <- DEG_analysis(file = DEG_reads,
                    rep.num = rep,
                    sig.level = pvalue,
                    exp.fold = log2fold,
                    file.name = "test",
                    plot = TRUE)

head(res$DEG_total) # The differential expression of all genes
head(res$DEG_sig) # The significant differential expression genes based on the threshold

#GO analysis
## Set data and parameters
sigDEG <- res$DEG_sig
geneid <- sigDEG$id
go_term <- GO_Term.example
mode <- "all"
algorithm <- "classic"
statistic <- "fisher"
num <- 20

GO_result <- GO_analysis(GO_term = go_term,
                         DEG_id = geneid,
                         type = mode,
                         algorithm = algorithm,
                         statistic = statistic,
                         top_num = num,
                         plot = TRUE)

head(GO_result) # Get the top significant gene ontology analysis result

write.table(GO_result,"GO_terms.txt",sep = "\t",quote = F,row.names = F)

#result annotation
## Set data and parameters
gtf <- "zma.v4.example.gtf"
chr <- 8
start <- 79943
end <- 3418075
BSA_position <- data.frame(chr,start,end) # Construct the candidate region of data frame format

sigDEG <- res$DEG_sig

                           
BSA_gene <- result_anno(BSA_data = BSA_position,
                        gtf_data = gtf,
                        DEG_data = sigDEG)

BSA_gene$genes_in_region # The genes fall in the candidate region
BSA_gene$DEGgenes_in_region # The differential expression genes fall in the candidate region

```
<p>&nbsp;</p>
<p>&nbsp;</p>

## Pipeline to get genotype file in VCF format
To help users got genotype and reads count files suitable for our package, we offer command line here:

For DNA sequencing data, the pipeline to get genotype data need install software BWA SAMtools Picard GATK baesd on a Linux system 

``` Linux
First create index for reference genome([Reference.Genome.fa]) by bwa SAMtools and Picard

bwa index [Reference.Genome.fa]

samtools faidx [Reference.Genome.fa]

java -jar picard.jar CreateSequenceDictionary R=[Reference.Genome.fa] O=Reference.Genome.dict
```

Second mapping the sequencing data to the reference genome by bwa and transform result to Binary format

``` Linux
bwa mem -t 16 -M [Reference.Genome.fa] *R1.fq.gz *R2.fq.gz > *.sam
# * R1. fq Gz * R2.fq.gz are the paired-end sequencing data

samtools view -b -S *.sam > *.bam
```


Third The pretreatment for mapping result,including clean sort and adding group information for mapping result, getting unique reads and create index

``` Linux
java -jar picard.jar CleanSam INPUT=*.bam OUTPUT=*.cleaned.bam

java -jar picard.jar FixMateInformation INPUT=*.cleaned.bam OUTPUT=*.cleaned_fixed.bam SO=coordinate

java -jar picard.jar AddOrReplaceReadGroups INPUT=*.cleaned_fixed.bam OUTPUT=*cleaned_fixed_group.bam LB=WH SO=coordinate RGPL=ILLUMINA PU=barcode SM=WH

samtools index *cleaned_fixed_group.bam

java -jar picard.jar MarkDuplicatesWithMateCigar INPUT=*cleaned_fixed_group.bam OUTPUT=*cleaned_fixed_group_DEDUP.bam M=*cleaned_fixed_group_DEDUP.mx AS=true REMOVE_DUPLICATES=true MINIMUM_DISTANCE=500

samtools index *cleaned_fixed_group_DEDUP.bam
```

Then,last perform SNP Calling by GATK for each pool
``` Linux
java -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R [Reference.Genome.fa] -I *cleaned_fixed_group_DEDUP.bam -o *.gvcf 
```

At last, combining the SNP Calling result of different pools to obtain the final file for BSA analysis by GATK
``` Linux
java -jar GenomeAnalysisTK.jar CombineGVCFs -R [Reference.Genome.fa] -V bulk1.gvcf -V bulk2.gvcf -O bsa.gvcf 
# bulk1.gvcf and bulk2.gvcf are SNP Calling result of different pool

java -jar GenomeAnalysisTK.jar GenotypeGVCFs -R [Reference.Genome.fa] -V bsa.gvcf -O bsa.vcf  
```


For RNA sequencing data, the pipeline to get genotype data need install software Hisat2 SAMtools Picard GATK baesd on a Linux system 

First,create index for reference genome by hisat2 
``` Linux
hisat2/extract_exons.py [Reference.Genome.gtf] > exons.txt 
hisat2/extract_splice_sites.py [Reference.Genome.gtf] > splicesites.txt
hisat2/hisat2-build -f [Reference.Genome.fa] --ss splicesites.txt --exon exons.txt [Reference.Genome.index]
# [Reference.Genome.fa] is reference genome，[Reference.Genome.gtf] is the reference genome annotation file，[Reference.Genome.index] is the created reference genome index file
```

Second,mapping RNA sequencing data to reference genome by hisat2
``` Linux
hisat2 --dta -x [Reference.Genome.index] -1 *R1.fq.gz -2 *R2.fq.gz -S *.sam --known-splicesite-infile splicesites.txt
# *R1.fq.gz *R2.fq.gz are the paired-end sequencing data
```

Third,The pretreatment for mapping result,including clean sort and adding group information for mapping result, getting unique reads and create index
``` Linux
samtools view -h -q 50 *.sam | samtools view -bS -> unique.*.bam

samtools sort unique.*.bam unique.*.sort

samtools index unique.*.sort.bam 

java -jar picard.jar AddOrReplaceReadGroups.jar I=unique.*.sort.bam O=unique.*.GR.bam RGID=unique.* RGLB=A RGPL=illumina RGPU=nounit RGSM=unique.*

samtools sort unique.*.GR.bam unique.sort.*.GR

samtools index unique.sort.*.GR.bam
```

At last,perform SNP calling by GATK
``` Linux
java -jar GenomeAnalysisTK.jar -T UnifiedGenotyper -R [Reference.Genome.fa] -I unique.sort.bulk1.GR.bam -I unique.sort.bulk2.GR.bam --heterozygosity 0.1 -stand_call_conf 50.0 -stand_emit_conf 20.0 -glm BOTH --num_threads 2 -ploidy 2 -U ALLOW_N_CIGAR_READS -o bsr_snp.vcf
# unique.sort.bulk1.GR.bam and unique.sort.bulk2.GR.bam is the pretreated mapping result of two different pools
```

count the reads fall in each gene region by HT-seqcount
``` Linux
htseq-count unique.*.sort.bam -f bam [Reference.Genome.gtf] > *.htseq.RC.txt
# unique.*.sort.bam is the mapping result after cleaning and sorting
```
<p>&nbsp;</p>
<p>&nbsp;</p>


















