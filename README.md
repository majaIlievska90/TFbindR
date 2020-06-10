# TFbindR

Counting the occurence of transcription factor binding motifs in a set of provided sequences upstream from the translation start site. This is used for performing motif enrichment analysis. First, the occurence matrix is constructed, by scanning FIMO files or a list of regular expressions representing the transcription factor binding motifs, and the set of sequences in FASTA format. The occurence is identified for both, the sense and antisense strands. The enrichment analysis performs a Fisher's exact test for identification of enriched motifs in a set of genes of interest.

## Installation

```
library(devtools)
install_github("https://github.com/majaIlievska90/TFbindR”)
library(TFbindR)
```

From here you can read the documentation for each function with  ``` ?TFbindR::function_name  ```, 
e.g  ```  ?TFbindR::makePromoterDBfromFimo() ```


Before using the functions unzip the package data, you will find the path for the TFbindR from  ``` system.file(package="TFbindR") ```, open the directory, go to the extdata folder and unpack the gz files. 


## Using the lbirary

Before computing the motif enrichment, you need to create an **R object** that contains information on the occurrence of the motifs in your set of sequences (usually promoter sequences). This is done based on FIMO motifs or Regular Expression.

**For the FIMO version**

In this version, you use the function ```makePromoterDBfromFimo()```, which implements a count of the Occurrence of Sequence Motifs from a file produced with the FIMO tool from the [MEME suit](http://meme-suite.org). 
If working with plant genome sequences, you can use our fimo file which is based on a collection of *Arabidopsis thaliana* motifs from 2 databases, [Plant Cistrome database](http://neomorph.salk.edu/dap_web/pages/index.php) and [JASPAR](http://jaspar.genereg.net). In this case you can skip step 1. and start directly with 2. If you want to use your custom collection, go through step 1.
 
1. The MEME suit motif-based sequence analysis tools have been developed by William Noble's lab at the University of Washington. Among the tools, FIMO scans a set of sequences for individual matches to a set of motifs that are provided as input.  You will need the [MEME formatted motifs](http://meme-suite.org/doc/meme-format.html), and the FASTA file with the sequences of interest.  To complete this step, visit the [MEME suit FIMO web version](http://meme-suite.org/tools/fimo) and upload your files. Run the tool with the default settings. The output will be saved in a **fimo_out folder** with a **fimo.tsv** file which is the one you will need further on.


2. Run the ```makePromoterDBfromFimo``` function which requires the following input arguments: **fimo** - the output of the fimo tool, **promoter.fasta** - fasta file of the promoter sequences and  **promoterDB** - the path to the directory where you want the output to be saved.  With our FIMO file, 
``` fimo.dir =system.file("extdata", "cistrome.fimo.txt", package="TFbindR", mustWork = TRUE) ``` 
If using your own FIMO file, replace the previous line with a simple line providing the path to your file.
```fimo.dir="path_to_your_FIMO_file/your_FIMO_file_name" ```
Next, provide the path to the FASTA file. If using the *Arabidopsis* promoter sequences, 1000 bp upstream of the translation start site, run this line
```promoter.dir =system.file("extdata","TAIR10_upstream_1000_translation_start_20101028.fa",package="TFbindR")```
If using your own promoter sequences, provide the path. 
```promoter.dir="path_to_your_FASTA_file/your_FASTA_file_name"```
Give the path for the result. Use the suffix to your file name .fimo.db.rda
```output.dir="path_to_result_directory/result_file_name.fimo.db.rda" ```.
**Finally, run the function**
```makePromoterDBfromFimo(fimo = fimo.dir, promoter.fasta=promoter.dir, promoterDB=file.path(output.dir))```


**For the Regular Expression (RE) version**

Here again, you can use the *Arabidopsis* collection of motifs defined with regular expression patters. The collection is based on [AtcisDB - Arabidopsis cis-regulatory element database](https://agris-knowledgebase.org) and secondary cell wall specific binding elements from literature. 

If you use your own motifs, replace the inputs accordingly. 

```
motifs.dir = system.file("extdata", "REcomplete.txt", package="TFbindR", mustWork = TRUE)
Motif.df=read.delim(file = motifs.dir ,header=TRUE, sep="\t", as.is=T)
promoter.dir = system.file("extdata","TAIR10_upstream_1000_translation_start_20101028.fa",package="TFbindR")

makePromoterDBfromRE(Motifs = Motif.df,promoter.fasta = promoter.dir, promoterDB=file.path(output.dir,"TAIR10_500bp_upstream.RE.db.rda"))
```

The result will be saved in the ```TAIR10_500bp_upstream.RE.db.rda``` object containing the frequency of the moitfs in each promoter sequence. 

**Enrichment**

From here, you can run the enrichment function. Example bellow for the RE version.  Here, you will use a text file with the gene IDs separated by new line, without header. 

```
genes=read.table("path_to_your_file/file_name")
enrichment=computeEnrichments(as.character(genes), "TAIR10_500bp_upstream.fimo.db.rda")
```


