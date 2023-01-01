# **************tfmice**************


"tfmice" is a R package allowing to build a mice combinatorial transcription factor network enriched from a mice gene list obtained by omics experiment.

### Synopsis
Combinatorial interactions among transcription factors are critical to directing tissue-specific gene expression [(Ravasi T. 2010)](https://pubmed.ncbi.nlm.nih.gov/20211142/). Tissue specificty is enabled by spatial and temporal patterns of gene expression which in turn are driven by transcriptional regulatory networks [(Naef and Huelsken 2005; Zhang 2004)](https://pubmed.ncbi.nlm.nih.gov/16030348/). Transcription factors do not act independently but form complexes with other transcription factors, chromatin modifiers and cofactor proteins, which bind together and assemble upon the regulatory regions of DNA to affect transcription [(Fedorova and Zink 2008)](https://pubmed.ncbi.nlm.nih.gov/18718493/). It has been estimated that approximatively 75 percent of metazoan transcription factors could heterodimerize with other factors [(Walhout 2006)](https://pubmed.ncbi.nlm.nih.gov/17053092/). It is important to understand which combination of transcription factors act together and how these combinations lead to modes of regulation that are not evident when each factor is considerate separately [(Ravasi T. 2010)](https://pubmed.ncbi.nlm.nih.gov/20211142/).  

### Package installation
```r
library(devtools)
install_github("cdesterke/tfmice")
```

### Description

Starting from an imput mice gene list obtained by omics experiments this tfmice package allows to perform a combinatorial transcription factor enrichment.

- In a first step compare your input genelist to mice TF-TF interaction database [(Ravasi Timothy 2010)](https://pubmed.ncbi.nlm.nih.gov/20211142/). This initial step of enrichment performs with "tfcrcalc" function allows to output a result table with significance of enrichments. 

- In second step you could output your enrichment results in barplot for "enrichment scores with "plotes" function and for negative log10 p-values with "plotnlp" function. 

- In last step you can performed a combinatorial transcription factor network with Louvain communitie detection: fonction "tfcrnet". Significant families during the enrichment were represented with red edges and not significant families with gray edges. 

## code for processing analysis

```r
##define a mice gene list
library(tfmice)
custom<-c("Hoxa10","Hoxa11","Hoxa13","Hoxa13","Hoxa2","Hoxa5","Hoxa9","Hoxb13","Hoxb4","Hoxb8","Hoxc13","Hoxd12","Hoxd13","Ikbkb")
## calculates the family enrichments
res<-tfcrcalc(custom)
head(res,n=15)
```
#![res](https://github.com/cdesterke/tfmice/blob/main/res.png)


```r
res<-head(res,n=15)
plotes(res)

```
#![es](https://github.com/cdesterke/tfmice/blob/main/es.png)


### barplot of negative log10 p-values

```r
res<-head(res,n=15)
plotnlp(res)

```
![nlp](https://github.com/cdesterke/tfmice/blob/main/nlp.png)


## network of enriched combinatorial transcription factors


### custom the network with some parameters

- cex parameter: change size of the vertex (nodes) label

- distance parameter: change the distance between the vertex and its label

- layout parameter: change the design of the network and have several options such as: (layout_as_star, layout_components, layout_in_circle, layout_nicely,layout_on_grid,
layout_on_sphere, layout_randomly, layout_with_dh, layout_with_drl, layout_with_fr, layout_with_gem,
layout_with_graphopt, layout_with_kk, layout_with_lgl, layout_with_mds)

output with nicely layout

```r
tfcrnet(custom,res,layout=layout_nicely,cex=1,distance=1.5)

```
#![nicely](https://github.com/cdesterke/tfmice/blob/main/netnicely.png)

output with star layout

```r
tfcrnet(custom,res,layout=layout_as_star,cex=1,distance=1.5)

```
#![star](https://github.com/cdesterke/tfmice/blob/main/netstar.png)




# REFERENCES

> Timothy Ravasi, Harukazu Suzuki, Carlo Vittorio Cannistraci, Shintaro Katayama, Vladimir B Bajic, Kai Tan, Altuna Akalin, Sebastian Schmeier, Mutsumi Kanamori-Katayama, Nicolas Bertin, Piero Carninci, Carsten O Daub, Alistair R R Forrest, Julian Gough, Sean Grimmond, Jung-Hoon Han, Takehiro Hashimoto, Winston Hide, Oliver Hofmann, Atanas Kamburov, Mandeep Kaur, Hideya Kawaji, Atsutaka Kubosaki, Timo Lassmann, Erik van Nimwegen, Cameron Ross MacPherson, Chihiro Ogawa, Aleksandar Radovanovic, Ariel Schwartz, Rohan D Teasdale, Jesper Tegnér, Boris Lenhard, Sarah A Teichmann, Takahiro Arakawa, Noriko Ninomiya, Kayoko Murakami, Michihira Tagami, Shiro Fukuda, Kengo Imamura, Chikatoshi Kai, Ryoko Ishihara, Yayoi Kitazume, Jun Kawai, David A Hume, Trey Ideker, Yoshihide Hayashizaki: Cell
. 2010 Mar 5;140(5):744-52. doi: 10.1016/j.cell.2010.01.044., An atlas of combinatorial transcriptional regulation in mouse and man (PMID: 20211142)

> Felix Naef, Joerg Huelsken: Nucleic Acids Res. 2005 Jul 19;33(13):e111. doi: 10.1093/nar/gni104., Cell-type-specific transcriptomics in chimeric models using transcriptome-based masks (PMID: 16030348)

> Elena Fedorova 1, Daniele Zink: Biochim Biophys Acta. 2008 Nov;1783(11):2174-84. doi: 10.1016/j.bbamcr.2008.07.018, Nuclear architecture and gene regulation (PMID: 18718493)

> Albertha J M Walhout: Genome Res. 2006 Dec;16(12):1445-54. doi: 10.1101/gr.5321506. Epub 2006 Oct 19., Unraveling transcription regulatory networks by protein-DNA and protein-protein interaction mapping (PMID: 17053092)


## complete use case

```r
library(devtools)
install_github("cdesterke/tfmice")


library(tfmice)

#save the imported gene list in a string vector
vector<-set$gene

#compute the score of enrichment
res<-tfcrcalc(vector)

#verified if the enriched tf are present in the initial genelist
library(dplyr)
df<-res
df$gene<-row.names(df)
set%>%inner_join(df,by="gene")->union
row.names(union)<-union$gene

#plot top 15 TF 
plotes(union)
plotnlp(union)

#remove TF with no significance
union%>%filter(pvalues!=1)->union

#adjust the table for the net
union$gene<-NULL

#perform the net
tfcrnet(vector,union,layout=layout_nicely,cex=1,distance=1.5)

```
