# InvenireSRNA
Small RNAs (sRNAs) constitute an important class of post-transcriptional regulators that control critical cellular processes in bacteria. While recent research has led to a dramatic increase in the discovery of bacterial sRNAs, it is generally believed that the currently identified sRNAs constitute a limited subset of the bacterial sRNA repertoire. In several cases, sRNAs belonging to a specific class are already known and the challenge is to identify additional sRNAs belonging to the same class. **InvenireSRNA** is an R package for learning a classification model for a given class of sRNA thus allowing for the discovery of additional sRNAs beloning to the same class. **InvenireSRNA** also provides a pretrained model for predicting RsmA/CsrA regulating sRNAs.

## Installation
Before installing this package, make sure you have the latest version of *Rstudio*, *R* and the *devtools* package. You also need to have C++11 available on your machine for the algorithm to run properly. The python package *Biopython* is also required for running the algorithm. Finally, you will need to install the [ViennaRNA](http://www.tbi.univie.ac.at/RNA/) package. You can install this R pacakge using the following:
```{R}
library(devtools)
install_github("carltonyfakhry/InvenireSRNA")
```

## Usage
For an introduction to *InvenireSRNA*, please see the Vignette for this package using the following:

```{R}
browseVignettes("InvenireSRNA")
```

## Webserver
An interactive webserver providing some of the functionality as this R package is available at [InvenireSRNA web server](http://markov.math.umb.edu/inveniresrna/).

## Citation
[1] Carl Tony Fakhry, Prajna Kulkarni, Ping Chen, Rahul Kulkarni and Kourosh Zarringhalam (2017). "Prediction of bacterial small RNAs in the RsmA (CsrA) and ToxT pathways: a machine learning approach." BMC Genomics, 18.

[2] Carl T. Fakhry, Kourosh Zarringhalam, and Rahul V. Kulkarni. "Bioinformatic Approach for Prediction of CsrA/RsmA-Regulating Small
RNAs in Bacteria." Methods in Molecular Biology, pp. 47-56. Humana
Press, New York, NY, 2018.
