# SciBet
A portable and fast single cell type identifier

## Contents

- [Overview](#overview)
- [Installation Guide](#installation-guide)
- [Tutorial](#tutorial)
- [Reproduction instructions](#Reproduction-instructions)
- [License](./LICENSE)
- [Citation](#citation)
- [Contact](#Contact)

## Overview
Fast, robust and technology-independent computational methods are needed for supervised cell type annotation of single-cell RNA sequencing data. We present SciBet, a single-cell classifier that accurately predicts cell identity for newly sequenced cells or cell clusters. We enable [web client deployment of SciBet](http://scibet.cancer-pku.cn/) for rapid local computation without uploading local data to the server. This user-friendly and cross-platform tool can be widely useful for single cell type identification (bioRxiv preprint link [here](https://www.biorxiv.org/content/10.1101/645358v2)).

## Installation Guide
**Installing dependency package**  
Before installing SciBet, the “tidyverse” package should be installed first:
```
install.packages("tidyverse")
```
**Installing ROGUE**  
To install ROGUE, run:
```
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("PaulingLiu/scibet")
```

## Tutorial
For more details and basic usage see following tutorials:
1.	[Guided Tutorial](http://scibet.cancer-pku.cn/document.html)

## Reproduction instructions
The scripts for producing all the quantitative results in our manuscript can be found in [scripts](./scripts).

## Citation
If you use SciBet in your research, please considering citing:
- [Li et al., bioRxiv 2019](https://www.biorxiv.org/content/10.1101/645358v2)

## Contact
Please contact us:  
Baolin Liu: pauling.liu@pku.edu.cn
Chenwei Li: lichenwei@pku.edu.cn
Zemin Zhang: zemin@pku.edu.cn

## Copyright
©2019 Chenwei Li, Baolin Liu. [Zhang Lab](http://cancer-pku.cn/). All rights reserved.
