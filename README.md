# CrossTx

## 1. Introduction

**CrossTx** is an R package designed to address a central challenge in translational cancer research: quantitatively evaluating the **transcriptional similarity between mouse models and human cancer samples**. It provides a streamlined pipeline for cross-species integration of RNA-seq data, enabling researchers to assess the fidelity of preclinical models in mimicking human tumor biology.

By integrating DEG analysis, PCA, TROM scoring, immune profiling, and GO enrichment in a unified workflow, CrossTx helps researchers systematically identify **which mouse models best resemble human tumor subtypes** at the transcriptomic level.

### CrossTx package workflow

![workflow](./figure/workflow.png)

---

## 2. Installation

You can install the package directly from GitHub:

```r
# Install devtools if necessary
if (!require("devtools")) install.packages("devtools")

# Install the package
devtools::install_github("wangdian-PKU/CrossTx")
```

---

## 3. User Manual

Detailed guidance of how to use CrossTx could be found in the https://wangdian-pku.github.io/CrossTxManual/.

---

## 4. Contact us

We are always open to your questions and encourage discussions here.

You can report potential bugs to the [issues](https://github.com/wangdian-PKU/CrossTx/issues).

Any question can directly contact the author by email: wangdian020803@gmail.com
