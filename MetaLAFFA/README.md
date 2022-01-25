# MetaLAFFA

### Overview

MetaLAFFA is a pipeline for annotating shotgun metagenomic data with abundances of functional orthology groups. This process consists of several steps to go from raw FASTQs (with sequencing adapters removed) to functional profiles:

1. Host read filtering (e.g. removing human DNA)
2. Duplicate read filtering
3. Quality trimming and filtering
4. Mapping reads to genes
5. Mapping genes to functional orthology groups (e.g. KOs)
6. Aggregating orthologs into higher-level groupings (e.g. pathways)

More information can be viewed directly on the official tutorial page for MetaLAFFA (see the link to the repository below).

Repository: [https://github.com/engal/MetaLAFFA](https://github.com/engal/MetaLAFFA)

Publication: [Eng, A., Verster, A. J., & Borenstein, E. (2020). MetaLAFFA: a flexible, end-to-end, distributed computing-compatible metagenomic functional annotation pipeline. BMC bioinformatics, 21(1), 1-9.](https://pubmed.ncbi.nlm.nih.gov/33087062/)

## Installation

An issue arose when using the installation instructions on the M3 cluster.
Therefore, please try using the commands below when installing this tool.

```bash
conda clean --all
conda create -n metalaffa python=3.6.10
conda activate metalaffa
conda config --remove channels conda-forge
conda install -c bioconda -c borenstein-lab metalaffa
```
