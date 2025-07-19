This repository contains the code and data to reproduce the bulk RNA-seq and promoter analyses presented in the
paper ["Single-cell analysis of signalling and transcriptional responses to type I interferons"](https://www.biorxiv.org/content/10.1101/2023.07.03.547491v1.full/).

-----

## Repository Overview

Here is a high-level overview of the repository's structure:

```
ifn-resource/
├── assemblies/GRCh38 # GRCh38 reference genome and processed annotations (populated during setup)
├── pixi.lock         # Pixi lock file for reproducible environments
├── pixi.toml         # Pixi project configuration, dependencies, and tasks
├── setup/            # Scripts for initial data download and setup
├── utils/            # Utility scripts and modules (e.g., for FASTA, BED, motif operations)
└── stories/          # Main analysis workflows
```

All analysis is organized into biological or technical "stories," such as Nextflow data preprocessing or differential
expression analysis. When necessary, stories are divided into sub-stories. Each story has a corresponding `ld` (local
data) directory that contains results, intermediate files, and other resources. The analysis steps, which involve
running various bash or Python scripts, are encapsulated as **Pixi** tasks for reproducibility and ease of use.

-----

## Installation & Setup

1. **Install Pixi:** The computational environment and all dependencies are managed by **Pixi**. Install it by following
   the instructions at [pixi.sh](https://pixi.sh/).
2. **Clone the Repository:**
   ```bash
   git clone https://github.com/alnfedorov/ifn-resource
   cd ifn-resource
   ```
3. **Set up the repository:**
   ```bash
   pixi run setup/download-gencode # Download the GRCh38 assembly and GENCODE annotation
   pixi run setup/parse-gencode    # Parse and index the annotation
   pixi run setup/terminus         # Compile Terminus
   pixi run setup/download-fastq   # Download the FASTQ files for bulk RNA-seq
   ```
   Other resources, such as ENCODE cCREs and JASPAR motifs, are commited directly to the repository.

**Note:** This repository has been tested on Linux (x86_64) and may not be compatible with other operating systems or
architectures.

-----

## Workflow

To reproduce the full analysis, execute the following **Pixi** tasks in order:

```shell
# Run the Nextflow preprocessing pipeline
pixi run nextflow/full          # Run using the full GENCODE annotation (except level 3)
pixi run setup/filter-gencode   # Filter the GENCODE annotation based on alignment results
pixi run nextflow/filtered      # Rerun using the filtered GENCODE annotation

# Run the downstream analysis
pixi run stories/terminus           # Collapse technically indistinguishable transcripts using Terminus
pixi run stories/DE                 # Perform differential expression analysis using DESeq2
pixi run stories/cCRE               # Derive promoter regions based on ENCODE cCREs
pixi run stories/STREME             # Run de novo motif discovery using STREME
pixi run stories/JASPAR/scoring     # Score promoter sequences using JASPAR motifs
pixi run stories/JASPAR/association # Identify associations between motifs and expression changes
```

-----

## Dependencies

All Python and system dependencies are defined in `pixi.toml` and managed by **Pixi**. They will be installed
automatically into an isolated environment when you run any **Pixi** command (e.g., `pixi run <task>`) or run
`pixi install` directly.
