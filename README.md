# BACDAS.
# Amplicon Sequence Analysis Tool Scripts

This repository contains a collection of Bash scripts for analyzing amplicon sequences using various bioinformatics tools. Each script is designed to perform a specific analysis or task, and they can be used individually or combined as part of a larger workflow.

## Table of Contents

- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
- [Scripts](#scripts)
- [Contributing](#contributing)
- [License](#license)

## Prerequisites

Before using these scripts, ensure that you have the following software and dependencies installed:

- [Bash](https://www.gnu.org/software/bash/) (version 4.0 or higher)
- [Usearch](https://www.drive5.com/usearch/download.html)
- [Mothur](https://mothur.org/wiki/download_mothur/)
- [Deblur](https://github.com/biocore/deblur)
- [DADA2](https://benjjneb.github.io/dada2/dada-installation.html)
- [MED](https://merenlab.org/software/oligotyping/)
- [prinseq]
- [cutPrimers](https://github.com/aakechin/cutPrimers)
- [blastn]

Each script may have additional prerequisites specific to the tool it utilizes. Refer to the documentation of the individual tools for more information.

## Installation

1. Clone this repository to your local machine or server:

   ```bash
   git clone https://github.com/MOFares-Bioinf/BACDAS.git
   ```

2. Navigate to the cloned repository:

   ```bash
   cd BACDAS
   ```

3. Install any necessary dependencies for the tools you plan to use (refer to the tool's documentation for instructions).

## Usage

Each script in this repository serves a specific purpose and can be executed independently. To run a script, you have to edit the file and change file paths to path of your working directory
save the changes and use the following command:

```bash
bash script_name.sh 
```

Replace `script_name.sh` with the actual script filename.

## Scripts

The following scripts are available in this repository:

- **preprocessing.sh:** Run the preprocessing pipline for quality filtering trimming and preparation of input reads.
- **run_usearch.sh:** Excutes Uparse and Unoise3 pipelinge on the input data.
- **run_mothur.sh:** Excutes Mothur pipeline using Vsearch-DGC, Average neighbrhood and opticlust clustering methods.
- **run_MED.sh:** Excutes MED oligotyping pipeline
- **run_deblur.sh:** Excutes MED oligotyping pipeline

You can find detailed information about each script in their respective files.

## Contributing

Contributions to this repository are welcome. If you have any improvements, bug fixes, or new scripts to contribute, please follow these steps:

1. Fork the repository.
2. Create a new branch for your feature or bug fix.
3. Make the necessary changes and commit your code.
4. Push your branch to your forked repository.
5. Open a pull request with a detailed description of your changes.

## License

This project is licensed under the [MIT License](LICENSE). Please review the license file for more information.
