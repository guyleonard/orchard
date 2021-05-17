# Orchard Pipeline Introduction
A Rapid Phylogenetic Tree Building Pipeline.

'Seed' Sequences --> Search DB --> Alignment --> Masking --> Phylogenetic Tree Reconstruction (+-> SVG/PDF Tree Images --> Taxon Annotation)

This package is intended to create a set of "pilot" trees, from a given a set of initial 'seed' amino acid sequences. It is separated into main two components:
 1. [orchard](#orchard)
 2. [orchard_accessories](#orchard_accessories)

The first program is responsible for all of the essential tree building stages, requiring as input:
 1. a set of amino acid sequences contained in a FASTA file,
 2. a list of taxon IDs from an [orchardDB](https://github.com/guyleonard/orchardDB) database,
 3. a parameters file to set the options for each stage/sub-program.

The second program is responsible for all the annotation of the trees. It will convert the computer-readable accessions from the [orchardDB](https://github.com/guyleonard/orchardDB) back to their original accessions and taxon names. Currently, it can also draw PDF and SVG trees providing Dendroscope and Inkscape have been installed. Eventually it will also annotate PFAM domains onto the PDFs.

This package relies on an [orchardDB](https://github.com/guyleonard/orchardDB) database. There are several to choose from or you can build your own.

[Insert Schematic]

Continue below to:
 * [Installation](#Installation)
 * [Program Execution](#Program-Execution)
 * [Example Parameters File](#Example_Parameters_File)
 * [Citation](#Citation)
 * [History](#History)

# Installation
The set of dependencies for the 'orchard' program can mostly be installed via conda. Dependencies for 'orchard_accessories', especially to produce PDF/SVG trees, may require some extra steps for dependencies outside of conda.

## Easy Steps
This is the minimal setup that assumes you want to run the 'orchard' script with default programs, and the basics of the 'orchard_accessories' script (e.g. no PDF/SVG output).

### orchardDB
Please set up an orchardDB from [here](https://github.com/guyleonard/orchardDB). You may like to use our *cider* DB to start with, or build your own.

### Conda
```bash
conda create -n orchard
conda activate orchard
conda install -c bioconda blast=2.11.0 mafft trimal fasttree perl-bioperl perl-io-prompt
```

## Advanced Steps
See below...

# Program Execution
As mentioned previously the 'orchard' package is separated in to two components:
 1. [orchard](#orchard)
 2. [orchard_accessories](#orchard_accessories)

## orchard
This program is responsible for all of the essential tree building stages, requiring as input:
 1. a set of amino acid sequences contained in a FASTA file,
 2. a list of taxon IDs from an [orchardDB](https://github.com/guyleonard/orchardDB) database,
 3. a parameters file to set the options for each stage/sub-program.

```bash
Usage: orchard -i seqs.fasta -t taxa.txt -p params.yaml -s -a -m -x

Required Parameters
	-i <sequences.fasta>	Input Sequences in FASTA Format
	-t <taxa_list.txt>	List of orchardDB Taxa IDs
	-p <parameters.yaml>	Parameters File in YAML Format
Optional Parameters
	-s	Run Searches e.g. BLASTp
	-a	Run Alignments e.g. MAFFT
	-m	Run Masking e.g trimAl
	-x	Run Trees e.g FastTree2
```

## orchard_accessories
This program is responsible for all the annotation of the trees. It will convert the computer-readable accessions from the [orchardDB](https://github.com/guyleonard/orchardDB) back to their original accessions and taxon names. Currently, it can also draw PDF and SVG trees providing Dendroscope and Inkscape have been installed. Eventually it will also annotate PFAM domains onto the PDFs.

```bash
Usage: orchard_accessories -p params.yaml -n

Required Parameters
	-p <parameters.yaml>
Optional Parameters
Renaming
	-s	Rename Taxa in Sequence Hits Files
	-a	Rename Taxa in Alignment Files
	-m	Rename Taxa in Masked Files
	-e	Rename Taxa in Excluded Files
	-n	Rename Taxa in Newick Trees
Tree Drawing
	-renamed	Convert Renamed Trees use with below
	-eps_tree	Draw a Phylogram in EPS Format (Basic)
	-svg_tree	Draw a Phylogram with Dendroscope in SVG Format
	-pdf_tree	Draw a Phylogram with Dendroscope in PDF Format
```

# Example Parameters File
This file contains all the default and changeable values for each sub-program included in 'orchard' and 'orchard_accessories'. You should make a new one for each analysis you conduct.

It is in the 'YAML' format, allowing for good readability for huamns and machines aliked. Each sub-program/step of 'orchard' has a section, and each setion has sub-sections which are the options. Anything with a '#' infront of it is a comment and will be ignored by the programs.

An example file is given below, for the most part you can leave the defaults as we have set them. Pay special attention to 'threads', 'e-value' and the values for 'directories' & 'databases' as these will be specific to your setup and analysis.

```yaml
user:
  results: testing    # folder name for results
  threads: 8          # number of cores to use

search:
  program: blast+     # blast+, blat
  subprogram: blastp  # blast: blastp, blastx. blat: dna, prot, dnax
  evalue: 1e-10       # 1e-10, 1e-05, etc
  tophits: 5          # integer
  maxlength: 3000     # ignore AA sequences larger than this value

special:
  taxa:               # comma separated list of genome IDs from orchardDB
  tophits:            # integer

alignment:
  program: mafft      # mafft, muscle
  options: --auto --quiet --reorder  # mafft: '--auto --quiet --reorder' or muscle: '-maxiters 2 -quiet -group'

masking:
  program: trimal     # trimal, divvier
  cutoff1: 50         # integer
  cutoff2: 30         # integer

trees:
  program: fasttreemp # fasttree, fasttreemp, iqtree
  options: -bionj -slow -lg -quiet  # fasttree: '-bionj -slow -lg -quiet' or iqtree: '-fast -alrt 1000 -quiet -mset WAG,LG,JTT -merit BIC'
  mintaxa: 3          # integer

directories:
  orchardDB: /home/cs02gl/Dropbox/git/orchardDB/testing/cider  # path to orchardDB SQL folder - .sql file must be the same name as the folder

database:
  username: test      # orchardDB username
  password: test      # orchardDB password
 ```

## Advanced Installation
This has a more indepth set of information pertaining to the installation of certain programs and other non-default options.

### Perl Modules
* Bio::DB::Fasta
* Bio::Perl
* Bio::SearchIO
* Bio::SeqIO
* Cwd
* DBI
* File::Basename
* File::Path
* File::Slurp
* Getopt::Long
* IO::Prompt
* IO::Tee
* YAML::XS

 e.g. you can install these by typing:
```bash
$ sudo cpanm Bio::DB::Fasta Bio::Perl Bio::SearchIO Bio::SeqIO Cwd DBI File::Basename File::Path File::Slurp Getopt::Long IO::Prompt IO::Tee YAML::XS
```

### Conda
```bash
conda create -n orchard
conda activate orchard
conda install -c bioconda blast=2.11.0 blat diamond mafft muscle trimal fasttree iqtree perl-bioperl perl-io-prompt
conda install -c conda-forge inkscape
```

### Manual Software Installation
One from each section, or all of the options below must be installed; the defaults/prefered are Blast+, MAFFT, trimAL, FastTree2.

#### Search Programs
 1. [Blast+](http://blast.ncbi.nlm.nih.gov/Blast.cgi/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
 2. [Blat](https://genome.ucsc.edu/FAQ/FAQblat.html)
 3. [Diamond](https://github.com/bbuchfink/diamond)
   1. Diamond DBs must be created in your orchardDB directory using the same version of Diamond that you use for searching.

#### Alignment Programs
 1. [MAFFT](http://mafft.cbrc.jp/alignment/software/)
 2. [MUSCLE](http://www.drive5.com/muscle/)

#### Masking Programs
 1. [trimAl v1.4](http://trimal.cgenomics.org/)
 2. [Divvier](https://github.com/simonwhelan/Divvier)

#### Tree Building Programs
 1. [FastTree2](http://meta.microbesonline.org/fasttree/)
 2. [IQTree](https://github.com/Cibiv/IQ-TREE)

#### SVG Tree Creation
 1. [Dendroscope 3](http://ab.inf.uni-tuebingen.de/software/dendroscope/)
 2. XVFB - Needed for running on a headless server.

```bash
sudo apt-get install xvfb
```

#### PDF Creation (SVG to PDF)
 1. [Inkscape](https://www.inkscape.org/en/)

# Citation
Please cite the original manuscript, and the updated repository DOI:

 * Thomas A. Richards, Darren M. Soanes, Peter G. Foster, Guy Leonard, Christopher R. Thornton and Nicholas J. Talbot. (2009). Phylogenomic Analysis Demonstrates a Pattern of Rare and Ancient Horizontal Gene Transfer between Plants and Fungi. The Plant Cell. 21(7). www.plantcell.org/cgi/doi/10.1105/tpc.109.065805
 * DOI: 

# History
This pipeline was first reported in Richards et. al. (2009) as a "an automated gene-by-gene phylogeny pipeline to generate a PhyML tree" and created by Darren Soannes. Since then it has gone through several iterations and changes, but remains similar in its goal. Between 2009 and 2013 there were quite a few changes in the processing of the script as we had moved to another institution and took development with us. None of these were particularly well documented, other than in internal code comments, although at all times the code was available on request or via github. What you see now is a total rewrite, and should be useful to other researchers without too many issues, I will accept pull requests and colaborations.