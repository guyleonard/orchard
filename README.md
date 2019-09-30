# Orchard Pipeline
A Quick Phylogenetic Tree Building Pipeline...

Seed Sequences --> DB Search --> Alignment --> Masking --> Phylogenetic Tree Reconstruction --> SVG/PDF Tree Images --> Taxon Annotation

This set of scripts is intended to create a group of "pilot" trees, given your initial gene(s) of interest. The resulting Phylogenies can be easily viewed, and then manually improved as needed. The aim is speed, in order to calculate 1000s of trees for a dataset in reasonable time.

# Dependencies
## orchardDB
Please set up an orchardDB from [here](https://github.com/guyleonard/orchardDB). You may like to use the *cider* DB scritps to start, or build your own.

## Perl Modules
All of: Bio::Perl Cwd DBI File::Basename File::Path File::Slurp Getopt::Long IO::Prompt IO::Tee YAML::XS

e.g. you may install these by typing:
```bash
$ sudo cpanm Bio::Perl Cwd DBI File::Basename File::Path File::Slurp Getopt::Long IO::Prompt IO::Tee YAML::XS
```

## Software
### Search Programs
One each, or all of the options below.

1. [Blast+](http://blast.ncbi.nlm.nih.gov/Blast.cgi/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
2. [Blat](https://genome.ucsc.edu/FAQ/FAQblat.html)
3. [Diamond](https://github.com/bbuchfink/diamond)
   1. Diamond DBs must be created in your orchardDB directory using the same version of Diamond for searching.

### Alignment Programs
1. [MAFFT](http://mafft.cbrc.jp/alignment/software/)
2. [MUSCLE](http://www.drive5.com/muscle/)

### Masking Programs
1. [trimAl v1.4](http://trimal.cgenomics.org/)
2. [Divvier](https://github.com/simonwhelan/Divvier)

### Tree Building Programs
1. [FastTree2](http://meta.microbesonline.org/fasttree/)
2. [IQTree](https://github.com/Cibiv/IQ-TREE)

### Renaming Taxa
1. [NCBI Taxonomy](ftp://ftp.ncbi.nih.gov/pub/taxonomy) - retrieve 'taxdump.tar.gz' from ftp://ftp.ncbi.nih.gov/pub/taxonomy

### SVG Tree Creation
1. [Dendroscope 3](http://ab.inf.uni-tuebingen.de/software/dendroscope/)
   1. XVFB - Needed for running on a headless server.

### PDF Creation (SVG to PDF)
1. [Inkscape](https://www.inkscape.org/en/)

e.g. You may install many of these with Conda (divvier, dendroscope and NCBI Taxdump will need to be manually installed)
```bash
conda install -c bioconda blast blat diamond mafft muscle trimal fasttree iqtree
conda install -c conda-forge inkscape
sudo apt-get install xvfb
```

# Program Execution
You will need four things to begin.

1) An 'orchardDB' folder and sqlite database.
2) Your initial search 'seed' sequence(s) in FASTA format.
3) A taxa list, based on the 'genome ID' from your orchardDB. One per line in a text file.
4) A parameters file, as seen below, with your options in YAML format.

You may run any Optional Parameter providing the previous output exists.

## orchard
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

### Example Output Directory Structure
```
.
└── output
    ├── seqs       # all hits to search seed
    ├── alignments # aligned hits
    ├── masking    # masked/trimmed hits
    ├── trees      # newick trees
    ├── excluded   # seeds that were excluded due to too few hits or poor alignment
    └── reports    # blast/search program reports
```

## orchard_accessories
```bash
Usage: orchard_accessories -p params.yaml -n
Required Parameters
  -p <parameters.yaml>
Optional Parameters
Renaming
  -s  Rename Taxa in Sequence Hits Files
  -a  Rename Taxa in Alignment Files
  -m  Rename Taxa in Masked Files
  -e  Rename Taxa in Excluded Files
  -n  Rename Taxa in Newick Trees
Tree Drawing
  -eps_tree Draw a Phylogram in EPS Format (Basic)
  -svg_tree Draw a Phylogram with Dendroscope in SVG Format
  -pdf_tree Draw a Phylogram with Dendroscope in PDF Format
Modify Tree Drawing
  -r  After -n use this with Tree Drawing options to use renamed tees
Annotation
  -c  Colourise taxa in SVG tree according to taxonomy colours in parameters file
SVG Conversion
  -p  Convert colourised SVG to PDF
Cite: https://github.com/guyleonard/orchard and doi: 10.1105/tpc.109.065805

```

### Example Output Directory Structure
```
.
└── testing
    ├── seqs
    │   └── renamed
    ├── alignments
    │   └── renamed
    ├── masking
    │   └── renamed
    ├── excluded
    │   └── renamed
    ├── reports
    └── trees
        ├── eps
        ├── pdf
        ├── svg
        └── renamed
            ├── eps 
            ├── pdf
            └── svg
                └── colourised
```

## Example Parameters File
Why a YAML file and not command line options?

1) Because it is human readable and,
2) then you have a record of all settings for each analysis.

```yaml
user:
  results: testing    # folder name for results
  threads: 8          # number of cores to use

search:
  program: blast+     # blast+, blat, diamond
  subprogram: blastp  # blast: blastp, blastx; blat: prot, dnax; diamond: blastp, blastx
  evalue: 1e-10       # 1e-10, 1e-05, etc
  tophits: 5          # integer
  maxlength: 3000     # ignore AA sequences larger than this value

special:
  taxa:               # comma separated list of genome IDs from orchardDB
  tophits:            # integer

alignment:
  program: mafft                     # mafft, muscle
  options: --auto --quiet --reorder  # mafft: '--auto --quiet --reorder' or muscle: '-maxiters 2 -quiet -group'

masking:
  program: trimal     # trimal, divvier
  cutoff1: 50         # higher integer
  cutoff2: 30         # lower integer

trees:
  program: fasttreemp               # fasttree, fasttreemp, iqtree
  options: -bionj -slow -lg -quiet  # fasttree: '-bionj -slow -lg -quiet' or iqtree: '-fast -alrt 1000 -quiet -mset WAG,LG,JTT -merit BIC'
  mintaxa: 3                        # integer

directories:
  orchardDB: /path/to/cider_db  # path to orchardDB SQL folder - .sql file must be the same name as the folder

database:
  username: test      # orchardDB username
  password: test      # orchardDB password

annotation:
  taxonomy_colours: Alveolata;#cc6677,Amoebozoa;#44aa99,Archaea;#696969,Bacteria;#343434,Capsaspora_owczarzaki;#005eb0,Choanoflagellatea;#88ccee,Cryptophyta;#ec660d,Discoba;#ff0000,Fonticula_alba;#005eb0,Fungi;#332288,Glaucophyta;#868510,Haptophyceae;#9f531b,Metamonada;#661100,Metazoa;#6699cc,Rhizaria;#aa4499,Rhodophyta;#ddcc77,Sphaeroforma_arctica;#005eb0,Stramenopiles;#aa4466,Thecamonas_trahens;#005eb0,Viridiplantae;#117733

 ```

# Citation
Thomas A. Richards, Darren M. Soanes, Peter G. Foster, Guy Leonard, Christopher R. Thornton and Nicholas J. Talbot. (2009). Phylogenomic Analysis Demonstrates a Pattern of Rare and Ancient Horizontal Gene Transfer between Plants and Fungi. The Plant Cell. 21(7). www.plantcell.org/cgi/doi/10.1105/tpc.109.065805

# History
This pipeline was first reported in Richards et. al. (2009) as a "an automated gene-by-gene phylogeny pipeline to generate a PhyML tree". Since then it has gone through several iterations and changes, but remains similar in its goal. Between 2009 and 2013 there were quite a few changes in the processing of the script as we had moved to another institution and took development with us. None of these were particularly well documented, other than in internal code comments, although at all times the code was available on request or via github. What you see now is a total rewrite and should be useful to other researchers without too many issues.