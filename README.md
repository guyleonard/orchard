# Orchard Pipeline
A Quick Phylogenetic Tree Building Pipeline...

Seed Sequences --> BLASTp --> Alignment --> Masking --> Phylogenetic Tree Reconstruction --> SVG/PDF Tree Images --> Taxon Annotation

This set of scripts is intended to create a set of "pilot" trees, given your initial gene of interest. The resulting trees can be easily viewed, and then manually improved as needed.

# Dependencies
## orchardDB
Please set up an orchardDB from [here](https://github.com/guyleonard/orchardDB). You may like to use the *cider* DB to start, or build your own.

## Perl Modules
* Bundle::BioPerl
* Bio::DB::Fasta
* DateTime::Format::Duration
* Digest::MD5
* IO::Prompt
* YAML::XS

e.g. you can install these by typing:
```
$ sudo cpanm Bundle::BioPerl Bio::DB::Fasta Digest::MD5 IO::Prompt YAML::XS
```
## Software
### Search Programs
One each, or all of the options below.

1. [Blast+](http://blast.ncbi.nlm.nih.gov/Blast.cgi/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
2. [Blat](https://genome.ucsc.edu/FAQ/FAQblat.html)

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

# Program Execution
## orchard
```
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
```
Usage: orchard_accessories -p params.yaml -n
Required Parameters:
	-p <parameters.yaml>
Other Parameters:
Renaming:
	-s Rename Taxa in Sequence Hits Files
	-a Rename Taxa in Alignment Files
	-m Rename Taxa in Masked Files
	-e Rename Taxa in Excluded Files
	-n Rename Taxa in Newick Trees
```
## Example Parameters File
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

# Citation
Thomas A. Richards, Darren M. Soanes, Peter G. Foster, Guy Leonard, Christopher R. Thornton and Nicholas J. Talbot. (2009). Phylogenomic Analysis Demonstrates a Pattern of Rare and Ancient Horizontal Gene Transfer between Plants and Fungi. The Plant Cell. 21(7). www.plantcell.org/cgi/doi/10.1105/tpc.109.065805

# History
This pipeline was first reported in Richards et. al. (2009) as a "an automated gene-by-gene phylogeny pipeline to generate a PhyML tree". Since then it has gone through several iterations and changes, but remains similar in its goal. Between 2009 and 2013 there were quite a few changes in the processing of the script as we had moved to another institution and took development with us. None of these were particularly well documented, other than in internal code comments, although at all times the code was available on request or via github. What you see now is a total rewrite and should be useful to other researchers without too many issues.