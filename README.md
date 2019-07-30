# Orchard Pipeline
A Quick Phylogenetic Tree Building Pipeline...

Seed Sequences --> BLASTp --> Alignment --> Masking --> Phylogenetic Tree Reconstruction --> SVG/PDF Tree Images --> Taxon Annotation

This set of scripts is intended to create a set of "pilot" trees, given your initial gene of interest. The resulting trees can be easily viewed, and then manually improved as needed.

# Dependencies
## orchardDB
Please set up an orchardDB from [here](https://github.com/guyleonard/orchardDB). You may like to use the *cider* DB to start or build your own.

## Perl Modules
* Bundle::BioPerl
* Bio::DB::Fasta
* DateTime::Format::Duration
* Digest::MD5
* IO::Prompt
* YAML::XS
You can install these by typing:
```
$ sudo cpanm Bundle::BioPerl Bio::DB::Fasta Digest::MD5 IO::Prompt YAML::XS
```
## Programs
### Search Programs
1. [Blast+](http://blast.ncbi.nlm.nih.gov/Blast.cgi/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
2. [Blat](https://genome.ucsc.edu/FAQ/FAQblat.html)

#### Alignment Programs
1. [MAFFT](http://mafft.cbrc.jp/alignment/software/)
2. [MUSCLE](http://www.drive5.com/muscle/)

#### Masking Programs
1. [trimAl v1.3](http://trimal.cgenomics.org/)
2. [Divvier](https://github.com/simonwhelan/Divvier) - Coming Soon.

#### Tree Building Programs
1. [FastTree2.1](http://meta.microbesonline.org/fasttree/)
2. [IQTree](https://github.com/Cibiv/IQ-TREE) - Coming Soon.

#### Renaming Taxa
1. [NCBI Taxonomy](ftp://ftp.ncbi.nih.gov/pub/taxonomy) - retrieve 'taxdump.tar.gz' from ftp://ftp.ncbi.nih.gov/pub/taxonomy

#### SVG Tree Creation
1. [Dendroscope 3](http://ab.inf.uni-tuebingen.de/software/dendroscope/)
  1. XVFB - If you are running this on a server, you will not have a graphics interface which Dendroscope needs (even in command line mode), so please install XVFB.
  2. e.g. sudo apt-get install xvfb
  3. Dendroscope is run: xvfb-run --auto-servernum --server-num=1 Dendroscope +g

#### PDF Creation (SVG to PDF)
1. [Inkscape](https://www.inkscape.org/en/)

# Program Execution
## orchard.pl
```
Required files for input:
	-s sequence(s) file
	-t taxa file
	-p paramaters file
Example: perl orchard.pl -s sequences.fasta -t taxa_list.txt -p paramaters.yaml
Other paramaters:
	-b blast only
	-a alignment only
	-m mask only
	-o tree building only
```
Order preference of parameters is ignored: i.e. specificy '-ab' will run blasts then alignments.

-s can also be a folder of multiple files with multiple fasta sequences in each, which can be used for either a) importing a bunch of files which only have one sequence in each (although you would be better concatenating them to one FASTA file and using the normal approach) or b) multiple related sequences grouped or clustered where you wish for only one tree to be generated but from all the blast hits to each seed sequence in one file - for example, we use this to look at orthologous groups of genes.

## orchard_accessories.pl
```
Required input:
	-p parameters file
Other parameters:
SVG Trees:
	-s Build SVG Trees (requires Dendroscope)
	-x Build SVG Trees (requires BioPerl)
Renaming:
	-n Rename taxa in newick trees
	-a Rename taxa in unmasked alignments
	-m Rename taxa in masked alignments
	-r Rename taxa in SVG trees
PDF Trees:
	-f Build PDF Trees of all newick trees (requires Dendroscope)
	-d Build PDFs of annotated SVG trees (requires Inkscape)
Annotating Taxonomy:
	-c Colourise taxon names in SVG trees
```
Order preference of parameters is ignored but you cannot colourise trees if you have not created them first!

## Example Parameters File
```yaml
user:
 run_id: test_ortho # give your run an ID
 reindex: n # force reindex of search db y/n
 retrieve: bioperl # bioperl or grep

search:
 program: blast+ # blast, blast+, blat, usearch
 subprogram: blastp # blastp, blastx, ublast, blat
 threads: 4 # it will be 1 if you don't specify
 ## values below are only supported by blast/blast+/usearch
 evalue: 1e-10
 top_hits: 5
 ## values below are only supported by blast/blast+
 max_length: 3000
 special_taxa: # must be a comma separated list of taxa names or blank
 special_top_hits: 5

alignment:
 program: mafft  # mafft, muscle
 # mafft = --auto --quiet --reorder
 # muscle = -maxiters 2 -quiet -group
 options: --auto --quiet --reorder
 threads: 4

masking:
 cutoff_1: 50
 cutoff_2: 30

trees:
 program: FastTreeMP # FastTree, FastTreeMP
 # e.g. -bionj -slow
 options: -bionj -slow -quiet
 min_taxa: 5 # min number of taxa to make a tree from
 node_colour: 255 0 0 # r g b (only, no hex)
 seed_colour: 76 175 80 # r g b (only, no hex)

directories:
 # the location of the *.fasta predicted protein files
 database: /home/genomes/cider
 taxdump: /home/genomes/taxonomy
 ```

# Citation
Thomas A. Richards, Darren M. Soanes, Peter G. Foster, Guy Leonard, Christopher R. Thornton and Nicholas J. Talbot. (2009). Phylogenomic Analysis Demonstrates a Pattern of Rare and Ancient Horizontal Gene Transfer between Plants and Fungi. The Plant Cell. 21(7). www.plantcell.org/cgi/doi/10.1105/tpc.109.065805

# History
This pipeline was first reported in Richards et. al. (2009) as a "an automated gene-by-gene phylogeny pipeline to generate a PhyML tree" and since then has gone through several iterations, but remains similar in its goal. Between 2009 and 2013 there were quite a few changes in the processing of the script as we had moved to another institution and took development with us, these included; new taxa added to the database along with a small redesign of the initial database, new programs added e.g. trimAl and FastTree2 instead of GBLocks and PhyML, some additional scripts to handle renaming of taxa, annotation and graphically addding PFAM domains, along with plenty of other code changes. None of these were particularly well documented, other than in internal code comments, although at all times the code were available on request or via github and it became quite messy to clone and to get other people to use.

This is an attempt to do better and make available a newer version of this process.

NB - In many of our publications we have used an accessory script to plot PFAM domains on to our trees at the relevant taxa leafs, this script is unfortunately not available in this version as it is under copyright of [Bill Wickstead](http://www.wicksteadlab.co.uk/).