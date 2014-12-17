# Darren's Orchard

A Quick Phylogenetic Tree Building Pipeline...

Seed Sequences --> BLASTp --> Alignment --> Masking --> Phylogenetic Tree Reconstruction --> SVG/PDF Tree Images --> Taxon Annotation

This pipeline was first reported in Richards et. al. (2009) as a "an automated gene-by-gene phylogeny pipeline to generate a PhyML tree" and since then has gone through several iterations, but remains similar in its goal. Between 2009 and 2013 there were quite a few changes in the processing of the script as we had moved to another institution and took development with us, these included; new taxa added to the database along with a small redesign of the initial database, new programs added e.g. trimAl and FastTree2 instead of GBLocks and PhyML, some additional scripts to handle renaming of taxa, annotation and graphically addding PFAM domains, along with plenty of other code changes. None of these were particularly well documented, other than in internal code comments, although at all times the code were available on request and it became quite messy to clone and to get other people to use.

This is an attempt to do better and make available a newer version of this process. It won't be entirely modular with easy switch in/out of other programs in the pipeline but it should be in a state where any other user can clone the repository, set up a databse of their own genomes, install the relevant programs, and start making 100s of phylogenetic guide trees very easily. An attempt has been made to reduce the number of accessory scripts, simplify the process of supplying options - by using YAML (human readable parameter files) and command line options instead of having to edit any code!

It should also be noted that [Finlay Maguire](https://github.com/fmaguire) 'forked' the initial code but ended up rewriting it in Python, when he grew his hipster beard, in an attempt to parallelise many of the steps, but still relying on the original database, which you may prefer to use.

NB - In many of our publications we have used an accessory script to plot PFAM domains on to our trees at the relevant taxa leafs, this script is unfortunately not available in this version as it is under copyright of [Bill Wickstead](http://www.wicksteadlab.co.uk/).

# Prerequisites

There are quite a few programs and other prerequisites you will need to have installed in your *environment path*, i.e. you will need to be able to call a program from the commandline without specifying the directory to run it from.

## Perl Modules

* Bundle::Bioperl
* Bio::DB::Fasta (sometimes I find this isn't installed with BioPerl !?)
* Digest::MD5
* IO::Prompt
* YAML::XS
* Perl4::CoreLibs ???

You can install these by typing:

```
$ sudo -MCPAN -e shell
cpan[1]> install prerequisit_name
```

## Standalone Programs

One each of the following, our preference is always #1 and if you wish to add your own to the line up then please feel free to do so and submit a pull request:

### orchard.pl

#### Search Programs

1. [Blast+](http://blast.ncbi.nlm.nih.gov/Blast.cgi/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
2. [Legacy Blast](http://blast.ncbi.nlm.nih.gov/Blast.cgi/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
3. [Blat](https://genome.ucsc.edu/FAQ/FAQblat.html)
4. [USEARCH](http://www.drive5.com/usearch/) - ublast only

~~X. [VSEARCH](https://github.com/torognes/vsearch) - usearch_global~~ No amino acid support :(


#### Alignment Programs

1. [MAFFT](http://mafft.cbrc.jp/alignment/software/)
2. [MUSCLE](http://www.drive5.com/muscle/)

#### Masking Programs

1. [trimAl v1.3](http://trimal.cgenomics.org/)
2. ~~[Gblocks](http://molevol.cmima.csic.es/castresana/Gblocks.html)~~ no longer supported. trimAl FTW!

#### Tree Programs

1. [FastTree2.1](http://meta.microbesonline.org/fasttree/)
2. ~~[PhyML 3.0](http://www.atgc-montpellier.fr/phyml/binaries.php)~~ no longer supported. FastTree FTW!

### orchard_accessories.pl

#### Renaming Taxa

1. [NCBI Taxonomy](ftp://ftp.ncbi.nih.gov/pub/taxonomy) - retrieve 'taxdump.tar.gz' from ftp://ftp.ncbi.nih.gov/pub/taxonomy and put expand it to a central directory (you will then add the lcoation to your parameters file)

#### SVG Tree Creation

1. [Dendroscope 3](http://ab.inf.uni-tuebingen.de/software/dendroscope/)

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

# Updates

This is going to be a somewhat full re-write of the "'Darren's Orchard'" Pipeline.

I have decided to split the pipeline in to two scripts, the main script "orchard.pl" will control the creation of the trees via blast searches, alignment and masking, and an accessory script "orchard_accessories.pl" which will handle all the renaming of taxon IDs to genus/species, the creation of SVG trees and the colorising and renaming of those SVGs. I have also created a parameters file option in YAML - this allows for a very easy to read settings file for the user and also allows to record better the options used on each run!

I am going to try and do away with the reliance on the mysql database. It doesn't do anything that I can't do with a directory of files and the taxdump files from NCBI. It's a pain to continue building it and doesn't help other users getting to use the pipeline. It's not really faster anymore either. It will rely on the user having up to date copies of the taxdump files though...I will also still support the 'remote' option.

I am also going to release this version with a Licence. For the moment, before public release, I have include GNU GPL V3.0, but this may change before the first public release...

# Citation

Thomas A. Richards, Darren M. Soanes, Peter G. Foster, Guy Leonard, Christopher R. Thornton and Nicholas J. Talbot. (2009). Phylogenomic Analysis Demonstrates a Pattern of Rare and Ancient Horizontal Gene Transfer between Plants and Fungi. The Plant Cell. 21(7). www.plantcell.org/cgi/doi/10.1105/tpc.109.065805
