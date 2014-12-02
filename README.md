# Darren's Orchard

A Quick Phylogenetic Tree Building Pipeline

# Prerequisites

There are quite a few programs and other prerequisites you will need to have installed in your *environment path*, i.e. you will need to be able to call a program from the commandline without specifying the directory to run it from.

## Perl Modules

* Bundle::Bioperl
* Bio::DB::Fasta (sometimes I find this isn't installed with BioPerl !?)
* Digest::MD5
* IO::Prompt
* YAML::XS

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
4. [USEARCH](http://www.drive5.com/usearch/)

#### Alignment Programs

1. [MAFFT](http://mafft.cbrc.jp/alignment/software/)
2. [MUSCLE](http://www.drive5.com/muscle/)

#### Masking Programs

1. [trimAl v1.3](http://trimal.cgenomics.org/)
2. ~~[Gblocks](http://molevol.cmima.csic.es/castresana/Gblocks.html)~~ *

* deprecated and no longer supported. trimAl FTW!

#### Tree Programs

1. [FastTree2.1](http://meta.microbesonline.org/fasttree/)
2. ~~[PhyML 3.0](http://www.atgc-montpellier.fr/phyml/binaries.php)~~ *

* deprecated and no longer supported. FastTree FTW!

### orchard_accessories.pl

#### Renaming Taxa

1. [NCBI Taxonomy](ftp://ftp.ncbi.nih.gov/pub/taxonomy) - FTP retrieve 'taxdump.tar.gz' from ftp://ftp.ncbi.nih.gov/pub/taxonomy

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

# Updates

This is going to be a somewhat full re-write of the "'Darren's Orchard'" Pipeline.

I have decided to split the pipeline in to two scripts, the main script "orchard.pl" will control the creation of the trees via blast searches, alignment and masking, and an accessory script "orchard_accessories.pl" which will handle all the renaming of taxon IDs to genus/species, the creation of SVG trees and the colorising and renaming of those SVGs. I have also created a parameters file option in YAML - this allows for a very easy to read settings file for the user and also allows to record better the options used on each run!

I am going to try and do away with the reliance on the mysql database. It doesn't do anything that I can't do with a directory of files and the taxdump files from NCBI. It's a pain to continue building it and doesn't help other users getting to use the pipeline. It's not really faster anymore either. It will rely on the user having up to date copies of the taxdump files though...I will also still support the 'remote' option.

I am also going to release this version with a Licence. For the moment, before public release, I have include GNU GPL V3.0, but this may change before the first public release...