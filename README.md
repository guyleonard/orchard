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

One each of the following (our preference is always 1.):

### Search Programs

1. [Blast+](http://blast.ncbi.nlm.nih.gov/Blast.cgi/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
2. [Legacy Blast](http://blast.ncbi.nlm.nih.gov/Blast.cgi/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
3. [Blat](https://genome.ucsc.edu/FAQ/FAQblat.html)
4. [USEARCH](http://www.drive5.com/usearch/)

### Alignment Programs

1. [MAFFT](http://mafft.cbrc.jp/alignment/software/)
2. [MUSCLE](http://www.drive5.com/muscle/)

### Masking Programs

1. [trimAl v1.3](http://trimal.cgenomics.org/)
2. ~~[Gblocks](http://molevol.cmima.csic.es/castresana/Gblocks.html)~~ *

* deprecated and no longer supported. trimAl FTW!

### Tree Programs

1. [FastTree2.1](http://meta.microbesonline.org/fasttree/)
2. ~~[PhyML 3.0](http://www.atgc-montpellier.fr/phyml/binaries.php)~~ *

* deprecated and no longer supported. FastTree FTW!

# Program Execution

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
	-t tree building only
	-q run sequentially
```

Order preference of parameters is ignored: i.e. specificy '-ab' will run blasts then alignments.

# Updates

This is going to be a somewhat full re-write of the "'Darren's Orchard'" Pipeline.

I am going to attempt to pull most of the scripts in to the main script.

I am also going to release this version with a Licence. For the moment, before public release, I have include GNU GPL V3.0, but this may change before the first public release...

I am going to try and do away with the reliance on the mysql database. It doesn't do anything that I can't do with a directory of files and the taxdump files from NCBI. It's a pain to continue building it and doesn't help other users getting to use the pipeline. It's not really faster either.