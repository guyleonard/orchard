Darren's Orchard
=======
A Quick Phylogenetic Tree Building Pipeline

Prerequisites
=============
Perl Modules
------------
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

Standalone Programs
-------------------
A choice of either:

1. [Legacy Blast](http://blast.ncbi.nlm.nih.gov/Blast.cgi/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
2. [Blast+](http://blast.ncbi.nlm.nih.gov/Blast.cgi/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
3. [Blat](https://genome.ucsc.edu/FAQ/FAQblat.html)
4. [USEARCH](http://www.drive5.com/usearch/)


Updates
=======
This is going to be a somewhat full re-write of the "'Darren's Orchard'" Pipeline.

I am going to attempt to pull most of the scripts in to the main script.

I am also going to release this version with a Licence. For the moment, before public release, I have include GNU GPL V3.0, but this may change before the first public release...

I am going to try and do away with the reliance on the mysql database. It doesn't do anything that I can't do with a directory of files and the taxdump files from NCBI. It's a pain to continue building it and doesn't help other users getting to use the pipeline. It's not really faster either.