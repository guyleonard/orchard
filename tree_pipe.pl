#!/usr/bin/perl
use strict;
use warnings;

#
use autodie;                       # bIlujDI' yIchegh()Qo'; yIHegh()!
use Cwd;                           # Gets pathname of current working directory
use English qw(-no_match_vars);    # No magic perl variables!
use File::Basename;                # Remove path information and extract 8.3 filename
use Getopt::Std;                   # Command line options, finally!
use YAML::XS qw/LoadFile/;         # for the parameters file, user friendly layout

#
our $WORKING_DIR = getcwd;
our $VERSION     = '2014-10-17';
###########################################################
#           Darren's Orchard Pipeline                     #
###########################################################
#
#     A Quick Phylogenetic Tree Building Pipeline
#
#     By	-	Darren Soannes & Thomas Richards
#               Guy Leonard
#               Finlay Maguire
#
#     This first version of this program was in 2009.
#     There have been many 'hacks', extensions and updates.
#	  This version attempts to refactor & fix those.
#
#     It is currently hosted here:
#     https://github.com/guyleonard/orchard
#
###########################################################

###########################################################
##           Main Program Flow                           ##
###########################################################
setup_main_directories();

# declare the perl command line flags/options we want to allow
my %options = ();
getopts( "s:t:p:hbv", \%options ) or display_help();

# Display the help message if the user invokes -h
if    ( $options{h} ) { display_help() }
elsif ( $options{v} ) { print "Orchard $VERSION\n"; }

###########################################################
##           Accessory Subroutines                       ##
###########################################################

sub display_help {

    print "You need 3 files for input:\n\t-s sequence(s) file\n\t-t taxa file\n\t-p paramaters file\n";
    print "Example: perl tree_pipe.pl -s sequences.fasta -t taxa_list.txt -p paramaters.yaml\n";
    print "Other paramaters:\n\t-b blast only\n\t-a alignment only\n\t-m mask only\n\t-t tree building only\n";

    return;
}

# this checks to see if the directories needed for file output
# are available, if not it creates them.
sub setup_main_directories {

    # main directories
    my $seqs_dir = "$WORKING_DIR/seqs";
    my $alig_dir = "$WORKING_DIR/alignments";
    my $mask_dir = "$WORKING_DIR/masks";
    my $tree_dir = "$WORKING_DIR/trees";
    my $excl_dir = "$WORKING_DIR/excluded";
    my $repo_dir = "$WORKING_DIR/report";

    # create directories if they don't exist!
    if ( !-d $seqs_dir ) { mkdir $seqs_dir }
    if ( !-d $alig_dir ) { mkdir $alig_dir }
    if ( !-d $mask_dir ) { mkdir $mask_dir }
    if ( !-d $tree_dir ) { mkdir $tree_dir }
    if ( !-d $excl_dir ) { mkdir $excl_dir }
    if ( !-d $repo_dir ) { mkdir $repo_dir }

    return;
}
