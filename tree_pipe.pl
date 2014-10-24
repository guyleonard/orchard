#!/usr/bin/perl
use strict;
use warnings;

#
use autodie;    # bIlujDI' yIchegh()Qo'; yIHegh()!
use Cwd;        # Gets pathname of current working directory
use Digest::MD5;# Generate random string for run ID
use English qw(-no_match_vars);    # No magic perl variables!
use File::Basename;                # Remove path information and extract 8.3 filename
use Getopt::Std;                   # Command line options, finally!
use YAML::XS qw/LoadFile/;         # for the parameters file, user friendly layout

#
use Data::Dumper;                  # temporary during rewrite to dump data nicely to screen

# remove ## comments before publication
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
#     The first version of this program is from 2009.
#     There have been many 'hacks', extensions and updates by GL & FM.
#	  This version attempts to refactor & fix those, mostly GL.
#
#     It is currently hosted here:
#     https://github.com/guyleonard/orchard
#
###########################################################

###########################################################
##           Main Program Flow                           ##
###########################################################

# declare the perl command line flags/options we want to allow
my %options = ();
getopts( "s:t:p:hvbamoq", \%options ) or croak display_help();    # or display_help();

# Display the help message if the user invokes -h
if ( $options{h} ) { display_help() }
if ( $options{v} ) { print "Orchard $VERSION\n"; }

if ( defined $options{p} && defined $options{t} && defined $options{s} ) {
    my $paramaters = LoadFile("$options{p}");
    ## print Dumper($paramaters);

    # read in parameters from YAML file and/or set defaults

    # user options
    my $user_options = $paramaters->{user}->{options} || substr(Digest::MD5::md5_hex(rand), 0, 10);
    my $run_directory = $WORKING_DIR . "/" . $user_options;
    
    # Now we can make the directories
    setup_main_directories($run_directory);
    # Report
    $output_report("Run ID: $user_options\nDirectory: $run_directory\n"); 

    # search options
    my $search_program       = $paramaters->{search}->{program}    || 'blast+';
    my $search_program_blast = $paramaters->{search}->{blast}      || 'blastp';
    my $search_evalue        = $paramaters->{search}->{evalue}     || '1e-10';
    my $search_tophits       = $paramaters->{search}->{top_hits}   || '1';
    my $search_maxlength     = $paramaters->{search}->{max_length} || '3000';
    my $search_special_taxa    = $paramaters->{search}->{special_taxa};        # no default
    my $search_special_tophits = $paramaters->{search}->{special_top_hits};    # no default
    my $search_threads         = $paramaters->{search}->{threads} || '1';

    # alignment options
    my $alignment_program = $paramaters->{alignment}->{program} || 'mafft';
    my $alignment_options = $paramaters->{alignment}->{options};               # no default
    my $alignment_threads = $paramaters->{alignment}->{threads} || '1';

    # masking options
    my $masking_program = $paramaters->{masking}->{program}  || 'trimal';
    my $masking_cutoff1 = $paramaters->{masking}->{cutoff_1} || '50';
    my $masking_cutoff2 = $paramaters->{masking}->{cutoff_2} || '20';

    # tree building options
    my $tree_program = $paramaters->{trees}->{program} || 'FastTreeMP';
    my $tree_options = $paramaters->{trees}->{options};                        # no default
    my $tree_mintaxa = $paramaters->{trees}->{min_taxa} || '4';

    # database options
    my $username  = $paramaters->{database}->{user};
    my $password  = $paramaters->{database}->{password};
    my $server_ip = $paramaters->{database}->{server};
    my $database  = $paramaters->{database}->{database};
    my $tablename = $paramaters->{database}->{tablename};

    # directory options
    my $programs = $paramaters->{directories}->{location} || '/usr/bin';
    my $seq_data = $paramaters->{directories}->{database};                     # no default

    # only run search (blast) step
    if ( $options{b} ) {
        print "Running: Search ($search_program) ONLY\n";
    }

    # only run alignment step
    if ( $options{a} ) {
        print "Running: Alignment ($alignment_program) ONLY\n";
    }

    # only run mask step
    if ( $options{m} ) {
        print "Running: Masking ($masking_program) ONLY\n";
    }

    # only run tree building step (o is for orchard)
    if ( $options{o} ) {
        print "Running: Tree Reconstruction ($tree_program) ONLY\n";
    }

    # run bamt for each sequence sequentially
    if ( $options{q} ) {
        ## experimental, is it really useful?
    }

    if ( !$options{b} && !$options{a} && !$options{m} && !$options{o} && !$options{q} ) {

        # run all steps but all blasts first then amt steps
        print "Running: ALL Steps, all searches ($search_program) first!\n";
    }
}
else {
    display_help();
}

###########################################################
##           Accessory Subroutines                       ##
###########################################################

sub output_report {

    my $message = shift;
    my $file_name = $WORKING_DIR\/report.txt;
    
    open my $report, ">>", $file_name;

    print $report $message;

    return;
}

sub display_help {

    print "You need 3 files for input:\n\t-s sequence(s) file\n\t-t taxa file\n\t-p paramaters file\n";
    print "Example: perl tree_pipe.pl -s sequences.fasta -t taxa_list.txt -p paramaters.yaml\n";
    print
"Other paramaters:\n\t-b blast only\n\t-a alignment only\n\t-m mask only\n\t-t tree building only\n\t-q run sequentially\n";

    exit(1);
}

# this checks to see if the directories needed for file output
# are available, if not it creates them.
sub setup_main_directories {

	my $run_directory = shift;

    # main directories
    my $seqs_dir = "$run_directory/seqs";
    my $alig_dir = "$run_directory/alignments";
    my $mask_dir = "$run_directory/masks";
    my $tree_dir = "$run_directory/trees";
    my $excl_dir = "$run_directory/excluded";
    my $repo_dir = "$run_directory/report";

    # create directories if they don't exist!
    if ( !-d $seqs_dir ) { mkdir $seqs_dir }
    if ( !-d $alig_dir ) { mkdir $alig_dir }
    if ( !-d $mask_dir ) { mkdir $mask_dir }
    if ( !-d $tree_dir ) { mkdir $tree_dir }
    if ( !-d $excl_dir ) { mkdir $excl_dir }
    if ( !-d $repo_dir ) { mkdir $repo_dir }

    return;
}
