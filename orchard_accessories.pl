#!/usr/bin/perl
use strict;
use warnings;

#
use autodie;    # bIlujDI' yIchegh()Qo'; yIHegh()!
use Bio::DB::Fasta;
use Cwd;                           # Gets pathname of current working directory
use DateTime;                      # Start and End times
use DateTime::Format::Duration;    # Duration of processes
use File::Basename;                # Remove path information and extract 8.3 filename
use Getopt::Std;                   # Command line options, finally!
use feature qw{ switch };
no if $] >= 5.017011, warnings => 'experimental::smartmatch';    # ignore experimental warning for 'when'
use IO::Prompt;                                                  # User prompts
use YAML::XS qw/LoadFile/;                                       # for the parameters file, user friendly layout

##
use Data::Dumper;                                                # temporary during rewrite to dump data nicely to screen

our $WORKING_DIR = getcwd();
our $VERSION     = '2014-10-17';
###########################################################
#           Darren's Orchard Pipeline - Accessories       #
###########################################################
#
#     A Quick Phylogenetic Tree Building Pipeline
#
#     By	-	Darren Soannes & Thomas Richards
#               Guy Leonard
#               Finlay Maguire
#
#     These are some tools to edit, rename & annotate
#     trees produced by Darren's Orchard.
#
#     It is currently hosted here:
#     https://github.com/guyleonard/orchard
#
#     Phylogenomic Analysis Demonstrates a Pattern of Rare
#     and Ancient Horizontal Gene Transfer between
#     Plants and Fungi
#     doi: 10.1105/tpc.109.065805
#
###########################################################

###########################################################
##           Global Variables                            ##
###########################################################
# These options are global and will be set from the user YAML
# file read in below, globals to avoid passing multiple values
# to sub routines, they should not be edited once set.
our $EMPTY       = q{};
our $USER_RUNID  = $EMPTY;
our $DIR_SEQS    = $EMPTY;
our $DIR_TAXDUMP = $EMPTY;

###########################################################
##           Main Program Flow                           ##
###########################################################

# declare the perl command line flags/options we want to allow
my %options = ();
getopts( "p:sxnamrc", \%options ) or display_help();    # or display_help();

# Display the help message if the user invokes -h
if ( $options{h} ) { display_help() }
if ( $options{v} ) { print "Orchard Accessories $VERSION\n"; }

if ( defined $options{p} ) {

    unless ($options{s} || $options{x} || $options{n}) {
        print "Missing option: what task do you want to do.\n";
        print display_help();
    }

    # read in parameters from YAML file and/or set defaults
    # user options
    # modify the md5_hex to ten chars from pos 0 if none given in YAML
    my $paramaters = LoadFile("$options{p}");

    $DIR_SEQS    = $paramaters->{directories}->{database};    # no default
    $DIR_TAXDUMP = $paramaters->{directories}->{taxdump};

    # Tell the user they haven't set the taxdump directory properly
    # Ideally we should also check if the directory exists...
    if ( $DIR_TAXDUMP eq $EMPTY ) { print "User must specify taxdump directory!\n"; display_help(); }

    # we don't want to generate a random ID here, the user must have one previously!
    # so therefore let's check and output help if not...
    $USER_RUNID = $paramaters->{user}->{run_id};
    if ( $USER_RUNID eq $EMPTY ) { print "User must specify a USER ID in paramaters files!\n"; display_help(); }

    my $run_directory = "$WORKING_DIR\/$USER_RUNID";

    # Now we can make the directories
    # setup_main_directories( $run_directory, $options{f} );

    # Report
    output_report("[INFO]\tRun ID: $USER_RUNID\n[INFO]\tDirectory: $run_directory\n");

    # directory options
    $DIR_SEQS = $paramaters->{directories}->{database};    # no default

    # rename newick trees
    if ( $options{n} ) {
        print "Renaming: Newick Trees\n";
        my $start_time = timing('start');
        rename_newick_trees();
        my $end_time = timing( 'end', $start_time );
    }
}
else {
    display_help();
}

###########################################################
##           Renaming Subroutines                        ##
###########################################################

# In order to rename any of the non-standard accessions
# used throughout the 'orchard' we will first have to lookup
# which taxa file they are located in, then the file name is
# the taxonomy
#
# for colouring based on taxonomy then go and lookup
# the taxon ID and then the taxonomy information from NCBI
#
# previously we could query the mysql database and get this
# information all together but the previous steps had to be
# done in the creation of the mysql db, and as I have alluded
# to elsewhere, that was getting tiresome to build, update and
# add new genomes to, compared to just adding a file to a folder

sub rename_newick_trees {

    # directory variables
    my $mask_directory  = "$WORKING_DIR\/$USER_RUNID\/mask";
    my $trees_directory = "$WORKING_DIR\/$USER_RUNID\/trees";

    # get list of all completed tree files
    my @trees_file_names = glob "$trees_directory\/*.tree";

    print Dumper @trees_file_names;

    retrieve_taxa_name();
}

# so I am going to do this with BioPerl at first
# but I may also try grep which may or may not be faster
sub retrieve_taxa_name {

    print "$DIR_SEQS and $DIR_TAXDUMP\n";

    # may need to warn users that a large DB will take some time to index
    # on the first run, if it hasn't already been done
    # check for 'directory.index'

    unless ( -e "$DIR_SEQS\/directory.index") { print "Indexing Directory: This may take some time!\n"}

    # read in directory of '.fas' files and create an index
    # 'glob' => '*.{fa,FA,fasta,FASTA,fast,FAST,dna,DNA,fna,FNA,faa,FAA,fsa,FSA}',
    # we have to modify it to include '.fas' as it's not a 'standard' FASTA extension
    # apparently
    my $sequences_db = Bio::DB::Fasta->new("$DIR_SEQS", -glob => "*.fas");
    my @ids          = $sequences_db->get_all_primary_ids;

    print Dumper $sequences_db;

}

###########################################################
##           Accessory Subroutines                       ##
###########################################################

sub output_report {

    # Append messages to report file
    my $message   = shift;
    my $file_name = "$WORKING_DIR\/$USER_RUNID\_accessory\_report.txt";
    open my $report, ">>", $file_name;
    print $report $message;
    close($report);
    return;
}

sub display_help {

    print "Required input:\n\t-p parameters file\n";
    print "Other parameters:\n";
    print "SVG Trees:\n\t-s Build SVG Trees (requires Dendroscope)\n\t-x Build SVG Trees (no Dendroscope)\n";
    print "Renaming:\n\t-n Rename taxa in newick trees\n\t-a Rename taxa in unmasked alignments\n\t-m Rename taxa in masked alignments\n\t-r Rename taxa in SVG trees\n";
    print "Colouring Taxonomy:\n\t-c Colourise taxon names in SVG trees\n";
    exit(1);
}

# this will need some reworking as the directories should already be set up
# and depending on which option you choose you won't need some of them anyway
sub setup_main_directories {

    my $run_directory = shift;
    my $force         = shift;

    # main directories
    my $seqs_dir = "$run_directory\/seqs";
    my $alig_dir = "$run_directory\/alignments";
    my $mask_dir = "$run_directory\/masks";
    my $tree_dir = "$run_directory\/trees";
    my $excl_dir = "$run_directory\/excluded";
    my $repo_dir = "$run_directory\/report";

    #if ( $force != 0 ) {
    if ( !-d $run_directory ) {
        output_report("[INFO]\tCreating Directory: $run_directory\n");
        mkdir $run_directory;

        # create sub-directories
        output_report("[INFO]\tCreating Subdirectories\n");

        # create directories if they don't exist!
        if ( !-d $seqs_dir ) { mkdir $seqs_dir }
        if ( !-d $alig_dir ) { mkdir $alig_dir }
        if ( !-d $mask_dir ) { mkdir $mask_dir }
        if ( !-d $tree_dir ) { mkdir $tree_dir }
        if ( !-d $excl_dir ) { mkdir $excl_dir }
        if ( !-d $repo_dir ) { mkdir $repo_dir }
    }
    else {
        print "Directory Already Exists!\nContinue anyway? (this may overwrite files) y/n\n";
        my $user_choice = prompt( " > : ", -yes_no1 );
        if ( $user_choice =~ m/n/ism ) {

            output_report("[INFO]\tTerminating: run directories already exist\n");
            print "Bye!\n";
            exit;
        }
        else {
            output_report("[WARN]\tContinuing: even though run directories already exist\n");

            # Let's just make sure we have everything we will need!

            # create directories if they don't exist!
            if ( !-d $seqs_dir ) { mkdir $seqs_dir }
            if ( !-d $alig_dir ) { mkdir $alig_dir }
            if ( !-d $mask_dir ) { mkdir $mask_dir }
            if ( !-d $tree_dir ) { mkdir $tree_dir }
            if ( !-d $excl_dir ) { mkdir $excl_dir }
            if ( !-d $repo_dir ) { mkdir $repo_dir }
        }
    }

    #}
    return;
}

# Handles output of start/end and duration times of different steps
sub timing {

    my $time_operation = shift;
    my $start_time     = shift;

    my $time_now = $EMPTY;

    if ( $time_operation eq 'start' ) {

        $time_now = DateTime->now;
        print "Start Time: $time_now\n";
        output_report("[INFO]\tStart Time: $time_now\n");
    }
    else {
        $time_now = DateTime->now;
        print "End Time: $time_now\n";
        output_report("[INFO]\tEnd Time: $time_now\n");
        my $total_time = $time_now->subtract_datetime_absolute($start_time);
        my $format = DateTime::Format::Duration->new( pattern => '%e days, %H hours, %M minutes, %S seconds' );
        print "Elapsed Time: " . $format->format_duration($total_time) . "\n";
        output_report( "[INFO]\tElapsed Time: " . $format->format_duration($total_time) . "\n" );
    }

    return $time_now;
}
