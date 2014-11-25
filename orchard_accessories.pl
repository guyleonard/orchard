#!/usr/bin/perl
use strict;
use warnings;

#
use autodie;                       # bIlujDI' yIchegh()Qo'; yIHegh()!
use Cwd;                           # Gets pathname of current working directory
use DateTime;                      # Start and End times
use DateTime::Format::Duration;    # Duration of processes
use Digest::MD5;                   # Generate random string for run ID
use File::Basename;                # Remove path information and extract 8.3 filename
use Getopt::Std;                   # Command line options, finally!
use feature qw{ switch };
no if $] >= 5.017011, warnings => 'experimental::smartmatch';    # ignore experimental warning for 'when'
use IO::Prompt;                                                  # User prompts
use YAML::XS qw/LoadFile/;                                       # for the parameters file, user friendly layout

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
###########################################################

###########################################################
##           Global Variables                            ##
###########################################################
# These options are global and will be set from the user YAML
# file read in below, globals to avoid passing multiple values
# to sub routines, they should not be edited once set.
our $EMPTY             = q{};
our $USER_RUNID   = $EMPTY;

###########################################################
##           Main Program Flow                           ##
###########################################################

# declare the perl command line flags/options we want to allow
my %options = ();
getopts( "s:t:p:hvbamof", \%options ) or croak display_help();    # or display_help();

# Display the help message if the user invokes -h
if ( $options{h} ) { display_help() }
if ( $options{v} ) { print "Orchard $VERSION\n"; }


###########################################################
##           Accessory Subroutines                       ##
###########################################################

sub output_report {

    # Append messages to report file
    my $message   = shift;
    my $file_name = "$WORKING_DIR\/$USER_RUNID\_report.txt";
    open my $report, ">>", $file_name;
    print $report $message;
    close($report);
    return;
}

sub display_help {

    print "Required files for input:\n\t-s sequence(s) file\n\t-t taxa file\n\t-p paramaters file\n";
    print "Example: perl orchard.pl -s sequences.fasta -t taxa_list.txt -p paramaters.yaml\n";
    print "Other paramaters:\n\t-b blast only\n\t-a alignment only\n\t-m mask only\n\t-o tree building only\n\t-q run sequentially\n\t-f force yes\n";
    exit(1);
}

# this checks to see if the directories needed for file output
# are available, if not it creates them all - for now.
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