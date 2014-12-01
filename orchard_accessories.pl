#!/usr/bin/perl
use strict;
use warnings;

#
use autodie;    # bIlujDI' yIchegh()Qo'; yIHegh()!
use Bio::DB::Fasta;
use Bio::TreeIO;
use Cwd;                           # Gets pathname of current working directory
use DateTime;                      # Start and End times
use DateTime::Format::Duration;    # Duration of processes
use File::Basename;                # Remove path information and extract 8.3 filename

#use File::Grep qw( fgrep fmap fdo );
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

#our $USER_RETRIEVE = $EMTPY;

###########################################################
##           Main Program Flow                           ##
###########################################################

# declare the perl command line flags/options we want to allow
my %options = ();
getopts( "p:sxnamrcfd", \%options ) or display_help();    # or display_help();

# Display the help message if the user invokes -h
if ( $options{h} ) { display_help() }
if ( $options{v} ) { print "Orchard Accessories $VERSION\n"; }

if ( defined $options{p} ) {

    unless ( $options{s} || $options{x} || $options{n} || $options{a} || $options{m} || $options{r} || $options{c} || $options{f}) {
        print "Missing option: what task do you want to do.\n";
        print display_help();
    }

    # read in parameters from YAML file and/or set defaults
    # user options
    # modify the md5_hex to ten chars from pos 0 if none given in YAML
    my $paramaters = LoadFile("$options{p}");

    $DIR_SEQS    = $paramaters->{directories}->{database};    # no default
    $DIR_TAXDUMP = $paramaters->{directories}->{taxdump};

    #$USER_RETRIEVE = $paramaters->{user}->{retrieve} || 'grep'; #default grep

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

    # SVG Trees
    if ( $options{s} ) {
        print "SVG: Using Dendroscope to Build Trees\n";
        my $start_time = timing('start');
        dendroscope_trees('svg');
        my $end_time = timing( 'end', $start_time );
    }

    # PDFs of Newick Trees
    if ( $options{f} ) {
        print "PDF: Using Dendroscope to Build Trees\n";
        my $start_time = timing('start');
        dendroscope_trees('pdf');
        my $end_time = timing( 'end', $start_time );
    }
}
else {
    display_help();
}

###########################################################
##           SVG Creation Subroutines                    ##
###########################################################

# the original pipeline focused heavily on the use of
# dendroscope, that is continued here, however, I would
# also like to include trees drawn by bioperl and possibly
# other newick to SVG style programs...
# the output should be editable as text and so it *should*
# be trivial as to which program has created the SVG

# folder layout is going to change from the original pipeline
# svg/ is now a folder under the main /trees directory
# pdfs/ will also be a sub-folder to /trees
# renamed trees will be named "renamed" instead of "fixed"
#
# pdfs can be created from original newick and 'renamed'
# newick trees via dendroscope, however, any *renamed_taxonomy
# must be converted with inkscape
#
# ~/ID/trees/
#           svg/
#               accession1_hits_FT.svg
#               accession1_hits_FT_renamed.svg
#               accession1_hits_FT_renamed_taxonomy.svg
#           pdf/
#               accession1_hits_FT.pdf
#               accession1_hits_FT_renamed.pdf
#               accession1_hits_FT_renamed_taxonomy.pdf
#           accession1_hits_FT.tree
#           accession1_hits_FT_renamed.tree

# the nodes are coloured red and the seed accession background
# is coloured green, these are hard set but could be user variables
# in the future

# the Dendroscope manual suggest running this
# xvfb-run --auto-servernum --server-num=1 Dendroscope +g
# to stop the window appearing or so it can be run on a server
# will require the user to install XVFB

sub dendroscope_trees {

    my $type_of_image = shift;

    # get list of all completed tree files
    my $trees_directory  = "$WORKING_DIR\/$USER_RUNID\/trees";
    my @trees_file_names = glob "$trees_directory\/*.tree";

    my $output_directory = "$trees_directory\/$type_of_image";

    if ( !-d $output_directory ) {
        output_report("[INFO]\tCreating $type_of_image Directory: $output_directory\n");
        print "Creating $type_of_image Directory at $output_directory\n";
        mkdir $output_directory;
    }

    foreach my $tree_file (@trees_file_names) {

        my ( $file, $dir, $ext ) = fileparse( $tree_file, '.tree' );

        # strip further filename information to get 'accession' to search in tree
        my $accession = $file;
        $accession =~ s/\_hits\_FT//ism;
        $accession =~ s/\_hits\_renamed//ism;

        # with version 3.0+ of Dendroscope I can't get it to accept this as one string
        # therefore I have to output it to a "command file" and then run the program
        # and then tidy up
        my $dendroscope_execute = "open file=$tree_file\;\n";    # open current tree file
        $dendroscope_execute .= "set drawer=RectangularPhylogram\;\nladderize=left\;\n";                              # set to phylogram tree drawer, ladderise the tree left
        $dendroscope_execute .= "zoom what=expand\;\nset sparselabels=false\;\n";                                     # expand zoom to full tree, show all BS labels
        $dendroscope_execute .= "select nodes=labeled\;\nset labelcolor=255 0 0\;\n";                                 # colour all labels red (includes BS results)
        $dendroscope_execute .= "deselect all\;\nselect nodes=leaves\;\nset labelcolor=0 0 0\;\ndeselect all\;\n";    # colour leaves back to black
        $dendroscope_execute .= "find searchtext=$accession\;\nset labelfillcolor=76 175 80\;\ndeselect all\;\n";        # find accession and colour red
        $dendroscope_execute .= "exportimage file=$output_directory/$file\.$type_of_image format=$type_of_image replace=true\;\n";    # export $type_of_image, overwrite
        $dendroscope_execute .= "quit\;\n";                                                                                           # finish up

        # output command file
        open my $command_file, '>', "$dir\/$file\.cmd";
        print $command_file $dendroscope_execute;
        close $command_file;

        # run dendroscope with previous command file
        my $dendroscope_command = "Dendroscope -g true";                                                                              # run without GUI
        $dendroscope_command .= " -c $dir\/$file\.cmd";                                                                               # execute following command
        system($dendroscope_command);

        # delete command file - consider moving to report?
        unlink "$dir\/$file\.cmd";
    }

    return;
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

    #my ( $file, $dir, $ext ) = fileparse $aligned_sequences, '\.tree';

    print Dumper @trees_file_names;

    foreach my $tree_file (@trees_file_names) {

        my $treeio = Bio::TreeIO->new(
            -format => 'newick',
            -file   => $tree_file
        );

        #my $tree = $input->next_tree;
        while ( my $tree = $treeio->next_tree ) {

            #my @tree_leafs = $tree->get_leaf_nodes;
            for my $node ( grep { $_->is_Leaf } $tree->get_nodes ) {

                print "Node is " . $node->id . "\t";
                my @new_node = retrieve_taxa_name( $node->id );
                print "from @new_node\n";
            }

            #print Dumper @tree_leafs;
        }
    }

    #retrieve_taxa_name_bioperl();

    return;
}

# Jeez this is also really really slow
sub retrieve_taxa_name {

    my $search_node  = shift;
    my @search_files = glob "$DIR_SEQS\/*.fas";

    #my $return_node = $EMPTY;

    my @return_node = system "grep $search_node $DIR_SEQS\/*.fas";

    #my $return_node = grep {$_ eq $search_node} @search_files;

    return @return_node;

}

# so I am going to do this with BioPerl at first
# but I may also try grep which may or may not be faster
# might abandom this method - it takes 5 hours to build an
# index and then a long time to load it back in to memory
sub retrieve_taxa_name_bioperl {

    #print "$DIR_SEQS and $DIR_TAXDUMP\n";

    # may need to warn users that a large DB will take some time to index
    # on the first run, if it hasn't already been done
    # check for 'directory.index'
    # with 1110 genomes it took 5 hours to index!

    unless ( -e "$DIR_SEQS\/directory.index" ) { print "Indexing Directory: This may take some time (1K genomes ~ 5 hours)!\n" }

    # read in directory of '.fas' files and create an index
    # 'glob' => '*.{fa,FA,fasta,FASTA,fast,FAST,dna,DNA,fna,FNA,faa,FAA,fsa,FSA}',
    # we have to modify it to include '.fas' as it's not a 'standard' FASTA extension
    # apparently
    my $sequences_db = Bio::DB::Fasta->new( "$DIR_SEQS", -glob => "*.fas" );
    my @ids = $sequences_db->get_all_primary_ids;

    #print Dumper $sequences_db;

    return;
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
    print "SVG Trees:\n\t-s Build SVG Trees (requires Dendroscope)\n\t-x Build SVG Trees (requires BioPerl)\n";
    print "Renaming:\n\t-n Rename taxa in newick trees\n\t-a Rename taxa in unmasked alignments\n\t-m Rename taxa in masked alignments\n\t-r Rename taxa in SVG trees\n";
    print "PDF Trees:\n\t-f Build PDF Trees of all newick trees (requires Dendroscope)\n\t-d Build PDFs of annotated SVG trees (requires Inkscape)\n";
    print "Annotating Taxonomy:\n\t-c Colourise taxon names in SVG trees\n";
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
