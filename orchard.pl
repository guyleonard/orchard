#!/usr/bin/perl
use strict;
use warnings;

use 5.012;
no warnings 'experimental::smartmatch';    # ignore experimental warning for given/when
use autodie;                               # bIlujDI' yIchegh()Qo'; yIHegh()!
use Bio::DB::Fasta;                        # for indexing fasta files for sequence retrieval
use Bio::SearchIO;                         # parse blast tabular format
use Bio::SeqIO;                            # Use Bio::Perl
use Bio::Taxon;
use Bio::Tree::Tree;
use Cwd;                                   # Gets pathname of current working directory
use DateTime::Format::Duration;            # Duration of processes
use DateTime;                              # Start and End times
use DBI;                                   # mysql database access
use Digest::MD5;                           # Generate random string for run ID
use Digest::MD5 'md5_hex';
use File::Basename;                        # Remove path information and extract 8.3 filename
use Getopt::Std;                           # Command line options, finally!
use IO::Prompt;                            # User prompts
use YAML::XS qw/LoadFile/;                 # for the parameters file, user friendly layout
our $VERSION = '2017-06-26';

##
use Data::Dumper;                          # temporary during rewrite to dump data nicely to screen
## remove ## comments before publication

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
my $Working_Dir = getcwd();

my $Alignment_Options;
my $Alignment_Program;
my $Alignment_Threads;

my $Masking_Cutoff1;
my $Masking_Cutoff2;

my $Search_Evalue = '1e-10';
my $Search_Maxlength = '3000';
my $Search_Program = 'blast+';
my $Search_Subprogram = 'blastp';
my $Search_Threads = '1';
my $Search_Tophits = '1';
my @Search_Special_Taxa;
my $Search_Special_Tophits = '1';

my $Seq_Data;

my $Tree_Mintaxa;
my $Tree_Options;
my $Tree_Program;
my $Tree_Model;
my $Tree_Bs;

my $User_Reindex;
my $User_Runid;

my $Filter_File;


###########################################################
##           Main Program Flow                           ##
###########################################################
# declare the perl command line flags/options we want to allow
my %options = ();
getopts( 's:t:p:hvbfamo', \%options ) or display_help();

# Display the help message if the user invokes -h
if ( $options{h} ) { display_help(); }

if ( $options{v} ) { print "Orchard $VERSION\n"; }

if ( defined $options{s} && defined $options{t} && defined $options{p} ) {
    my $taxa_list_fname = "$options{t}";
    my @taxa_array      = taxa_list_to_array($taxa_list_fname);

    my $input_seqs_fname = "$options{s}";

    # read in parameters from YAML file overiding set defaults
    # user options
    # modify the md5_hex to ten chars from pos 0 if none given in YAML
    my $paramaters = LoadFile("$options{p}");
    $User_Runid = $paramaters->{user}->{run_id} || substr( Digest::MD5::md5_hex(rand), 0, 10 );
    $User_Reindex = $paramaters->{user}->{reindex} || 'n';    # default no
    if ( $User_Reindex eq 'y' ) {
        output_report("[INFO]\tUser requested reindexing of database directory files - slow!\n");
    }
    my $run_directory = "$Working_Dir\/$User_Runid";

    # Now we can make the directories
    setup_main_directories( $run_directory, $options{f} );

    # Report
    output_report("[INFO]\tRun ID: $User_Runid\n[INFO]\tDirectory: $run_directory\n");

    # search options
    $Search_Program    = $paramaters->{search}->{program};
    $Search_Subprogram = $paramaters->{search}->{subprogram};
    $Search_Evalue     = $paramaters->{search}->{evalue};
    $Search_Tophits    = $paramaters->{search}->{top_hits};
    $Search_Maxlength  = $paramaters->{search}->{max_length};
    @Search_Special_Taxa    = split /,/, $paramaters->{search}->{special_taxa};
    $Search_Special_Tophits = $paramaters->{search}->{special_top_hits};
    $Search_Threads         = $paramaters->{search}->{threads};

    #$Search_Other = $paramaters->{search}->{other};                                # no default

    # alignment options
    $Alignment_Program = $paramaters->{alignment}->{program} || 'mafft';
    $Alignment_Options = $paramaters->{alignment}->{options};              # no default
    $Alignment_Threads = $paramaters->{alignment}->{threads} || '1';

    # masking options
    $Masking_Cutoff1 = $paramaters->{masking}->{cutoff_1} || '50';
    $Masking_Cutoff2 = $paramaters->{masking}->{cutoff_2} || '20';

    # tree building options
    $Tree_Program = $paramaters->{trees}->{program} || 'FastTreeMP';
    $Tree_Options = $paramaters->{trees}->{options};                       # no default
    $Tree_Mintaxa = $paramaters->{trees}->{min_taxa} || '4';
    $Tree_Model   = $paramaters->{trees}->{model} || 'PROTGAMMAAUTO';
    $Tree_Bs      = $paramaters->{trees}->{bootstraps} || '100';

    # directory options
    $Seq_Data = $paramaters->{directories}->{database};                    # no default

    # filtering options
    $Filter_File = $paramaters->{filter}->{filename};                      # no default

    # only run search (blast) step
    if ( $options{b} ) {

        # Run the user selected sequence search option
        my $start_time = timing('start');

        if ( -d $input_seqs_fname ) {
            print "Running: Search ($Search_Program) from Ortho Groups Directory ONLY\n";

            search_step_ortho_groups( \@taxa_array, $input_seqs_fname );
        }
        else {
            print "Running: Search ($Search_Program) ONLY\n";
            search_step( \@taxa_array, $input_seqs_fname );
        }

        my $end_time = timing( 'end', $start_time );
    }

    # Filtering step. This options removes trees that do not have hits to any
    # members of a group of taxa or higher level classification
    # e.g. if the potential tree does not have "plants" exclude it
    if ( $options{f} ) {
        print "Running: Filtering un-aligned seqs that don't contain taxa/groups from $Filter_File\n";
        my $start_time = timing('start');
        filtering_step();
        my $end_time = timing( 'end', $start_time );
    }

    # only run alignment step
    if ( $options{a} ) {
        print "Running: Alignment ($Alignment_Program) ONLY\n";
        my $start_time = timing('start');
        alignment_step();
        my $end_time = timing( 'end', $start_time );
    }

    # only run mask step
    if ( $options{m} ) {
        print "Running: Masking with trimal ONLY\n";
        my $start_time = timing('start');
        masking_step();
        my $end_time = timing( 'end', $start_time );
    }

    # only run tree building step (o is for orchard)
    if ( $options{o} ) {
        print "Running: Tree Reconstruction with $Tree_Program ONLY\n";
        my $start_time = timing('start');
        tree_step();
        my $end_time = timing( 'end', $start_time );
    }

    # run default
    if ( !$options{b} && !$options{a} && !$options{m} && !$options{o} && !$options{f} ) {

        # run all steps, but all blasts first then amt steps
        print "Running: ALL Steps; all searches ($Search_Program) first!\n";

        # Run the user selected sequence search option
        my $start_time = timing('start');

        if ( -d $input_seqs_fname ) {
            print "Running: Search ($Search_Program) from Ortho Groups Directory ONLY\n";

            search_step_ortho_groups( \@taxa_array, $input_seqs_fname );
        }
        else {
            print "Running: Search ($Search_Program) ONLY\n";
            search_step( \@taxa_array, $input_seqs_fname );
        }

        my $end_time = timing( 'end', $start_time );

        # Run the user selected alignment option
        print "Running: Alignment with $Alignment_Program\n";
        $start_time = timing('start');
        alignment_step();
        $end_time = timing( 'end', $start_time );

        # Run the user selected masking option
        print "Running: Masking with trimal\n";
        $start_time = timing('start');
        masking_step();
        $end_time = timing( 'end', $start_time );

        # Run the user selected tree building options
        print "Running: Tree Reconstruction with $Tree_Program\n";
        $start_time = timing('start');
        tree_step();
        $end_time = timing( 'end', $start_time );
    }
}
else {
    display_help();
}

###########################################################
##           Main Subroutines 	                         ##
###########################################################

#################################################
##           Filtering Subroutine              ##
#################################################

sub filtering_step {

    # directories
    my $sequence_directory = "$Working_Dir\/$User_Runid\/seqs";
    my $excluded_directory = "$Working_Dir\/$User_Runid\/excluded";

    # get list of sequences that were not excluded in previous stage
    my @seq_file_names = glob "$sequence_directory\/*.fas";

    # read in the list of taxa that must be in set of seqs
    open my $Filter_File, '<', $Filter_File;
    my @filter_list = <$Filter_File>;
    close $Filter_File;

    # Let's not continue if there are no files!
    if ( $#seq_file_names == -1 ) {
        print "There are no files in the seqs directory, did you run blast?\n";
        exit;
    }

    # iterate through the stream of sequences and perform each search
    output_report("[INFO]\tFILTERING: $#seq_file_names sequence files using $Filter_File\n");

    for my $i ( 0 .. $#seq_file_names ) {

        my $current_sequences = $seq_file_names[$i];
        my ( $file, $dir, $ext ) = fileparse $current_sequences, '.fas';

        print "Filtering $current_sequences with $Filter_File\n";

        my @taxa_accessions = get_accession_list("$sequence_directory\/$file$ext");
        my $missing_count   = 0;
        my $matched         = 0;
        print "MC= $missing_count\t TA= $#taxa_accessions\n";
        foreach my $accession (@taxa_accessions) {

            my $taxon_name = get_taxon_names($accession);

            #print "$taxon_name\t";

            # if the accession is the same as the taxon name
            # it's either a seed or not in the DB, ignore it
            if ( $accession eq $taxon_name ) {

                #print "Probaly the 'seed' sequence, ignoring\n";
                $missing_count++;
                next;
            }

            # modify the taxon name to be Genus species only
            # this is because the ncbi taxonomy is EXACT match
            # only, so if your taxa has a strain like:
            # 'MAD 698-R' but in NCBI it is 'Mad-698-R' you
            # won't get a match, but it will cause this script
            # to false positive a match...
            # this goes for spelling mistakes too....!
            my ( $genus, $species ) = split ' ', $taxon_name;
            $taxon_name = "$genus $species";

            # get the taxonomic lineage and put it in to an array
            # then split it up with pipes, and it works in a search
            # not sure why compared to just the returned comma
            # separated string...hey ho!
            my @taxonomy = split ',', get_taxonomy($taxon_name);
            my $taxonomy_search = join( '|', @taxonomy );

            #print "$taxonomy_search\n";

            # check for species name in the filter list
            if ( my ($matched_taxa) = grep { $_ eq $taxon_name } @filter_list ) {

                print "\tMATCH1: $matched_taxa\tTaxa: $taxon_name\n";
                $matched = 1;
                print "MC1 $missing_count\tM: $matched\n";
                output_report("[INFO]\tFILTERING1: Matched $accession - $taxon_name to $matched_taxa, retained.\n");
                last;
            }

            # then check for taxonomic lineages
            elsif ( my ($matched_taxonomy) = grep { /$taxonomy_search/ } @filter_list ) {

                print "\t\tMATCH2:\tTaxonomy: @taxonomy\n";
                print "MC2 $missing_count\tM: $matched\n";
                $matched = 1;
                output_report(
                    "[INFO]\tFILTERING2: Matched $accession - $taxon_name within $matched_taxonomy, retained.\n");
                last;
            }
            else {
                $missing_count++;
                print "MC: $missing_count\tM: $matched\n";
                print "No matched filtered taxa/taxonomy\n";
            }
        }

        # at the end of all the searching if our count
        # equals the number of accession or is larger
        # then we likely haven't hit a taxon from the filter list and doubly
        # so if out match status is still 0
        if ( $missing_count >= $#taxa_accessions || $matched eq 0 ) {
            my $command = "mv $current_sequences $excluded_directory";
            system($command);
            output_report("[INFO]\tFILTERING3: No matches for $current_sequences, excluding\n");
            print "MC3 $missing_count\tM: $matched\n";
            print "Moving file to excluded...\n";
        }

        #elsif ( $matched eq 0 ) {
        #    my $command = "mv $current_sequences $excluded_directory";
        #    system($command);
        #    output_report("[INFO]\tFILTERING4: No matches for $current_sequences, excluding\n");
        #    print "MC4 $missing_count\tM: $matched\n";
        #    print "Moving file to excluded...\n";
        #}
        #}
    }
}

sub get_taxonomy {

    my $taxon_name = shift;

    ## Eventually I will deprecate the web interface
    ## as it is slow and go to using local caches files
    ## but not today!

    ### Bio::Taxon
    ## Local Files from ftp://ftp.ncbi.nih.gov/pub/taxonomy/ - taxcat
    ##my $dbh = Bio::DB::Taxonomy->new(
    ##    -source => 'flatfile',
    ##    #    -directory => '/home/cs02gl/Desktop/genomes/taxonomy',
    ##    -nodesfile => '/home/cs02gl/Desktop/genomes/taxonomy/nodes.dmp',
    ##    -namesfile => '/home/cs02gl/Desktop/genomes/taxonomy/names.dmp'
    ##);

    ## Entrez providing stable connection.
    my $dbh = Bio::DB::Taxonomy->new( -source => 'entrez' );

    # Retreive taxon_name
    my $unknown = $dbh->get_taxon( -name => "$taxon_name" );

    # build an empty tree
    my $tree_functions = Bio::Tree::Tree->new();

    # and get the lineage of the taxon_name
    my @lineage = $tree_functions->get_lineage_nodes($unknown);

    # Then we can extract the name of each node, which will give us the Taxonomy lineages...
    my $taxonomy;
    foreach my $item (@lineage) {
        my $name = $item->node_name;

        #my $rank = $item->rank;
        $taxonomy = "$taxonomy$name,";
    }

    return $taxonomy;
}

# this gets all accessions from a a file
# returned in an array
sub get_accession_list {

    my $filename   = shift;
    my @taxa_array = ();

    my $input_stream = Bio::SeqIO->new( '-file' => $filename );
    while ( my $seq = $input_stream->next_seq() ) {
        my $seq_id = $seq->id;
        push( @taxa_array, $seq_id );
    }
    return @taxa_array;
}

# return the taxa name from the given accession
# this currently runs with the mysql database
# in future this will become depracated
# but I can't give it up just yet
sub get_taxon_names {

    my $seqID = shift;

    #my $taxon_name = "";

    # Database Settings
    my $USER      = "orchard";
    my $PASSWORD  = "jcsy4s8b";
    my $SERVER    = "144.173.27.211";
    my $TABLENAME = "cider";
    my $DATABASE  = "dbi:mysql:new_proteins:$SERVER";    # do not edit

    # Database handle
    my $dbh = DBI->connect( $DATABASE, $USER, $PASSWORD ) or die "\nError ($DBI::err):$DBI::errstr\n";

    my $statement = $dbh->prepare("SELECT species FROM $TABLENAME where protein_ID='$seqID'")
      or die "\nError ($DBI::err):$DBI::errstr\n";
    $statement->execute or die "\nError ($DBI::err):$DBI::errstr\n";
    my ($taxon_name) = $statement->fetchrow_array;

    if ( !defined $taxon_name ) {
        $taxon_name = "$seqID";
    }

    return "$taxon_name";
}

#################################################
##           Tree Subroutines                  ##
#################################################

sub tree_step {

    my $masks_directory = "$Working_Dir\/$User_Runid\/masks";
    my $tree_directory  = "$Working_Dir\/$User_Runid\/trees";

    # get list of sequences that were not excluded in previous stage
    my @mask_file_names = glob "$masks_directory\/*.afa-tr";

    output_report("[INFO]\tTREE BUILDING: $#mask_file_names alignments using $Tree_Program\n");

    if ( $Tree_Program =~ /fasttree/ism ) {
        for my $i ( 0 .. $#mask_file_names ) {
            my ( $file, $dir, $ext ) = fileparse $mask_file_names[$i], '\.afa-tr';

            run_fasttree("$masks_directory\/$file\.afa-tr");
        }
    }
    elsif ( $Tree_Program =~ /iqtree/ism ) {
        for my $i ( 0 .. $#mask_file_names ) {
            my ( $file, $dir, $ext ) = fileparse $mask_file_names[$i], '\.afa-tr';

            run_iqtree("$masks_directory\/$file\.afa-tr");
        }
    }
    else {
        for my $i ( 0 .. $#mask_file_names ) {
            my ( $file, $dir, $ext ) = fileparse $mask_file_names[$i], '\.afa-tr';

            run_raxml_ml_rapid_bd("$masks_directory\/$file\.afa-tr");
        }
    }
}

sub run_iqtree {
    my $masked_sequences = shift;
    my $iqtree_command;

    my ( $file, $dir, $ext ) = fileparse $masked_sequences, '\.afa\-tr';

    $iqtree_command = "iqtree -s $masked_sequences $Tree_Options";

}

sub run_raxml_ml_rapid_bd {
    my $masked_sequences = shift;
    my $raxmltree_command;

    my ( $file, $dir, $ext ) = fileparse $masked_sequences, '\.afa\-tr';

    if ( $Tree_Program =~ /raxmlHPC-PTHREADS/ism ) {
        $raxmltree_command = "raxmlHPC-PTHREADS -T $Search_Threads";
        $raxmltree_command = " -f a";                                                           # fast ML
        $raxmltree_command = " -m $Tree_Model";                                                 # model
        $raxmltree_command = " -p " . int( rand(10000) ) . " -x " . int( rand(10000) ) . "";    # rapid BS mode
        $raxmltree_command = " -# $Tree_Bs";                                                    # e.g. 100 BS
        $raxmltree_command = " -s $masked_sequences";                                           # masked alignment
    }
    else {
        $raxmltree_command = "raxmlHPC";
        $raxmltree_command = " -f a";                                                           # fast ML
        $raxmltree_command = " -m $Tree_Model";                                                 # model
        $raxmltree_command = " -p " . int( rand(10000) ) . " -x " . int( rand(10000) ) . "";    # rapid BS mode
        $raxmltree_command = " -# $Tree_Bs";                                                    # e.g. 100 BS
        $raxmltree_command = " -s $masked_sequences";                                           # masked alignment
    }

    system($raxmltree_command);
}

sub run_fasttree {
    my $masked_sequences = shift;
    my $fasttree_command;

    my ( $file, $dir, $ext ) = fileparse $masked_sequences, '\.afa\-tr';

    print "\tCalculating FastTree on $file\n";

    if ( $Tree_Program =~ /fasttreemp/ism ) {
        $fasttree_command = "FastTreeMP $Tree_Options";
        $fasttree_command .= " $masked_sequences";
        $fasttree_command .= " > $Working_Dir\/$User_Runid\/trees\/$file\_FT.tree";
    }
    else {
        $fasttree_command = "FastTree $Tree_Options";
        $fasttree_command .= " $masked_sequences";
        $fasttree_command .= " > $Working_Dir\/$User_Runid\/trees\/$file\_FT.tree";
    }

    system($fasttree_command);
}

#################################################
##           Masking Subroutines               ##
#################################################

sub masking_step {

    # directory variables
    my $alignments_directory = "$Working_Dir\/$User_Runid\/alignments";
    my $masks_directory      = "$Working_Dir\/$User_Runid\/masks";

    # get list of sequences that were not excluded in previous stage
    my @algn_file_names = glob "$alignments_directory\/*.afa";

    # iterate through the stream of sequences and perform each search
    output_report("[INFO]\tMASKING: $#algn_file_names alignments using trimal\n");

    for my $i ( 0 .. $#algn_file_names ) {
        my $current_sequences = $algn_file_names[$i];

        my ( $file, $dir, $ext ) = fileparse $current_sequences, '\.afa';

        # run trimal and report mask length with -nogaps option
        my $mask_length = run_trimal( $current_sequences, '-nogaps' );

        # if the length is less than the first limit
        if ( $mask_length <= $Masking_Cutoff1 ) {

            # then re-run trimal and report mask length with -automated1 option
            $mask_length = run_trimal( $current_sequences, '-automated1' );

            print "Mask Length is ($mask_length) and is";

            # if the length is less than the second limit (always the smaller)
            if ( $mask_length <= $Masking_Cutoff2 ) {
                print " not OK. Excluding sequence.\n";

                # then abandon this sequence, report it, move to excluded
                output_report(
"[WARN]\t$file does not satisfy trimal cutoffs ($Masking_Cutoff1 or $Masking_Cutoff2). Moved to excluded directory.\n"
                );
                system
"mv $Working_Dir\/$User_Runid\/masks\/$file\.afa-tr $Working_Dir\/$User_Runid\/excluded\/$file\.afa-tr";
            }
            else {
                print " OK.\n";
            }
        }
        else {
            print "\t\tMask Length: $mask_length is OK\n";
        }
    }
}

sub run_trimal {

    my $aligned_sequences = shift;
    my $option            = shift;
    my $masks_directory   = "$Working_Dir\/$User_Runid\/masks";

    my ( $file, $dir, $ext ) = fileparse $aligned_sequences, '\.afa';

    print "\tRunning trimal $option on $aligned_sequences\n";

    my $trimal_command = 'trimal';
    $trimal_command .= " -in $aligned_sequences";
    $trimal_command .= " -out $masks_directory\/$file\.afa\-tr";
    $trimal_command .= " $option";

    system($trimal_command);

    my $length = mask_check("$masks_directory\/$file\.afa\-tr");

    return $length;
}

sub mask_check {

    my $aligned_sequences = shift;

    # I want to silence the "MSG: Got a sequence without letters. Could not guess alphabet"
    # when trimal/gblock produce an alignment with 0 letters...as I deal with it above
    # anybody have any ideas as verbose -1 isn't working!
    my $phylipstream = Bio::SeqIO->new( '-file' => $aligned_sequences, -verbose => -1 );
    my $seqobj       = $phylipstream->next_seq;
    my $length       = $seqobj->length;

    return $length;
}

#################################################
##           Alignment Subroutines             ##
#################################################

sub alignment_step {

    ## need to add a stage to copy files from orthogroup dir to seqs
    ## if we are going to be bypassing the search step...

    # directory variables
    my $sequence_directory   = "$Working_Dir\/$User_Runid\/seqs";
    my $alignments_directory = "$Working_Dir\/$User_Runid\/alignments";

    # get list of sequences that were not excluded in previous stage
    my @seq_file_names = glob "$sequence_directory\/*.fas";

    # iterate through the stream of sequences and perform each search
    output_report("[INFO]\tALIGNING: $#seq_file_names alignments using $Alignment_Program\n");

    for my $i ( 0 .. $#seq_file_names ) {

        my $current_sequences = $seq_file_names[$i];

        print "Aligning $current_sequences with $Alignment_Program\n";

        my ( $file, $dir, $ext ) = fileparse $current_sequences, '\.fas';

        given ($Alignment_Program) {
            when (/mafft/ism) {
                my $mafft_command = "mafft $Alignment_Options";
                $mafft_command .= " --thread $Alignment_Threads";
                $mafft_command .= " $current_sequences > $alignments_directory\/$file\.afa";

                system($mafft_command);
            }
            when (/muscle/ism) {

                my $muscle_command = "muscle $Alignment_Options";
                $muscle_command .= " -in $current_sequences -out $alignments_directory\/$file\.afa";

                system($muscle_command);
            }
        }
    }
}

#################################################
##           Search Subroutines                ##
#################################################
sub search_step_ortho_groups {

    # get all of the values passed to the sub...
    my ( $taxa_array_ref, $input_seqs_fname ) = @_;

    # files from orthogroup/multi-seq directory
    my @file_names        = glob "$input_seqs_fname\/*.fasta $input_seqs_fname\/*.fas";
    my $ortho_files_total = @file_names;

    # dereference the array
    my @taxa_array = @{$taxa_array_ref};

    # we should always have one sequence, right!?
    my $ortho_files_count = 1;

    foreach my $orthofile (@file_names) {

        # we should always have one sequence, right!?
        my $input_seqs_count = 1;

        print "Starting: $ortho_files_count of $ortho_files_total = $orthofile\n";
        output_report("[INFO]\tStarting $ortho_files_count of $ortho_files_total = $orthofile\n");

        my ( $ortho_file, $ortho_dir, $ortho_ext ) = fileparse( $orthofile, qr'\..*' );

        #
        my $num_hit_seqs = 0;
        my $sequence_name;

        # open bioperl seqio object with user input sequences
        my $seq_in = Bio::SeqIO->new( -file => "<$orthofile" );

        # I am still relying on 'grep' to count the number of sequences
        # there is no way to get this directly from the Bio::Seq object
        # without needless iteration. Anyone?
        chomp( my $input_seqs_total = `grep -c ">" $orthofile` );

        # iterate through the stream of sequences and perform each search
        output_report("[INFO]\tStarting $input_seqs_total searches using $Search_Program($Search_Subprogram)\n");

        while ( my $seq = $seq_in->next_seq() ) {

            print "\tProcessing: $input_seqs_count of $input_seqs_total\n";
            output_report("[INFO]\tProcessing: $input_seqs_count of $input_seqs_total\n");

            # Get Sequence Name
            $sequence_name = $seq->id;

            # orthofroups have a pipe symbol that throws off commands
            # replace it with an underscore
            $sequence_name =~ s/\|/\_/;

            # Output Query Sequence
            my $query_out = Bio::SeqIO->new(
                -file   => ">$sequence_name\_query.fas",
                -format => 'fasta'
            );
            $query_out->write_seq($seq);

            # Running total of hits, note '>>' append
            my $hits_out = Bio::SeqIO->new(
                -file   => ">>$sequence_name\_hits.fas",
                -format => 'fasta'
            );
            $hits_out->write_seq($seq);

            # This has been added to catch the results of all sequence blasts
            # for each of the orthogroups, look for the sort/uniq step below too...
            my $ortho_seqs_out = Bio::SeqIO->new(
                -file   => ">>$ortho_file\_ortho\_hits.fas",
                -format => 'fasta'
            );
            $ortho_seqs_out->write_seq($seq);

            for ($Search_Program) {
                when (/BLAST[+]/ism) {

                    print "\tRunning: blast+ on $sequence_name\n";
                    run_blast_plus( \@taxa_array, "$sequence_name\_query.fas", $sequence_name );
                }
                when (/^BLAST$/ism) {

                    print "\tRunning: legacy blast\n";
                    run_blast_legacy( \@taxa_array, "$sequence_name\_query.fas", $sequence_name );
                }
                when (/BLAT/ism) {

                    print "\tRunning: blat\n";
                    run_blat( \@taxa_array, "$sequence_name\_query.fas", $sequence_name );
                }

                #when (/VSEARCH/ism) {
                #    print "\tRunning: vsearch\n";
                #    run_vsearch( \@taxa_array, "$sequence_name\_query.fas", $sequence_name );
                #}
                when (/USEARCH/ism) {

                    print "\tRunning: usearch\n";
                    run_usearch( \@taxa_array, "$sequence_name\_query.fas", $sequence_name );
                }
                default {

                    print "\tRunning: (default) blast+ on $sequence_name\n";

                    # lets just make sure we are using the defaults and not user passed variables
                    $Search_Program    = 'blast+';
                    $Search_Subprogram = 'blastp';
                    run_blast_plus( \@taxa_array, "$sequence_name\_query.fas", $sequence_name );
                }
            }
            $input_seqs_count++;

            # remove the query hits after each run as they are now stored
            # in the *ortho_seqs.fas file
            unlink "$Working_Dir\/$sequence_name\_query.fas";

            # append the current seed sequence's hits to the ortho_seqs.fas file
            system "cat $Working_Dir\/$sequence_name\_hits.fas >> $Working_Dir\/$ortho_file\_ortho\_hits.fas";

            # then get rid of the query file as we don't need it anymore
            unlink "$Working_Dir\/$sequence_name\_hits.fas";
        }
        $ortho_files_count++;

        # now we can remove duplicate sequences based on accession only
        remove_duplicate_sequences("$Working_Dir\/$ortho_file\_ortho\_hits.fas");

        # get the number of sequences in the file
        chomp( $num_hit_seqs = `grep -c ">" $Working_Dir\/$ortho_file\_ortho\_hits.fas` );

        if ( $num_hit_seqs <= $Tree_Mintaxa ) {
            output_report("[WARN]\t$ortho_file: Too few hits\n");
            unlink "$Working_Dir\/$ortho_file\_ortho\_hits.fas";
            system
              "mv $Working_Dir\/$ortho_file\_ortho\_non\_redundant\_hits.fas $Working_Dir\/$User_Runid\/excluded\/";
        }
        else {
            unlink "$Working_Dir\/$ortho_file\_ortho\_hits.fas";
            system
"mv $Working_Dir\/$ortho_file\_ortho\_non\_redundant\_hits.fas $Working_Dir\/$User_Runid\/seqs\/$ortho_file\_hits.fas";

            # remove selenocystein and non-alpha numeric characters as they cause BLAST/MAFFT
            # to complain/crash (and the --anysymbol option produces terrible alignments)
            # I have to check for this as some genome projects are just full of junk!
            system "sed -i \'/^>/!s/U/X/g\' $Working_Dir\/$User_Runid\/seqs\/$ortho_file\_hits.fas";
        }
    }
}

sub remove_duplicate_sequences {

    # based on IDs, this keeps our ortholog/multi-seq original seqs
    # in the file, as duplicates, but shouldn't affect align/tree
    my $file_name = shift;
    my ( $FILE, $DIR, $EXT ) = fileparse( $file_name, "\_hits.fas" );

    my %digests;
    my $in = Bio::SeqIO->new(
        - file   => "$file_name",
        - format => 'fasta'
    );

    my $out = Bio::SeqIO->new(
        - file   => ">$Working_Dir\/$FILE\_non_redundant$EXT",
        - format => 'fasta'
    );

    while ( my $seq = $in->next_seq ) {
        $digests{ md5_hex( $seq->id ) }++;
        if ( $digests{ md5_hex( $seq->id ) } == 1 ) {
            $out->write_seq($seq);
        }
    }
}

sub search_step {

    # get all of the values passed to the sub...
    my ( $taxa_array_ref, $input_seqs_fname ) = @_;

    # dereference the array
    my @taxa_array = @{$taxa_array_ref};

    # we should always have one sequence, right!?
    my $input_seqs_count = 1;

    #
    my $num_hit_seqs = 0;
    my $sequence_name;

    # open bioperl seqio object with user input sequences
    my $seq_in = Bio::SeqIO->new( -file => "<$input_seqs_fname" );

    # I am still relying on 'grep' to count the number of sequences
    # there is no way to get this directly from the Bio::Seq object
    # without needless iteration. Anyone?
    chomp( my $input_seqs_total = `grep -c ">" $input_seqs_fname` );

    # iterate through the stream of sequences and perform each search
    output_report("[INFO]\tStarting $input_seqs_total searches using $Search_Program($Search_Subprogram)\n");
    while ( my $seq = $seq_in->next_seq() ) {

        print "Processing: $input_seqs_count of $input_seqs_total\n";

        # Get Sequence Name
        $sequence_name = $seq->id;

        # Output Query Sequence
        my $query_out = Bio::SeqIO->new(
            -file   => ">$sequence_name\_query.fas",
            -format => 'fasta'
        );
        $query_out->write_seq($seq);

        # Running total of hits, note '>>' append
        my $hits_out = Bio::SeqIO->new(
            -file   => ">>$sequence_name\_hits.fas",
            -format => 'fasta'
        );
        $hits_out->write_seq($seq);

        for ($Search_Program) {
            when (/BLAST[+]/ism) {

                print "\tRunning: blast+ on $sequence_name\n";
                $num_hit_seqs = run_blast_plus( \@taxa_array, $input_seqs_fname, $sequence_name );
            }
            when (/^BLAST$/ism) {

                print "\tRunning: legacy blast\n";
                $num_hit_seqs = run_blast_legacy( \@taxa_array, $input_seqs_fname, $sequence_name );
            }
            when (/BLAT/ism) {

                print "\tRunning: blat\n";
                $num_hit_seqs = run_blat( \@taxa_array, $input_seqs_fname, $sequence_name );
            }

            #when (/VSEARCH/ism) {
            #    print "\tRunning: vsearch\n";
            #    run_vsearch( \@taxa_array, "$sequence_name\_query.fas", $sequence_name );
            #}
            when (/USEARCH/ism) {

                print "\tRunning: usearch\n";
                $num_hit_seqs = run_usearch( \@taxa_array, $input_seqs_fname, $sequence_name );
            }
            default {
                # lets just make sure we are using the defaults and not user passed variables
                $Search_Program    = 'blast+';
                $Search_Subprogram = 'blastp';
                print "\tRunning: (default) blast+ on $sequence_name\n";
                $num_hit_seqs = run_blast_plus( \@taxa_array, $input_seqs_fname, $sequence_name );
            }
        }
        $input_seqs_count++;

        if ( $num_hit_seqs <= $Tree_Mintaxa ) {
            output_report("[WARN]\t$sequence_name: Too few hits\n");
            unlink "$Working_Dir\/$sequence_name\_query.fas";
            system "mv $sequence_name\_hits.fas $Working_Dir\/$User_Runid\/excluded\/";
        }
        else {
            unlink "$Working_Dir\/$sequence_name\_query.fas";
            system "mv $sequence_name\_hits.fas $Working_Dir\/$User_Runid\/seqs\/";

            # remove selenocystein and non-alpha numeric characters as they cause BLAST/MAFFT
            # to complain/crash (and the --anysymbol option produces terrible alignments)
            # I have to check for this as some genome projects are just full of junk!
            system "sed -i \'/^>/!s/U|\\w/X/g\' $Working_Dir\/$User_Runid\/seqs\/$sequence_name\_hits.fas";
        }
    }
    return;
}

sub run_blast_plus {

    my ( $taxa_array_ref, $input_seqs_fname, $sequence_name ) = @_;

    my @taxa_array = @{$taxa_array_ref};
    my $taxa_total = @taxa_array;
    my $taxa_count = 1;

    my $sequence_name_for_blast = $sequence_name;
    $sequence_name_for_blast =~ s/\s+/\_/gms;    # Replace spaces with '_'

    while (@taxa_array) {

        # Current Taxa Name
        my $taxa_name = shift(@taxa_array);

        my $taxa_name_for_blast = $taxa_name;
        $taxa_name_for_blast =~ s/\s+/\_/gms;    # Replace spaces with '_'

        if ( $taxa_name =~ m/^#/sm ) {
            print "\t\tSkipping commented out $taxa_name\n";
            output_report("[INFO]\t$sequence_name: Skipping commented out $taxa_name\n");
            $taxa_count++;
        }
        else {

            # Blast Output Filename
            my $search_output = "$sequence_name\_v\_$taxa_name_for_blast\.$Search_Subprogram";

            # BLAST Output Progress
            printf
              "\t\t>: $Search_Subprogram: $taxa_count of $taxa_total\n\e[A";    # - $taxa_name\n\e[A";    # Progress...

            my $database = $taxa_name_for_blast . '.fas';

            # In the future I might consider replacing grep with List:MoreUtils 'any' to save
            # large searches
            if ( grep { $_ eq $taxa_name } @Search_Special_Taxa ) {

                # blast(x) from blast+ package command
                # we will use tabulated output as it's smaller than XML
                # and we don't really need much information other than the hit ID
                my $blast_command = "$Search_Subprogram -task $Search_Subprogram";
                $blast_command .= " -db $Seq_Data\/$database";
                $blast_command .= " -query $sequence_name_for_blast\_query.fas";
                $blast_command .= " -out $search_output";
                $blast_command .= " -evalue $Search_Evalue";
                $blast_command .= ' -outfmt 6';
                $blast_command .= " -max_target_seqs $Search_Special_Tophits";
                $blast_command .= " -num_threads $Search_Threads";

                #$blast_command .= " $Search_Other";

                system($blast_command);
                parse_search_output( \@taxa_array, $input_seqs_fname, $sequence_name, $taxa_name, $database );

                output_report(
                    "[INFO]\t$sequence_name: Using special -max_target_seqs $Search_Special_Tophits for $taxa_name\n");
            }
            else {
                # blast(x) from blast+ package command
                # we will use tabulated output as it's smaller than XML
                # and we don't really need much information other than the hit ID
                my $blast_command = "$Search_Subprogram -task $Search_Subprogram";
                $blast_command .= " -db $Seq_Data\/$database";
                $blast_command .= " -query $sequence_name_for_blast\_query.fas";
                $blast_command .= " -out $search_output";
                $blast_command .= " -evalue $Search_Evalue";
                $blast_command .= ' -outfmt 6';
                $blast_command .= " -max_target_seqs $Search_Tophits";
                $blast_command .= " -num_threads $Search_Threads";

                #$blast_command .= " $Search_Other";

                system($blast_command);
                parse_search_output( \@taxa_array, $input_seqs_fname, $sequence_name, $taxa_name, $database );
            }

            $taxa_count++;
        }
    }

    # Find total number of hits
    # I am still relying on 'grep' to count the number of sequences
    # there is no way to get this directly from the Bio::Seq object
    # without needless iteration. Anyone?
    chomp( my $hit_seqs_total = `grep -c ">" $sequence_name\_hits.fas` );

    #
    print "\n\tNumber of Sequences found = $hit_seqs_total\n";
    return $hit_seqs_total;
}

sub run_blast_legacy {

    my ( $taxa_array_ref, $input_seqs_fname, $sequence_name ) = @_;

    my @taxa_array = @{$taxa_array_ref};
    my $taxa_total = @taxa_array;
    my $taxa_count = 1;

    my $sequence_name_for_blast = $sequence_name;
    $sequence_name_for_blast =~ s/\s+/\_/gms;    # Replace spaces with '_'

    while (@taxa_array) {

        # Current Taxa Name
        my $taxa_name = shift(@taxa_array);

        my $taxa_name_for_blast = $taxa_name;
        $taxa_name_for_blast =~ s/\s+/\_/gms;    # Replace spaces with '_'

        if ( $taxa_name =~ m/^#/sm ) {
            print "\t\tSkipping commented out $taxa_name\n";
            output_report("[INFO]\t$sequence_name: Skipping commented out $taxa_name\n");
            $taxa_count++;
        }
        else {

            # Blast Output Filename
            my $search_output = "$sequence_name\_v\_$taxa_name_for_blast\.$Search_Subprogram";

            # BLAST Output Progress
            printf "\t\t>: $Search_Subprogram: $taxa_count of $taxa_total\n\e[A";    # Progress...

            my $database = $taxa_name_for_blast . '.fas';

            # In the future I might consider replacing grep with List:MoreUtils 'any' to save
            # large searches
            if ( grep { $_ eq $taxa_name } @Search_Special_Taxa ) {

                # blast(x) from legacy blast package command
                # we will use tabulated output as it's smaller than XML
                # and we don't really need much information other than the hit ID
                my $blast_command = "blastall -p $Search_Subprogram";
                $blast_command .= " -d $Seq_Data\/$database";
                $blast_command .= " -i $sequence_name_for_blast\_query.fas";
                $blast_command .= " -o $search_output";
                $blast_command .= " -e $Search_Evalue";
                $blast_command .= ' -m 8';
                $blast_command .= " -b $Search_Special_Tophits";

                #$blast_command .= " -v $Search_Tophits";
                $blast_command .= " -a $Search_Threads";

                #$blast_command .= " $Search_Other";

                system($blast_command);
                parse_search_output( \@taxa_array, $input_seqs_fname, $sequence_name, $taxa_name, $database );

                output_report(
                    "[INFO]\t$sequence_name: Using special -max_target_seqs $Search_Special_Tophits for $taxa_name\n");

            }
            else {
                # blast(x) from legacy blast package command
                # we will use tabulated output as it's smaller than XML
                # and we don't really need much information other than the hit ID
                my $blast_command = "blastall -p $Search_Subprogram";
                $blast_command .= " -d $Seq_Data\/$database";
                $blast_command .= " -i $sequence_name_for_blast\_query.fas";
                $blast_command .= " -o $search_output";
                $blast_command .= " -e $Search_Evalue";
                $blast_command .= ' -m 8';
                $blast_command .= " -b $Search_Tophits";

                #$blast_command .= " -v $Search_Tophits";
                $blast_command .= " -a $Search_Threads";

                #$blast_command .= " $Search_Other";

                system($blast_command);
                parse_search_output( \@taxa_array, $input_seqs_fname, $sequence_name, $taxa_name, $database );
            }

            $taxa_count++;
        }
    }

    # Find total number of hits
    # I am still relying on 'grep' to count the number of sequences
    # there is no way to get this directly from the Bio::Seq object
    # without needless iteration. Anyone?
    chomp( my $hit_seqs_total = `grep -c ">" $sequence_name\_hits.fas` );

    #
    print "\n\tNumber of Sequences found = $hit_seqs_total\n";
    return $hit_seqs_total;
}

sub run_blat {

    my ( $taxa_array_ref, $input_seqs_fname, $sequence_name ) = @_;

    my @taxa_array = @{$taxa_array_ref};
    my $taxa_total = @taxa_array;
    my $taxa_count = 1;

    my $sequence_name_for_blast = $sequence_name;
    $sequence_name_for_blast =~ s/\s+/\_/gms;    # Replace spaces with '_'

    while (@taxa_array) {

        # Current Taxa Name
        my $taxa_name = shift(@taxa_array);

        my $taxa_name_for_blast = $taxa_name;
        $taxa_name_for_blast =~ s/\s+/\_/gms;    # Replace spaces with '_'

        if ( $taxa_name =~ m/^#/sm ) {
            print "\t\tSkipping commented out $taxa_name\n";
            output_report("[INFO]\t$sequence_name: Skipping commented out $taxa_name\n");
            $taxa_count++;
        }
        else {

            # Blast Output Filename
            my $search_output = "$sequence_name\_v\_$taxa_name_for_blast\.$Search_Subprogram";

            # BLAST Output Progress
            printf "\t\t>: $Search_Subprogram: $taxa_count of $taxa_total\n\e[A";    # Progress...

            my $database = $taxa_name_for_blast . '.fas';

            # blat from blat package command
            # we will use tabulated output as it's smaller than XML
            # and we don't really need much information other than the hit ID
            my $blat_command = 'blat';
            $blat_command .= ' -prot';                                  # this should be a user option eventually
            $blat_command .= " $Seq_Data\/$database";
            $blat_command .= " $sequence_name_for_blast\_query.fas";    # query file
            $blat_command .= ' -out=blast8';
            $blat_command .= " $search_output";                         # output filename

            system($blat_command);
            parse_search_output( \@taxa_array, $input_seqs_fname, $sequence_name, $taxa_name, $database );

            $taxa_count++;
        }
    }

    # Find total number of hits
    # I am still relying on 'grep' to count the number of sequences
    # there is no way to get this directly from the Bio::Seq object
    # without needless iteration. Anyone?
    chomp( my $hit_seqs_total = `grep -c ">" $sequence_name\_hits.fas` );

    #
    print "\n\tNumber of Sequences found = $hit_seqs_total\n";
    return $hit_seqs_total;
}

sub run_vsearch {

    my ( $taxa_array_ref, $input_seqs_fname, $sequence_name ) = @_;

    my @taxa_array = @{$taxa_array_ref};
    my $taxa_total = @taxa_array;
    my $taxa_count = 1;

    my $sequence_name_for_blast = $sequence_name;
    $sequence_name_for_blast =~ s/\s+/\_/gms;    # Replace spaces with '_'

    while (@taxa_array) {

        # Current Taxa Name
        my $taxa_name = shift(@taxa_array);

        my $taxa_name_for_blast = $taxa_name;
        $taxa_name_for_blast =~ s/\s+/\_/gms;    # Replace spaces with '_'

        if ( $taxa_name =~ m/^#/sm ) {
            print "\t\tSkipping commented out $taxa_name\n";
            output_report("[INFO]\t$sequence_name: Skipping commented out $taxa_name\n");
            $taxa_count++;
        }
        else {

            # Blast Output Filename
            my $search_output = "$sequence_name\_v\_$taxa_name_for_blast\.$Search_Subprogram";

            # BLAST Output Progress
            printf "\t\t>: $Search_Subprogram: $taxa_count of $taxa_total\n\e[A";    # Progress...

            my $database = $taxa_name_for_blast . '.fas';

            # ublast from usearch package command
            # we will use tabulated output as it's smaller than XML
            # and we don't really need much information other than the hit ID
            my $usearch_command = 'vsearch --usearch_global';
            $usearch_command .= " $sequence_name_for_blast\_query.fas";    # query file
            $usearch_command .= " --db $Seq_Data\/$database";
            $usearch_command .= " --id 0.9";
            $usearch_command .= " --blast6out $search_output";             # output filename
            $usearch_command .= " --threads $Search_Threads";
            $usearch_command .= " --maxaccepts $Search_Tophits";

            system($usearch_command);
            parse_search_output( \@taxa_array, $input_seqs_fname, $sequence_name, $taxa_name, $database );

            $taxa_count++;
        }
    }

    # Find total number of hits
    # I am still relying on 'grep' to count the number of sequences
    # there is no way to get this directly from the Bio::Seq object
    # without needless iteration. Anyone?
    chomp( my $hit_seqs_total = `grep -c ">" $sequence_name\_hits.fas` );

    #
    print "\n\tNumber of Sequences found = $hit_seqs_total\n";
    return $hit_seqs_total;
}

sub run_usearch {

    my ( $taxa_array_ref, $input_seqs_fname, $sequence_name ) = @_;

    my @taxa_array = @{$taxa_array_ref};
    my $taxa_total = @taxa_array;
    my $taxa_count = 1;

    my $sequence_name_for_blast = $sequence_name;
    $sequence_name_for_blast =~ s/\s+/\_/gms;    # Replace spaces with '_'

    while (@taxa_array) {

        # Current Taxa Name
        my $taxa_name = shift(@taxa_array);

        my $taxa_name_for_blast = $taxa_name;
        $taxa_name_for_blast =~ s/\s+/\_/gms;    # Replace spaces with '_'

        if ( $taxa_name =~ m/^#/sm ) {
            print "\t\tSkipping commented out $taxa_name\n";
            output_report("[INFO]\t$sequence_name: Skipping commented out $taxa_name\n");
            $taxa_count++;
        }
        else {

            # Blast Output Filename
            my $search_output = "$sequence_name\_v\_$taxa_name_for_blast\.$Search_Subprogram";

            # BLAST Output Progress
            printf "\t\t>: $Search_Subprogram: $taxa_count of $taxa_total\n\e[A";    # Progress...

            my $database = $taxa_name_for_blast . '.fas';

            # ublast from usearch package command
            # we will use tabulated output as it's smaller than XML
            # and we don't really need much information other than the hit ID
            my $usearch_command = 'usearch -ublast';
            $usearch_command .= " $sequence_name_for_blast\_query.fas";    # query file
            $usearch_command .= " -db $Seq_Data\/$database";
            $usearch_command .= " -evalue $Search_Evalue";                 # this should be a user option eventually
            $usearch_command .= " -blast6out $search_output";              # output filename
            $usearch_command .= " -threads $Search_Threads";
            $usearch_command .= " -maxaccepts $Search_Tophits";

            system($usearch_command);
            parse_search_output( \@taxa_array, $input_seqs_fname, $sequence_name, $taxa_name, $database );

            $taxa_count++;
        }
    }

    # Find total number of hits
    # I am still relying on 'grep' to count the number of sequences
    # there is no way to get this directly from the Bio::Seq object
    # without needless iteration. Anyone?
    chomp( my $hit_seqs_total = `grep -c ">" $sequence_name\_hits.fas` );

    #
    print "\n\tNumber of Sequences found = $hit_seqs_total\n";
    return $hit_seqs_total;
}

# this should be able to read and parse the output from blast/blast+, blat and usearch
# as we are forcing them all to output tabulated data...
sub parse_search_output {

    my ( $taxa_array_ref, $input_seqs_fname, $sequence_name, $taxa_name, $database ) = @_;

    my @taxa_array = @{$taxa_array_ref};

    # blast cannot handle spaces in names (it is also bad practice)
    # but sequence retreival needs an unmodified accession
    my $taxa_name_for_blast = $taxa_name;
    $taxa_name_for_blast =~ s/\s+/\_/gms;    # Replace spaces with '_'

    # search output filename
    my $search_output = "$sequence_name\_v\_$taxa_name_for_blast\.$Search_Subprogram";

    # open file for seq retrieval
    # this will create an index file on the first run, this may slow things down once
    # it may also cause issues if the index is not rebuilt for updated *.fas files
    my $sequence_file;
    if ( $User_Reindex eq 'y' ) {
        $sequence_file = Bio::DB::Fasta->new( "$Seq_Data\/$taxa_name_for_blast\.fas", -reindex );
    }
    else {
        $sequence_file = Bio::DB::Fasta->new("$Seq_Data\/$taxa_name_for_blast\.fas");
    }

    # open file in to searchio stream
    my $read_search_output = Bio::SearchIO->new( -format => 'blasttable', -file => $search_output );

    # searchio for blasttable returns an empty object but not undef
    # so to count 'no hits' - an empty search output - we need a boolean
    my $hit_count = 0;

    # iterate through the stream
    while ( my $result = $read_search_output->next_result ) {
        while ( my $hit = $result->next_hit ) {
            my $hit_name = $hit->name;

            # get the sequence from the file stream
            my $sequence = $sequence_file->seq($hit_name);

            # if there is a problem, don't output an empty sequence
            # report sequence retreival problems to output report
            if ( defined $sequence ) {

                # get the length of the sequence - this was a feature added by Fin originally
                # he was getting some ridiculous sequence lengths 3k++ which were causing
                # many alignment issues
                my $sequence_length = length $sequence;

                if ( $sequence_length <= $Search_Maxlength ) {

                    # output hits to hits files
                    open my $hits_seq_out, '>>', "$sequence_name\_hits.fas";
                    print {$hits_seq_out} ">$hit_name \n$sequence\n";
                    close $hits_seq_out;
                }
                else {
                    output_report("[WARN]\t$sequence_name: sequence hit too long for: $hit_name -> $taxa_name\n");
                }
            }
            else {
                output_report("[WARN]\t$sequence_name: sequence retreival problem for: $hit_name -> $taxa_name\n");
            }

            # if there is a hit = true
            $hit_count = 1;
        }
    }

    # move the files once we are finished to the report directory
    # we may need to make it first
    my $report_dir = "$Working_Dir\/$User_Runid\/report\/$sequence_name";
    if ( !-d $report_dir ) { mkdir $report_dir }
    if ( !-d "$report_dir\/$Search_Subprogram" ) { mkdir "$report_dir\/$Search_Subprogram" }

    system "mv $search_output $report_dir\/$Search_Subprogram";

    return;
}

###########################################################
##           Accessory Subroutines                       ##
###########################################################

sub taxa_list_to_array {

    my $taxa_list_fname = shift;
    open my $infile, '<', $taxa_list_fname;
    chomp( my @data = <$infile> );
    close($infile);

    return @data;
}

sub output_report {

    # Append messages to report file
    my $message   = shift;
    my $file_name = "$Working_Dir\/$User_Runid\_report.txt";
    open my $report, '>>', $file_name;
    print {$report} $message;
    close($report);
    return;
}

sub display_help {

    print "Required files for input:\n\t-s sequence(s) file\n\t-t taxa file\n\t-p paramaters file\n";
    print "Example: perl orchard.pl -s sequences.fasta -t taxa_list.txt -p paramaters.yaml\n";
    print
"Other paramaters:\n\t-b blast only\n\t-a alignment only\n\t-m mask only\n\t-o tree building only\n\t-q run sequentially\n\t-f filter file\n";
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
        my $user_choice = prompt( ' > : ', -yes_no1 );
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

    return;
}

# Handles output of start/end and duration times of different steps
sub timing {

    my $time_operation = shift;
    my $start_time     = shift;

    my $time_now;

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
        print 'Elapsed Time: ' . $format->format_duration($total_time) . "\n";
        output_report( "[INFO]\tElapsed Time: " . $format->format_duration($total_time) . "\n" );
    }

    return $time_now;
}
