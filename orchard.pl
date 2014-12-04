#!/usr/bin/perl
use strict;
use warnings;

# Given/when instead of switch - warns in 5.18 eventually I will switch this to "for()" http://www.effectiveperlprogramming.com/2011/05/use-for-instead-of-given/
#use experimental 'smartmatch';
no if $] >= 5.017011, warnings => 'experimental::smartmatch';    # ignore experimental warning for 'when'
use autodie;                                                     # bIlujDI' yIchegh()Qo'; yIHegh()!
use Bio::DB::Fasta;                                              # for indexing fasta files for sequence retrieval
use Bio::SearchIO;                                               # parse blast tabular format
use Bio::SeqIO;                                                  # Use Bio::Perl
use Cwd;                                                         # Gets pathname of current working directory
use DateTime::Format::Duration;                                  # Duration of processes
use DateTime;                                                    # Start and End times
use Digest::MD5;                                                 # Generate random string for run ID
use Digest::MD5 'md5_hex';
use feature qw{ switch };
use File::Basename;                                              # Remove path information and extract 8.3 filename
use Getopt::Std;                                                 # Command line options, finally!
use IO::Prompt;                                                  # User prompts
use YAML::XS qw/LoadFile/;                                       # for the parameters file, user friendly layout

##
use Data::Dumper;                                                # temporary during rewrite to dump data nicely to screen

# remove ## comments before publication
#
our $WORKING_DIR = getcwd();
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
our $EMPTY             = q{};
our $ALIGNMENT_OPTIONS = $EMPTY;
our $ALIGNMENT_PROGRAM = $EMPTY;
our $ALIGNMENT_THREADS = $EMPTY;

our $MASKING_CUTOFF1 = $EMPTY;
our $MASKING_CUTOFF2 = $EMPTY;

our $SEARCH_EVALUE    = $EMPTY;
our $SEARCH_MAXLENGTH = $EMPTY;

#our $SEARCH_OTHER           = $EMPTY;
our $SEARCH_PROGRAM         = $EMPTY;
our @SEARCH_SPECIAL_TAXA    = $EMPTY;
our $SEARCH_SPECIAL_TOPHITS = $EMPTY;
our $SEARCH_SUBPROGRAM      = $EMPTY;
our $SEARCH_THREADS         = $EMPTY;
our $SEARCH_TOPHITS         = $EMPTY;

our $SEQ_DATA = $EMPTY;

our $TREE_MINTAXA = $EMPTY;
our $TREE_OPTIONS = $EMPTY;
our $TREE_PROGRAM = $EMPTY;

our $USER_REINDEX = $EMPTY;
our $USER_RUNID   = $EMPTY;

###########################################################
##           Main Program Flow                           ##
###########################################################

# declare the perl command line flags/options we want to allow
my %options = ();
getopts( 's:t:p:hvbamof', \%options ) or display_help();    # or display_help();

# Display the help message if the user invokes -h
if ( $options{h} ) { display_help(); }
if ( $options{v} ) { print "Orchard $VERSION\n"; }

if ( defined $options{p} && defined $options{t} && defined $options{s} ) {
    my $taxa_list_fname = "$options{t}";
    my @taxa_array      = taxa_list_to_array($taxa_list_fname);

    my $input_seqs_fname = "$options{s}";

    # read in parameters from YAML file and/or set defaults
    # user options
    # modify the md5_hex to ten chars from pos 0 if none given in YAML
    my $paramaters = LoadFile("$options{p}");
    $USER_RUNID = $paramaters->{user}->{run_id} || substr( Digest::MD5::md5_hex(rand), 0, 10 );
    $USER_REINDEX = $paramaters->{user}->{reindex} || 'n';    # default no
    if ( $USER_REINDEX eq 'y' ) { output_report("[INFO]\tUser requested reindexing of database directory files\n"); }
    my $run_directory = "$WORKING_DIR\/$USER_RUNID";

    # Now we can make the directories
    setup_main_directories( $run_directory, $options{f} );

    # Report
    output_report("[INFO]\tRun ID: $USER_RUNID\n[INFO]\tDirectory: $run_directory\n");

    # search options
    $SEARCH_PROGRAM    = $paramaters->{search}->{program}    || 'blast+';
    $SEARCH_SUBPROGRAM = $paramaters->{search}->{subprogram} || 'blastp';
    $SEARCH_EVALUE     = $paramaters->{search}->{evalue}     || '1e-10';
    $SEARCH_TOPHITS    = $paramaters->{search}->{top_hits}   || '1';
    $SEARCH_MAXLENGTH  = $paramaters->{search}->{max_length} || '3000';
    @SEARCH_SPECIAL_TAXA = split( /,/ms, $paramaters->{search}->{special_taxa} );    # no default
    $SEARCH_SPECIAL_TOPHITS = $paramaters->{search}->{special_top_hits} || $SEARCH_TOPHITS;
    $SEARCH_THREADS         = $paramaters->{search}->{threads}          || '1';

    #$SEARCH_OTHER = $paramaters->{search}->{other};                                # no default

    # alignment options
    $ALIGNMENT_PROGRAM = $paramaters->{alignment}->{program} || 'mafft';
    $ALIGNMENT_OPTIONS = $paramaters->{alignment}->{options};                        # no default
    $ALIGNMENT_THREADS = $paramaters->{alignment}->{threads} || '1';

    # masking options
    $MASKING_CUTOFF1 = $paramaters->{masking}->{cutoff_1} || '50';
    $MASKING_CUTOFF2 = $paramaters->{masking}->{cutoff_2} || '20';

    # tree building options
    $TREE_PROGRAM = $paramaters->{trees}->{program} || 'FastTreeMP';
    $TREE_OPTIONS = $paramaters->{trees}->{options};                                 # no default
    $TREE_MINTAXA = $paramaters->{trees}->{min_taxa} || '4';

    # directory options
    $SEQ_DATA = $paramaters->{directories}->{database};                              # no default

    # only run search (blast) step
    if ( $options{b} ) {

        # Run the user selected sequence search option
        my $start_time = timing('start');

        if ( -d $input_seqs_fname ) {
            print "Running: Search ($SEARCH_PROGRAM) from Ortho Groups Directory ONLY\n";

            search_step_ortho_groups( \@taxa_array, $input_seqs_fname );
        }
        else {
            print "Running: Search ($SEARCH_PROGRAM) ONLY\n";
            search_step( \@taxa_array, $input_seqs_fname );
        }
        my $end_time = timing( 'end', $start_time );
    }

    # only run alignment step
    if ( $options{a} ) {
        print "Running: Alignment ($ALIGNMENT_PROGRAM) ONLY\n";
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
        print "Running: Tree Reconstruction with $TREE_PROGRAM ONLY\n";
        my $start_time = timing('start');
        tree_step();
        my $end_time = timing( 'end', $start_time );
    }

    # run default
    if ( !$options{b} && !$options{a} && !$options{m} && !$options{o} && !$options{q} ) {

        # run all steps, but all blasts first then amt steps
        print "Running: ALL Steps; all searches ($SEARCH_PROGRAM) first!\n";

        # Run the user selected sequence search option
        my $start_time = timing('start');
        search_step( \@taxa_array, $input_seqs_fname );
        my $end_time = timing( 'end', $start_time );

        # Run the user selected alignment option
        print "Running: Alignment with $ALIGNMENT_PROGRAM\n";
        $start_time = timing('start');
        alignment_step();
        $end_time = timing( 'end', $start_time );

        # Run the user selected masking option
        print "Running: Masking with trimal\n";
        $start_time = timing('start');
        masking_step();
        $end_time = timing( 'end', $start_time );

        # Run the user selected tree building options
        print "Running: Tree Reconstruction with $TREE_PROGRAM\n";
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
##           Tree Subroutines                  ##
#################################################

sub tree_step {

    my $masks_directory = "$WORKING_DIR\/$USER_RUNID\/masks";
    my $tree_directory  = "$WORKING_DIR\/$USER_RUNID\/trees";

    # get list of sequences that were not excluded in previous stage
    my @mask_file_names = glob "$masks_directory\/*.afa-tr";

    output_report("[INFO]\tTREE BUILDING: $#mask_file_names alignments using $TREE_PROGRAM\n");

    for my $i ( 0 .. $#mask_file_names ) {
        my ( $file, $dir, $ext ) = fileparse $mask_file_names[$i], '\.afa-tr';

        run_fasttree("$masks_directory\/$file\.afa-tr");
    }

    return;
}

sub run_fasttree {
    my $masked_sequences = shift;
    my $fasttree_command = $EMPTY;

    my ( $file, $dir, $ext ) = fileparse $masked_sequences, '\.afa\-tr';

    print "\tCalculating FastTree on $file\n";

    if ( $TREE_PROGRAM =~ /fasttreemp/ism ) {
        $fasttree_command = "FastTreeMP $TREE_OPTIONS";
        $fasttree_command .= " $masked_sequences";
        $fasttree_command .= " > $WORKING_DIR\/$USER_RUNID\/trees\/$file\_FT.tree";
    }
    else {
        $fasttree_command = "FastTree $TREE_OPTIONS";
        $fasttree_command .= " $masked_sequences";
        $fasttree_command .= " > $WORKING_DIR\/$USER_RUNID\/trees\/$file\_FT.tree";
    }

    system($fasttree_command);
    return;
}

#################################################
##           Masking Subroutines               ##
#################################################

sub masking_step {

    # directory variables
    my $alignments_directory = "$WORKING_DIR\/$USER_RUNID\/alignments";
    my $masks_directory      = "$WORKING_DIR\/$USER_RUNID\/masks";

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
        if ( $mask_length <= $MASKING_CUTOFF1 ) {

            # then re-run trimal and report mask length with -automated1 option
            $mask_length = run_trimal( $current_sequences, '-automated1' );

            print "Mask Length is ($mask_length) and is";

            # if the length is less than the second limit (always the smaller)
            if ( $mask_length <= $MASKING_CUTOFF2 ) {
                print " not OK. Excluding sequence.\n";

                # then abandon this sequence, report it, move to excluded
                output_report("[WARN]\t$file does not satisfy trimal cutoffs ($MASKING_CUTOFF1 or $MASKING_CUTOFF2). Moved to excluded directory.\n");
                system "mv $WORKING_DIR\/$USER_RUNID\/masks\/$file\.afa-tr $WORKING_DIR\/$USER_RUNID\/excluded\/$file\.afa-tr";
            }
            else {
                print " OK.\n";
            }
        }
        else {
            print "\t\tMask Length: $mask_length is OK\n";
        }
    }

    return;
}

sub run_trimal {

    my $aligned_sequences = shift;
    my $option            = shift;
    my $masks_directory   = "$WORKING_DIR\/$USER_RUNID\/masks";

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

    # directory variables
    my $sequence_directory   = "$WORKING_DIR\/$USER_RUNID\/seqs";
    my $alignments_directory = "$WORKING_DIR\/$USER_RUNID\/alignments";

    # get list of sequences that were not excluded in previous stage
    my @seq_file_names = glob "$sequence_directory\/*.fas";

    # iterate through the stream of sequences and perform each search
    output_report("[INFO]\tALIGNING: $#seq_file_names alignments using $ALIGNMENT_PROGRAM\n");

    for my $i ( 0 .. $#seq_file_names ) {
        my $current_sequences = $seq_file_names[$i];

        my ( $file, $dir, $ext ) = fileparse $current_sequences, '\.fas';

        for ($ALIGNMENT_PROGRAM) {
            when (/mafft/ism) {
                my $mafft_command = "mafft $ALIGNMENT_OPTIONS";
                $mafft_command .= " --thread $ALIGNMENT_THREADS";
                $mafft_command .= " $current_sequences > $alignments_directory\/$file\.afa";

                system($mafft_command);
            }
            when (/muscle/ism) {

                my $muscle_command = "muscle $ALIGNMENT_OPTIONS";
                $muscle_command .= " -in $current_sequences -out $alignments_directory\/$file\.afa";

                system($muscle_command);
            }
        }
    }

    return;
}

#################################################
##           Search Subroutines                ##
#################################################
sub search_step_ortho_groups {

    # get all of the values passed to the sub...
    my ( $taxa_array_ref, $input_seqs_fname ) = @_;

    # files from orthogroup/multi-seq directory
    my @file_names        = glob "$input_seqs_fname\/*.fasta";
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

        our ( $ortho_file, $ortho_dir, $ortho_ext ) = fileparse( $orthofile, ".fasta" );

        #
        my $num_hit_seqs  = 0;
        my $sequence_name = $EMPTY;

        # open bioperl seqio object with user input sequences
        my $seq_in = Bio::SeqIO->new( -file => "<$orthofile" );

        # I am still relying on 'grep' to count the number of sequences
        # there is no way to get this directly from the Bio::Seq object
        # without needless iteration. Anyone?
        chomp( my $input_seqs_total = `grep -c ">" $orthofile` );

        # iterate through the stream of sequences and perform each search
        output_report("[INFO]\tStarting $input_seqs_total searches using $SEARCH_PROGRAM($SEARCH_SUBPROGRAM)\n");

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

            for ($SEARCH_PROGRAM) {
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
                when (/USEARCH/ism) {

                    print "\tRunning: usearch\n";
                    run_usearch( \@taxa_array, "$sequence_name\_query.fas", $sequence_name );
                }
                default {

                    print "\tRunning: (default) blast+ on $sequence_name\n";
                    run_blast_plus( \@taxa_array, "$sequence_name\_query.fas", $sequence_name );
                }
            }
            $input_seqs_count++;

            # remove the query hits after each run as they are now stored
            # in the *ortho_seqs.fas file
            unlink "$WORKING_DIR\/$sequence_name\_query.fas";

            # append the current seed sequence's hits to the ortho_seqs.fas file
            system "cat $WORKING_DIR\/$sequence_name\_hits.fas >> $WORKING_DIR\/$ortho_file\_ortho\_hits.fas";

            # then get rid of the query file as we don't need it anymore
            unlink "$WORKING_DIR\/$sequence_name\_hits.fas";
        }
        $ortho_files_count++;

        # now we can remove duplicate sequences based on accession only
        remove_duplicate_sequences("$WORKING_DIR\/$ortho_file\_ortho\_hits.fas");

        # get the number of sequences in the file
        chomp( $num_hit_seqs = `grep -c ">" $WORKING_DIR\/$ortho_file\_ortho\_hits.fas` );

        if ( $num_hit_seqs <= $TREE_MINTAXA ) {
            output_report("[WARN]\t$ortho_file: Too few hits\n");
            unlink "$WORKING_DIR\/$ortho_file\_ortho\_hits.fas";
            system "mv $WORKING_DIR\/$ortho_file\_ortho\_non\_redundant\_hits.fas $WORKING_DIR\/$USER_RUNID\/excluded\/";
        }
        else {
            unlink "$WORKING_DIR\/$ortho_file\_ortho\_hits.fas";
            system "mv $WORKING_DIR\/$ortho_file\_ortho\_non\_redundant\_hits.fas $WORKING_DIR\/$USER_RUNID\/seqs\/$ortho_file\_hits.fas";

            # remove selenocystein and non-alpha numeric characters as they cause BLAST/MAFFT
            # to complain/crash (and the --anysymbol option produces terrible alignments)
            # I have to check for this as some genome projects are just full of junk!
            system "sed -i \'/^>/! s/U|\\w/X/g\' $WORKING_DIR\/$USER_RUNID\/seqs\/$ortho_file\_hits.fas";
        }
    }

    return;
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
        - file   => ">$WORKING_DIR\/$FILE\_non_redundant$EXT",
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

    print "XXXX $input_seqs_fname XXXX\n";

    # dereference the array
    my @taxa_array = @{$taxa_array_ref};

    # we should always have one sequence, right!?
    my $input_seqs_count = 1;

    #
    my $num_hit_seqs  = 0;
    my $sequence_name = $EMPTY;

    # open bioperl seqio object with user input sequences
    my $seq_in = Bio::SeqIO->new( -file => "<$input_seqs_fname" );

    # I am still relying on 'grep' to count the number of sequences
    # there is no way to get this directly from the Bio::Seq object
    # without needless iteration. Anyone?
    chomp( my $input_seqs_total = `grep -c ">" $input_seqs_fname` );

    # iterate through the stream of sequences and perform each search
    output_report("[INFO]\tStarting $input_seqs_total searches using $SEARCH_PROGRAM($SEARCH_SUBPROGRAM)\n");
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

        for ($SEARCH_PROGRAM) {
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
            when (/USEARCH/ism) {

                print "\tRunning: usearch\n";
                $num_hit_seqs = run_usearch( \@taxa_array, $input_seqs_fname, $sequence_name );
            }
            default {

                print "\tRunning: (default) blast+ on $sequence_name\n";
                $num_hit_seqs = run_blast_plus( \@taxa_array, $input_seqs_fname, $sequence_name );
            }
        }
        $input_seqs_count++;

        if ( $num_hit_seqs <= $TREE_MINTAXA ) {
            output_report("[WARN]\t$sequence_name: Too few hits\n");
            unlink "$WORKING_DIR\/$sequence_name\_query.fas";
            system "mv $sequence_name\_hits.fas $WORKING_DIR\/$USER_RUNID\/excluded\/";
        }
        else {
            unlink "$WORKING_DIR\/$sequence_name\_query.fas";
            system "mv $sequence_name\_hits.fas $WORKING_DIR\/$USER_RUNID\/seqs\/";

            # remove selenocystein and non-alpha numeric characters as they cause BLAST/MAFFT
            # to complain/crash (and the --anysymbol option produces terrible alignments)
            # I have to check for this as some genome projects are just full of junk!
            system "sed -i \'/^>/! s/U|\\w/X/g\' $WORKING_DIR\/$USER_RUNID\/seqs\/$sequence_name\_hits.fas";
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
            my $search_output = "$sequence_name\_v\_$taxa_name_for_blast\.$SEARCH_SUBPROGRAM";

            # BLAST Output Progress
            printf "\t\t>: $SEARCH_SUBPROGRAM: $taxa_count of $taxa_total\n\e[A";    # - $taxa_name\n\e[A";    # Progress...

            my $database = $taxa_name_for_blast . '.fas';

            # In the future I might consider replacing grep with List:MoreUtils 'any' to save
            # large searches
            if ( grep { $_ eq $taxa_name } @SEARCH_SPECIAL_TAXA ) {

                # blast(x) from blast+ package command
                # we will use tabulated output as it's smaller than XML
                # and we don't really need much information other than the hit ID
                my $blast_command = "$SEARCH_SUBPROGRAM -task $SEARCH_SUBPROGRAM";
                $blast_command .= " -db $SEQ_DATA\/$database";
                $blast_command .= " -query $sequence_name_for_blast\_query.fas";
                $blast_command .= " -out $search_output";
                $blast_command .= " -evalue $SEARCH_EVALUE";
                $blast_command .= ' -outfmt 6';
                $blast_command .= " -max_target_seqs $SEARCH_SPECIAL_TOPHITS";
                $blast_command .= " -num_threads $SEARCH_THREADS";

                #$blast_command .= " $SEARCH_OTHER";

                system($blast_command);
                parse_search_output( \@taxa_array, $input_seqs_fname, $sequence_name, $taxa_name, $database );

                output_report("[INFO]\t$sequence_name: Using special -max_target_seqs $SEARCH_SPECIAL_TOPHITS for $taxa_name\n");
            }
            else {
                # blast(x) from blast+ package command
                # we will use tabulated output as it's smaller than XML
                # and we don't really need much information other than the hit ID
                my $blast_command = "$SEARCH_SUBPROGRAM -task $SEARCH_SUBPROGRAM";
                $blast_command .= " -db $SEQ_DATA\/$database";
                $blast_command .= " -query $sequence_name_for_blast\_query.fas";
                $blast_command .= " -out $search_output";
                $blast_command .= " -evalue $SEARCH_EVALUE";
                $blast_command .= ' -outfmt 6';
                $blast_command .= " -max_target_seqs $SEARCH_TOPHITS";
                $blast_command .= " -num_threads $SEARCH_THREADS";

                #$blast_command .= " $SEARCH_OTHER";

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
            my $search_output = "$sequence_name\_v\_$taxa_name_for_blast\.$SEARCH_SUBPROGRAM";

            # BLAST Output Progress
            printf "\t\t>: $SEARCH_SUBPROGRAM: $taxa_count of $taxa_total\n\e[A";    # Progress...

            my $database = $taxa_name_for_blast . '.fas';

            # blast(x) from legacy blast package command
            # we will use tabulated output as it's smaller than XML
            # and we don't really need much information other than the hit ID
            my $blast_command = "blastall -p $SEARCH_SUBPROGRAM";
            $blast_command .= " -d $SEQ_DATA\/$database";
            $blast_command .= " -i $sequence_name_for_blast\_query.fas";
            $blast_command .= " -o $search_output";
            $blast_command .= " -e $SEARCH_EVALUE";
            $blast_command .= ' -m 8';
            $blast_command .= " -b $SEARCH_TOPHITS";

            #$blast_command .= " -v $SEARCH_TOPHITS";
            $blast_command .= " -a $SEARCH_THREADS";

            #$blast_command .= " $SEARCH_OTHER";

            system($blast_command);
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
            my $search_output = "$sequence_name\_v\_$taxa_name_for_blast\.$SEARCH_SUBPROGRAM";

            # BLAST Output Progress
            printf "\t\t>: $SEARCH_SUBPROGRAM: $taxa_count of $taxa_total\n\e[A";    # Progress...

            my $database = $taxa_name_for_blast . '.fas';

            # blat from blat package command
            # we will use tabulated output as it's smaller than XML
            # and we don't really need much information other than the hit ID
            my $blat_command = 'blat';
            $blat_command .= ' -prot';                                               # this should be a user option eventually
            $blat_command .= " $SEQ_DATA\/$database";
            $blat_command .= " $sequence_name_for_blast\_query.fas";                 # query file
            $blat_command .= ' -out=blast8';
            $blat_command .= " $search_output";                                      # output filename

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
            my $search_output = "$sequence_name\_v\_$taxa_name_for_blast\.$SEARCH_SUBPROGRAM";

            # BLAST Output Progress
            printf "\t\t>: $SEARCH_SUBPROGRAM: $taxa_count of $taxa_total\n\e[A";    # Progress...

            my $database = $taxa_name_for_blast . '.fas';

            # ublast from usearch package command
            # we will use tabulated output as it's smaller than XML
            # and we don't really need much information other than the hit ID
            my $usearch_command = 'usearch -ublast';
            $usearch_command .= " $sequence_name_for_blast\_query.fas";              # query file
            $usearch_command .= " -db $SEQ_DATA\/$database";
            $usearch_command .= " -evalue $SEARCH_EVALUE";                           # this should be a user option eventually
            $usearch_command .= " -blast6out $search_output";                        # output filename
            $usearch_command .= " -threads $SEARCH_THREADS";
            $usearch_command .= " -maxaccepts $SEARCH_TOPHITS";

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
    my $search_output = "$sequence_name\_v\_$taxa_name_for_blast\.$SEARCH_SUBPROGRAM";

    # open file for seq retrieval
    # this will create an index file on the first run, this may slow things down once
    # it may also cause issues if the index is not rebuilt for updated *.fas files
    my $sequence_file;
    if ( $USER_REINDEX eq 'y' ) {
        $sequence_file = Bio::DB::Fasta->new( "$SEQ_DATA\/$taxa_name_for_blast\.fas", -reindex );
    }
    else {
        $sequence_file = Bio::DB::Fasta->new("$SEQ_DATA\/$taxa_name_for_blast\.fas");
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

            ##print "\t\t\tHit: " . $hit->name . " :\n" if $VERBOSE == 1;

            # get the sequence from the file stream
            my $sequence = $sequence_file->seq($hit_name);

            # if there is a problem, don't output an empty sequence
            # report sequence retreival problems to output report
            if ( defined $sequence ) {

                # get the length of the sequence - this was a feature added by Fin originally
                # he was getting some ridiculous sequence lengths 3k++ which were causing
                # many alignment issues
                my $sequence_length = length $sequence;

                if ( $sequence_length <= $SEARCH_MAXLENGTH ) {

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
    my $report_dir = "$WORKING_DIR\/$USER_RUNID\/report\/$sequence_name";
    if ( !-d $report_dir ) { mkdir $report_dir }
    if ( !-d "$report_dir\/$SEARCH_SUBPROGRAM" ) { mkdir "$report_dir\/$SEARCH_SUBPROGRAM" }

    system "mv $search_output $report_dir\/$SEARCH_SUBPROGRAM";

    # return
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
    my $file_name = "$WORKING_DIR\/$USER_RUNID\_report.txt";
    open my $report, '>>', $file_name;
    print {$report} $message;
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
        print 'Elapsed Time: ' . $format->format_duration($total_time) . "\n";
        output_report( "[INFO]\tElapsed Time: " . $format->format_duration($total_time) . "\n" );
    }

    return $time_now;
}
