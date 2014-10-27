#!/usr/bin/perl
use strict;
use warnings;

#
use autodie;                       # bIlujDI' yIchegh()Qo'; yIHegh()!
use Cwd;                           # Gets pathname of current working directory
use Bio::SeqIO;                    # Use Bio::Perl
use Digest::MD5;                   # Generate random string for run ID
use English qw(-no_match_vars);    # No magic perl variables!
use File::Basename;                # Remove path information and extract 8.3 filename
use Getopt::Std;                   # Command line options, finally!
use feature qw{ switch };          # Given/when instead of switch
use IO::Prompt;                    # User prompts
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
    my $taxa_list_fname = "$options{t}";
    my @taxa_array      = taxa_list_to_array($taxa_list_fname);

    my $input_seqs_fname = "$options{s}";

    my $paramaters = LoadFile("$options{p}");

    # read in parameters from YAML file and/or set defaults
    # user options
    # modify the md5_hex to ten chars from pos 0 if none given in YAML
    my $user_options = $paramaters->{user}->{run_id} || substr( Digest::MD5::md5_hex(rand), 0, 10 );
    my $run_directory = "$WORKING_DIR\/$user_options";

    # Now we can make the directories
    setup_main_directories($run_directory);

    # Report
    output_report("Run ID: $user_options\nDirectory: $run_directory\n");

    # search options
    my $search_program    = $paramaters->{search}->{program}    || 'blast+';
    my $search_subprogram = $paramaters->{search}->{subprogram} || 'blastp';
    my $search_evalue     = $paramaters->{search}->{evalue}     || '1e-10';
    my $search_tophits    = $paramaters->{search}->{top_hits}   || '1';
    my $search_maxlength  = $paramaters->{search}->{max_length} || '3000';
    my $search_special_taxa    = $paramaters->{search}->{special_taxa};        # no default
    my $search_special_tophits = $paramaters->{search}->{special_top_hits};    # no default
    my $search_threads         = $paramaters->{search}->{threads} || '1';
    my $search_other           = $paramaters->{search}->{other};               # no default

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
    # this needs to quit the program if they are not set...
    my $username  = $paramaters->{database}->{user};
    my $password  = $paramaters->{database}->{password};
    my $server_ip = $paramaters->{database}->{server};
    my $database  = $paramaters->{database}->{database};
    my $tablename = $paramaters->{database}->{tablename};

    # directory options
    my $programs = $paramaters->{directories}->{location} || '/usr/bin';
    my $seq_data = $paramaters->{directories}->{database};    # no default

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

    # run b-a-m-t for each sequence sequentially
    if ( $options{q} ) {
        ## experimental, is it really useful?
    }

    # run default
    if ( !$options{b} && !$options{a} && !$options{m} && !$options{o} && !$options{q} ) {

        # run all steps but all blasts first then amt steps
        print "Running: ALL Steps, all searches ($search_program) first!\n";
        search_step(
            $search_program,         $search_subprogram, \@taxa_array,      $input_seqs_fname,
            $search_evalue,          $search_tophits,    $search_maxlength, $search_special_taxa,
            $search_special_tophits, $search_threads,    $search_other
        );
    }
}
else {
    display_help();
}

###########################################################
##           Main Subroutines 	                         ##
###########################################################

sub search_step {

    my (
        $search_program,         $search_subprogram, $taxa_array_ref,   $input_seqs_fname,
        $search_evalue,          $search_tophits,    $search_maxlength, $search_special_taxa,
        $search_special_tophits, $search_threads,    $search_other
    ) = @_;
    my @taxa_array = @{$taxa_array_ref};

    my $input_seqs_count = 1;

    # open bioperl seqio object with user input sequences
    my $seq_in = Bio::SeqIO->new( -file => "<$input_seqs_fname" );

    # I am still relying on 'grep' to count the number of sequences
    # there is no way to get this directly from the Bio::Seq object
    # without needless iteration. Anyone?
    chomp( my $input_seqs_total = `grep -c ">" $input_seqs_fname` );

    # iterate through the stream of sequences and perform searches
    output_report("Starting $input_seqs_total searches using $search_program($search_subprogram)\n");
    while ( my $seq = $seq_in->next_seq() ) {

        print "Processing: $input_seqs_count of $input_seqs_total\n";

        given ($search_program) {
            when (/BLAST\+/ism) {

                my $num_seqs = run_blast_plus(
                    $search_program,         $search_subprogram, \@taxa_array,      $input_seqs_fname,
                    $search_evalue,          $search_tophits,    $search_maxlength, $search_special_taxa,
                    $search_special_tophits, $search_threads,    $search_other
                );
                print "\tRunning: blast+\n";
            }
            when (/BLAST/ism) {

                #$num_seqs = run_blast( "$sequence_name", "$sequence_name\_test.fas", "$sequence_name\_seqs.fas" );
                print "\tRunning: legacy blast\n";
            }
            when (/BLAT/ism) {

                #$num_seqs = run_blat( "$sequence_name", "$sequence_name\_test.fas", "$sequence_name\_seqs.fas" );
            }
            when (/USEARCH/ism) {

                #$num_seqs = run_ublast( "$sequence_name", "$sequence_name\_test.fas", "$sequence_name\_seqs.fas" );
            }
            default {

                #$num_seqs = run_blast( "$sequence_name", "$sequence_name\_test.fas", "$sequence_name\_seqs.fas" );
            }
        }
        $input_seqs_count++;
    }

    return;
}

sub run_blast_plus {

    my (
        $search_program,         $search_subprogram, $taxa_array_ref,   $input_seqs_fname,
        $search_evalue,          $search_tophits,    $search_maxlength, $search_special_taxa,
        $search_special_tophits, $search_threads,    $search_other
    ) = @_;
    my @taxa_array = @{$taxa_array_ref};

    # blast(x) from blast+ package command
    #my $blast_command = "blastp -task $search_program_blast";
    #$blast_command .= " -db $seq_data\/XXX";
    #$blast_command .= " -query XXX";
    #$blast_command .= " -out XXX";
    #$blast_command .= " -evalue $search_evalue";
    #$blast_command .= " -outfmt XXX";
    #$blast_command .= " -num_alignments 0";
    #$blast_command .= " -max_target_seqs $search_tophits";
    #$blast_command .= " -num_threads $search_threads";
    #$blast_command .= " $search_other";

    #system($blast_command);
}

sub run_legacy_blast {

    # blast(x) from legacy blast package command
    #my $blast_command = "blastall -p $search_program_blast";
    #$blast_command .= " -d $seq_data\/XXX";
    #$blast_command .= " -i XXX";
    #$blast_command .= " -o XXX";
    #$blast_command .= " -e $search_evalue";
    #$blast_command .= " -m XXX";
    #$blast_command .= " -b 0";
    #$blast_command .= " -v $search_tophits";
    #$blast_command .= " -a $search_threads";
    #$blast_command .= " $search_other";

    #system($blast_command);
}

sub run_blat {

    #my $blat_command = "blat";
    #$blat_command .= " -prot"; # this should be a user option eventually
    #$blat_command .= " $seq_data\/XXX"; #database
    #$blat_command .= " XXX"; # query file
    #$blat_command .= " -out=blast";
    #$blat_command .= " XXX"; # output filename
}

sub run_usearch {

    #my $usearch_command = "usearch -ublast";
    #$blat_command .= " XXX"; # query file
    #$blat_command .= " -db $seq_data\/XXX";
    #$blat_command .= " -evalue $search_evalue"; # this should be a user option eventually
    #$blat_command .= " -blast6out XXX; # output filename
    #$blat_command .= " -threads $search_threads";
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
    my $file_name = "$WORKING_DIR\/report.txt";
    open my $report, ">>", $file_name;
    print $report $message;
    close($report);
    return;
}

sub display_help {

    print "Required files for input:\n\t-s sequence(s) file\n\t-t taxa file\n\t-p paramaters file\n";
    print "Example: perl tree_pipe.pl -s sequences.fasta -t taxa_list.txt -p paramaters.yaml\n";
    print
"Other paramaters:\n\t-b blast only\n\t-a alignment only\n\t-m mask only\n\t-t tree building only\n\t-q run sequentially\n";
    exit(1);
}

# this checks to see if the directories needed for file output
# are available, if not it creates them all - for now.
sub setup_main_directories {

    my $run_directory = shift;

    # main directories
    my $seqs_dir = "$run_directory\/seqs";
    my $alig_dir = "$run_directory\/alignments";
    my $mask_dir = "$run_directory\/masks";
    my $tree_dir = "$run_directory\/trees";
    my $excl_dir = "$run_directory\/excluded";
    my $repo_dir = "$run_directory\/report";

    if ( !-d $run_directory ) {
        output_report("Creating Directory $run_directory\n");
        mkdir $run_directory;

        # create sub-directories
        output_report("Creating Subdirectories\n");

        # create directories if they don't exist!
        if ( !-d $seqs_dir ) { mkdir $seqs_dir }
        if ( !-d $alig_dir ) { mkdir $alig_dir }
        if ( !-d $mask_dir ) { mkdir $mask_dir }
        if ( !-d $tree_dir ) { mkdir $tree_dir }
        if ( !-d $excl_dir ) { mkdir $excl_dir }
        if ( !-d $repo_dir ) { mkdir $repo_dir }
    }
    else {
        print "Directory Already Exists!\nContinue anyway? y/n\n";
        my $user_choice = prompt( ">: ", -yes_no1 );
        if ( $user_choice =~ m/n/ism ) {

            output_report("Terminating, run directories already exist\n");
            exit;
        }
        else {
            output_report("Continuing, even though run directories already exist\n");

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
