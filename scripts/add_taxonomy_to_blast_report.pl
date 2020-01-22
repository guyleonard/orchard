#!/usr/bin/env perl
use strict;
use warnings;

use DBI;
use File::Basename;

my $current_file = $ARGV[0];
my $db           = $ARGV[1];

my $orchardDB_username = 'orcharddb';
my $orchardDB_password = 'jcsy4s8b';

open( my $fh_in, '<', $current_file )
    or die "Cannot open $current_file: $!\n";

my $outfile = fileparse $current_file;
$outfile = "$outfile\_taxonomy.out";

open( my $fh_out, '>', $outfile )
    or die "Cannot open $outfile: $!\n";

while ( my $line = <$fh_in> ) {
    chomp($line);

    my $accession = ( split /\t/, $line )[1];

    my $genus_species = rename_sequences( $accession, $db );

    print $fh_out "$line\t$genus_species\n";
}

#
# Renaming
#
sub rename_sequences {
    my $accession = shift;
    my $db        = shift;
    my $dbh       = sqllite($db);

    my $statement
        = qq(SELECT odb_accessions.extracted_accession, odb_maintable.genus_species, odb_maintable.source FROM odb_maintable INNER JOIN odb_accessions ON odb_maintable.genome_id = odb_accessions.lookup_id WHERE hashed_accession='$accession');
    my $prepare = $dbh->prepare($statement);
    $prepare->execute();

    my ( $extracted_accession, $genus_species, $source )
        = $prepare->fetchrow_array();

    return $genus_species;
}

#
# SQLLite
#
sub sqllite {
    my $db_name = shift;
    $db_name = fileparse $db_name;

    my $driver   = 'SQLite';
    my $database = "$db\/$db_name\.sqlite";
    my $dsn      = "DBI:$driver:dbname=$database";

    my $dbh
        = DBI->connect( $dsn, $orchardDB_username, $orchardDB_password,
        { RaiseError => 1 } )
        or die $DBI::errstr;

    return $dbh;
}
