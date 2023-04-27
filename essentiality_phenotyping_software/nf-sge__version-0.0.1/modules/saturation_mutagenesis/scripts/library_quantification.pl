#! /usr/bin/env perl

use strict;
use Data::Dumper;
use Getopt::Long;
use File::Basename;
use Bio::SeqIO;

my ( $library_file, $output_directory, $sample_fastq, $wildtype_seq );

GetOptions ( "library=s"    => \$library_file,
             "sample=s"     => \$sample_fastq,
             "wildtype=s"    => \$wildtype_seq,
             "outdir=s"     => \$output_directory )
    or die("Error in command line arguments\n");

die "Cannot find library file: $library_file\n" unless ( -e $library_file );
die "Cannot find sample file: $sample_fastq\n" unless ( -e $sample_fastq );
die "Cannot find output_directory: $output_directory\n" unless ( -d $output_directory );

#########################################
#          GET SAMPLE NAME              #
#########################################

my $sample_name = basename( $sample_fastq );
$sample_name =~ s/\.fastq\.gz//;

#########################################
#          OUTPUT FILENAMES             #
#########################################

my $count_file = $sample_name . ".counts.tsv";
my $stats_file = $sample_name . ".stats.tsv";
my $ambiguous_file = $sample_name . ".ambiguous.tsv";
my $unmapped_file = $sample_name . ".unmapped.tsv";

if ( defined $output_directory && $output_directory ne "" ) {
    $count_file = $output_directory . "/" . $count_file;
    $stats_file = $output_directory . "/" . $stats_file;
    $ambiguous_file = $output_directory . "/" . $ambiguous_file;
    $unmapped_file = $output_directory . "/" . $unmapped_file;
}

#########################################
#           PARSE LIBRARY               #
#########################################

print "Parsing library: $library_file\n";

my %library;
open( my $REF, "<", $library_file ) or die "Can't open $library_file:$\n";
while ( my $line = <$REF> ) {
    chomp $line;
    my ( $id, $seq ) = split( "\t", $line );
	$seq = uc($seq);
    if ( exists $library{ $seq } ) {
        push( @{ $library { $seq }{'id'} }, $id );
    } else {
        $library{$seq}{'id'}[0] = $id;
        $library{$seq}{'count'} = 0;
    }
}
close( $REF );

#########################################
#            VALIDATE LIBRARY           #
#########################################

print "Checking wildtype sequence not in library\n";

if ( defined $wildtype_seq ) {
    if ( exists $library{ $wildtype_seq } ) {
        die "Library oligo matches wildtype sequence. Please remove wildtype sequence from library.\n$wildtype_seq\n";
    }
}

#########################################
#           PARSE READS                 #
#########################################
my %unmapped_seqs;
my $reads_mapped_to_wildtype = 0;

print "Mapping sample: $sample_fastq\n";

open my $gz, "gunzip -c $sample_fastq |" or die $!;
my $fq_in = Bio::SeqIO->new( -fh=> $gz, -format=>'fastq' );
my $reads_mapped_to_oligos = 0;
my $total_reads = 0;
while ( my $seq_obj = $fq_in->next_seq ) {
	$total_reads++;

	if ( exists $library{ $seq_obj->seq } ) {
		$library{ $seq_obj->seq }{'count'}++;
		$reads_mapped_to_oligos++;
	} elsif ( defined $wildtype_seq && $wildtype_seq eq $seq_obj->seq ) {
        $reads_mapped_to_wildtype++;
  } else {
		if ( exists $unmapped_seqs{ $seq_obj->seq } ) {
        $unmapped_seqs{ $seq_obj->seq }++
    } else {
        $unmapped_seqs{ $seq_obj->seq } = 1;
    }
	}
}
close ( $gz );

#########################################
#         WRITE UNMAPPED READS          #
#########################################
print "Writing unmapped reads: $unmapped_file\n";

my $total_unmapped_reads = 0;
open ( my $UNMAPPED, '>', $unmapped_file ) or die "Cannot open unmapped file ($unmapped_file): $!\n";
if ( scalar( keys %unmapped_seqs ) > 0 ) {
    foreach my $unmapped_seq ( keys %unmapped_seqs ) {
        $total_unmapped_reads += $unmapped_seqs{ $unmapped_seq };
        print $UNMAPPED join( "\t", $unmapped_seq, $unmapped_seqs{ $unmapped_seq } ) . "\n";
    }
}
close( $UNMAPPED );

#########################################
#         WRITE AMBIGUOUS READS         #
#########################################
print "Writing ambiguous reads: $ambiguous_file\n";

my $total_ambiguous_oligo_ids = 0;
my $total_ambiguous_oligo_seqs = 0;
my $total_ambiguous_reads = 0;

open ( my $AMBIGUOUS, '>', $ambiguous_file ) or die "Cannot open ambiguous file ($ambiguous_file): $!\n";
foreach my $seq ( keys %library ) {
    if ( scalar( @{ $library{ $seq }{'id'} } ) > 1 ) {
        my $ids = join( ",", @{ $library{ $seq }{'id'} } );
        print $AMBIGUOUS join( "\t", $seq, $ids, $library{ $seq }{'count'} ) . "\n";
        $total_ambiguous_oligo_seqs++;
        $total_ambiguous_oligo_ids += scalar( @{ $library{ $seq }{'id'} } );
        $total_ambiguous_reads += $library{ $seq }{'count'};
    }
}
close( $AMBIGUOUS );

#########################################
#              WRITE COUNTS             #
#########################################
print "Writing counts: $count_file\n";

# These will include ambiguous read counts
my $zero_counts = 0;
open ( my $COUNTS, '>', $count_file ) or die "Cannot open count file ($count_file): $!\n";
print $COUNTS join( "\t", 'id', $sample_name ) . "\n";
foreach my $seq ( keys %library ) {
	foreach my $id ( @{ $library{ $seq }{'id'} } ) {
		print $COUNTS join( "\t", $id, $library{ $seq }{'count'} ) . "\n";
	}
	if ( $library{ $seq }{'count'} == 0 ) {
		$zero_counts++;
	}
}
close( $count_file );

#########################################
#             WRITE STATISTICS          #
#########################################
print "Writing statistics: $stats_file\n";

my $total_library_seqs = scalar( keys %library );
my $unmapped_reads =  $total_unmapped_reads;

my $pct_reads_mapped_to_oligos = sprintf("%.2f", ( $reads_mapped_to_oligos / $total_reads ) * 100 ) . "%";
my $pct_reads_mapped_to_wildtype = sprintf("%.2f", ( $reads_mapped_to_wildtype / $total_reads ) * 100 ) . "%";
my $pct_unmapped_reads = sprintf("%.2f", ( $unmapped_reads / $total_reads ) * 100 ) . "%";
my $pct_ambiguous_reads = sprintf("%.2f", ( $total_ambiguous_reads / $total_reads ) * 100 ) . "%";

my $library_mapped = $total_library_seqs - $zero_counts;
my $pct_library_mapped = sprintf("%.2f", ( ($total_library_seqs - $zero_counts) / $total_library_seqs ) * 100 );
my $pct_zero_counts = sprintf("%.2f", ( $zero_counts / $total_library_seqs ) * 100 );

open ( my $STATS, '>', $stats_file ) or die "Cannot open count file ($stats_file): $!\n";

print $STATS join( "\t", 'sample',
                        'total_reads',
                        'library_oligos',
                        'library_mapped',
                        'unmapped_reads',
                        'ambiguous_reads',
                        'reads_mapped_to_wildtype',
                        'reads_mapped_to_oligos',
                        'ambiguous_library_ids',
                        'ambiguous_library_sequences',
                        'pct_library_mapped',
                        'pct_unmapped_reads',
                        'pct_ambiguous_reads',
                        'pct_reads_mapped_to_wildtype',
                        'pct_reads_mapped_to_oligos' ) . "\n";

print $STATS join( "\t", $sample_name,
                        $total_reads,
                        $total_library_seqs,
                        $library_mapped,
                        $unmapped_reads,
                        $total_ambiguous_reads,
                        $reads_mapped_to_wildtype,
                        $reads_mapped_to_oligos,
                        $total_ambiguous_oligo_ids,
                        $total_ambiguous_oligo_seqs,
                        $pct_library_mapped,
                        $pct_unmapped_reads,
                        $pct_ambiguous_reads,
                        $pct_reads_mapped_to_wildtype,
                        $pct_reads_mapped_to_oligos );
close ( $STATS );

print "Counts written to: $count_file\n";
print "Percentage library mapped: $pct_library_mapped ($library_mapped : $total_library_seqs)\n";
print "Percentage mapped reads (wildtype): $pct_reads_mapped_to_wildtype ($reads_mapped_to_wildtype : $total_reads)\n";
print "Percentage mapped reads (library): $pct_reads_mapped_to_oligos ($reads_mapped_to_oligos : $total_reads)\n";
