#!/usr/bin/env perl

our @BOOL = qw(false true);
our $VERBOSE = 0;
our (%DEPTH, %GENES);

# usage: perl region_depth.pl coord file_depth > output
# coord = tab-delimited file, first column contig, second column region start, third column region end
# file_depth = output from $samtools depth -a file.bam > file_depth
# author: Danny W. Rice

use warnings;
use strict;
use Getopt::Long;

my ($provirus_file, @depth_files) = parseArgs();
getGenes($provirus_file);
getDepth(\@depth_files);
findGeneDepth();

sub findGeneDepth {
  my ($i, $orient, $gene, $s, $e, $contig, $depth, $name);
  foreach $contig ( sort keys %GENES ) {
    foreach $gene ( @{$GENES{$contig}} ) {
      ($s, $e, $orient, $name) = @{$gene};
      $depth = 0;
      for( $i = $s; $i <= $e; $i++ ) {
        if(defined $DEPTH{$contig}[$i] ) {
          $depth += $DEPTH{$contig}[$i];
        }
      }
      $depth /= ($e - $s + 1);
      if( $orient < 0 ) {
        printf "%.5f\t%s\n", $depth, join "\t", $contig, $e, $s;
      } else {
        printf "%.5f\t%s\n", $depth, join "\t", $contig, $s, $e;
      }
    }
  }
}

sub getGenes {
  my $file = shift;
  open my $fh, $file or die "error: can't open file=$file\n";
  my ($orient, $gene, $s, $e, $contig);
  while (<$fh>) {
    /(\S+)\t(\d+)\t(\d+)/ or next;
    if ($2 > $3) {
      $orient = -1;
      ($s, $e) = ($3, $2);
    } else {
      $orient = 1;
      ($s, $e) = ($2, $3);
    }
    $gene = $1;
    ($contig = $gene) =~ s/_\d+$//;
    push @{$GENES{$contig}}, [$s, $e, $orient, $gene];
  }
}

sub getDepth {
  my $files = shift;
  foreach my $file ( @$files ) {
    open my $fh, $file or die "error: can't open file=$file\n";
    my @F;
    while (<$fh>) {
      /^(\S+)\t(\d+)\t(\d+)/ or next;
      exists $GENES{$1} or next;
      $DEPTH{$1}[$2] = $3;
    }
  }
}

sub parseArgs {
  my @args = qw(provirus-file depth-file);
  my $usage = qq{usage: $0 @args [depth-file ...] [options]

options
-------
-verbose (default: $BOOL[$VERBOSE])
};

  my $result = GetOptions
    (
     'verbose!' => \$VERBOSE,
    );
  $result or print $usage and exit;
  @ARGV >= @args or print $usage and exit;
  return @ARGV;
}
