#!/usr/bin/env perl
use strict;
use warnings;

use File::Basename;
use File::Path;
use Getopt::Std;
use Scalar::Util qw(looks_like_number);
use Carp ();			# To track uninitialized variables using the local var below
local $SIG{__WARN__} = \&Carp::cluck;

# default input
my $fasta_file = 'input.fasta';
my $slip = 5;
my $temp = 37;
my $verbose = 0;

#check to see if all parameters are provided and read them
my $prog = basename($0);

getopt('hvi:d:t:'); 
use vars qw( $opt_h $opt_i $opt_d $opt_t $opt_v);

if( $opt_h ) 	{ usage( $prog ); exit(); } 
if( $opt_i ) 	{ $fasta_file = $opt_i; }
else 			{ usage( $prog ); exit(); }

if( $opt_d ) 	{ $slip = $opt_d; }
if( $opt_t ) 	{ $temp = $opt_t; }
if( $opt_v ) 	{ $verbose = 1; }

# check that input parameters are in reasonable ranges
die "ERROR -t <temperature> [ $temp ] must be non-negative number not exceeding 200." 
    unless 0 <= $temp && $temp <= 200;

# read the fasta file
my ($sequence, $length) = read_fasta($fasta_file); 

# check that slip makes sense
die "ERROR -d <distance threshold> [ $slip ] must not exceed the sequence length [ $length ]."
    unless $slip <= $length;

my $base = $fasta_file; 
$base =~ s/^.+\///; # remove dir path

#enforce base pairs
#one pair:
#cat ${base}Cent${iter}.37.plot |tail -n +2 |awk -v max=0 '{if($3>max){i=$1;j=$2;max=$3}}END{print "F",i,j,"1"}' > ${base}MFE${iter}.aux
#percentile of pairs:
#cat ${base}Cent${iter}.37.plot |tail -n +2 |awk '{if ($3>0.8) print "F",$1,$2,"1"}' > ${base}MFE${iter}.aux
#stable stem:
#run_command("java -cp . selectcentroidbp.Main " . $fasta_file . " " . $base . "Cent" . $iter . ".37.plot 0 > " . $base . "MFE" . $iter. ".aux");

run_command("hybrid-ss-min -t 37 -T 37 " . $fasta_file . " -o " . $base . "MFE0");

# remove $output_file if it exists
my $output_file = $base . ".out";
if(-e $output_file) { unlink $output_file; }

my $energy = get_energy($temp,$base . "MFE0.ct");

open(OUT, ">>", $output_file) || die "ERROR unable to open $output_file for writing"; 
print OUT (">" . $energy . "\n");
print OUT ($sequence . "\n");
close(OUT);
sleep(0.1); 

run_command("java -cp . unafold2vienna.Main " . $fasta_file . " " . $base . "MFE0.37.plot >> $output_file");

#cat ${base}Cent${prev}.37.plot |tail -n +2 |awk '{if ($3>0.5) print "P",$1,$2,"1"}' >> ${base}Cent${iter}.aux
system("java -cp . slip.Main " . $base . "MFE0.37.plot $slip $length >> " . $base . "Cent1.aux");
#cat ${base}MFE${prev}.37.plot |tail -n +2 |awk '{print "P",$1,$2,$3}' >> ${base}Cent${iter}.aux
system("hybrid-ss -t " . $temp . " -T " . $temp . " " . $fasta_file  . " -c -o " . $base . "Cent1");

#enforce base pair
#one pair:
#cat ${base}Cent${iter}.37.plot |tail -n +2 |awk -v max=0 '{if($3>max){i=$1;j=$2;max=$3}}END{print "F",i,j,"1"}' > ${base}MFE${iter}.aux
#percentile of pairs:
#cat ${base}Cent${iter}.37.plot |tail -n +2 |awk '{if ($3>0.8) print "F",$1,$2,"1"}' > ${base}MFE${iter}.aux
#stable stem:
system("java -cp . selectcentroidbp.Main " . $fasta_file . " " . $base . "Cent1." . $temp . ".plot 0 > " . $base . "MFE1.aux");
system("hybrid-ss-min -t " . $temp . " -T " . $temp . " " . $fasta_file . " -c -o " . $base . "MFE1");

$energy = get_energy($temp,$base . "MFE1.ct");

open(OUT, ">>", $output_file) || die "ERROR unable to open $output_file for writing"; 
print OUT (">" . $energy . "\n");
print OUT ($sequence . "\n");
close(OUT);
sleep(0.1); 

run_command("java -cp . unafold2vienna.Main " . $fasta_file . " " . $base . "MFE1." . $temp . ".plot >> $output_file");

if (! $verbose) { clean($fasta_file); }

#################################################################
# SUBROUTINES 
#################################################################

sub clean
{
  my $fasta_file = shift;

  opendir( PWD, '.' );
  my @fns = grep( /^$fasta_file/, readdir( PWD ) );

  foreach my $fn (@fns) {
  
    if ($fn =~ /^($fasta_file)(Cent1|MFE0|MFE1)(\.37\.(ext|plot))/) { unlink $fn; }
    if ($fn =~ /^($fasta_file)(Cent1)\.(aux|dG|run)/) { unlink $fn; }
	if ($fn =~ /^$fasta_file(MFE0|MFE1)\.(ct|dG|run)/) { unlink $fn; }
  }

  closedir( PWD );
}

sub run_command {
  my $sub_name = "run_command()";
  my $nargs_expected = 1;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($cmd) = @_;
  
  system($cmd);

  if($? != 0) { 
    die "ERROR, the following command failed:\n$cmd\n", $?;
  }

  return;
}

sub read_fasta { 
  my $sub_name = "read_fasta()";
  my $nargs_expected = 1;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($fasta_file) = @_;
  my $line;
  
  open(FASTA, $fasta_file) || die "ERROR unable to open $fasta_file"; 
  my $header = <FASTA>;
  if($header !~ /^\>/) { 
    die "ERROR first line of fasta file $fasta_file does not begin with \">\":\n$line\n";
  }

  my $seq = "";
  my $seqlen = 0;
  my $keep_going = 1;  # set to 0 if we find another header line
  my $seqline; 
  while((defined ($seqline = <FASTA>)) && ($keep_going)) { 
    chomp $seqline;
    if($seqline =~ m/^\>/) { 
      $keep_going = 0;
    }
    else { 
	  $seqline =~ s/\s//g;
      my @char_A = split("", $seqline);
      foreach my $char (@char_A) { 
        if($char !~ m/[ACGUacgu]/) { 
          die "ERROR found non-RNA sequence character (A|C|G|U) in sequence: $char\n";
        }
      }
      $seq .= $seqline;
      $seqlen += length($seqline);
    }
  }
  close(FASTA);

  if($seqlen < 1) { 
    die "ERROR first sequence in fasta file has length below one: $seqlen\n";
  }
  if($seqlen > 500) { 
    die "ERROR first sequence in fasta file has length above 500: $seqlen\n";
  }

  return ($seq, $seqlen);
}

sub get_energy { 
  my $sub_name = "get_energy()";
  my $nargs_expected = 2;

  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($temp,$infile) = @_;

  $energy = `ct-energy -t $temp $infile`;
  chomp $energy;
  return $energy;
}

#########################################################################################

# Usage display

sub usage 
# program_name_
{
    my $program_name_ = shift;
    my $version = '2.1';
    
    my $p = $program_name_ . ' ' . $version;
    print STDERR <<USAGE
    
$program_name_ (Version $version)

Usage: perl $program_name_ [options]

This program finds the alternative RNA secondary structures based on conditional base-pair probabilities.

The program inputs: 
    (1) an RNA sequence in fasta format in [input_file];
    (2) the distance threshold for base-pair exclusion; and
    (3) temperature in degrees Celsius for folding the alternative structure.

The output files are: 
    (1) [input_file].out
            contains dominant and alternative structure predictions, 
			with energy, sequence, and structure; and
    (2) [input_file]MFE1.aux
            contains the base-pairs of the stem seed L*.

Program options are arguments with \'-\' followed by a letter.
An option requiring further input(s) appears with a colon.
The default for an option, if any, is indicated in parentheses.

-i : input fasta file
-d : (5) tau = the distance threshold for excluded base-pairs from the MFE structure
-t : (37) T = temperature in degrees Celsius for folding the alternative structure
USAGE
}
#########################################################################################

1