#!/usr/bin/perl -w
################################################################################
# Copyright or © or Copr. Anthony Mathelier and Alessandra Carbone (02-05-2010) 
# anthony.mathelier@gmail.com, alessandra.carbone@lip6.fr
# 
# This software is a computer program whose purpose is to provide a platform to
# check several questions around miRNAs and pre-miRNAs prediction and play
# between sensitivity and specificity by using parameters variability.
# 
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use, 
# modify and/or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info". 
# 
# As a counterpart to the access to the source code and rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author, the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability. 
# 
# In this respect, the user's attention is drawn to the risks associated
# with loading, using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean that it is complicated to manipulate, and that  also
# therefore means that it is reserved for developers and experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or 
# data to be ensured and,  more generally, to use and operate it in the 
# same conditions as regards security. 
# 
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
################################################################################

################################################################################
# Construct the file corresponding to the signature of potential precursors by
# using the signatures coming from miRDeep script and by adding the information
# of putative mature sequences.
################################################################################

use strict;
use warnings;
use Getopt::Std;

################################################################################
# Main
################################################################################

my $usage = 
"
Usage:
$0 -s <miRDeep signatures> -p <precursors> -o <sig output> -a <prec output>

";
my %inputs = ();
getopts("s:p:o:a:", \%inputs);

my $sig_file;
if($inputs{s}){
  $sig_file = $inputs{s};
}else{
  die $usage;
}
my $prec_file;
if($inputs{p}){
  $prec_file = $inputs{p};
}else{
  die $usage;
}
my $sig_output;
if($inputs{o}){
  $sig_output = $inputs{o};
}else{
  die $usage;
}
my $prec_output;
if($inputs{a}){
  $prec_output = $inputs{a};
}else{
  die $usage;
}

## Parsing of the miRDeep signature file

open(STREAM, "<$sig_file") or die "Cannot open $sig_file\n";
my %subjectHash;
while(my $line = <STREAM>){
  if($line =~ /^(\S+)\s+(\S+)\s+(\d+)\.+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\.+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(.+)$/){
    if(exists($subjectHash{$5}{"$1$7$8"})){
      die "Sig error with: $5 $1$7$8\n";
    }
    $subjectHash{$5}{"$1$7$8"} = $line;
  }
}
close(STREAM);

## Parsing of the file containing potential precursors

open(STREAM, "<$prec_file") or die "Cannot open $prec_file\n";
my %precHash;
my $prec;
my $query;
my $begin;
my $end;
while(my $line = <STREAM>){
  if($line =~ /^>(.+_x\d+)_(.+)_\d+:(\d+)_(\d+)_\d+ before:\d+ after:\d+$/){
    $prec = $2;
    $query = $1;
    $begin = $4 + 1;
    $end = $4 + $3;
  }elsif(not $line =~ /[\(\)]/){
    if(exists($precHash{$prec}{"$query$begin$end"})){
      die "Prec error with: $prec $query$begin$end\n";
    }
    $precHash{$prec}{"$query$begin$end"} = $line;
  }else{
    if(not exists($precHash{$prec}{"$query$begin$end"})){
      die "Struct without seq: $prec $query$begin$end\n";
    }
    $precHash{$prec}{"$query$begin$end"} = "$precHash{$prec}{\"$query$begin$end\"}$line";
  }
}
close(STREAM);

## Print the file corresponding to the signatures of the potential precursors,
## taking into account the putative miRNA mature. Notice that a potential
## precursors may contain several putative miRNAs, hence we have to give a
## specific name to each precursor/miRNA pair.

open(SIGOUT, ">$sig_output") or die "Cannot open $sig_output\n";
open(PRECOUT, ">$prec_output") or die "Cannot open $prec_output\n";
foreach my $prec (keys(%precHash)){
  my $prec_name = $prec;
  my $nb_miRNA = scalar(keys(%{$precHash{$prec}}));
  if($nb_miRNA > 1){
    my $count = -1;
    foreach my $mir_name (keys(%{$precHash{$prec}})){
      $count++;
      $prec_name = $prec . "_$count";
      foreach my $query (keys(%{$subjectHash{$prec}})){
        my $sig = $subjectHash{$prec}{$query};
        $sig =~ /^(\S+)\s+(\S+)\s+(\d+)\.+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\.+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(.+)$/;
        if($query eq $mir_name){
          print SIGOUT "$1\t$2\t$3..$4\t$prec_name\t$6\t$7..$8\t$9\t$10\t$11\t$12\tmature\n";
          my $strand = find_strand($sig);
          print PRECOUT ">$prec_name strand:$strand excise_beg:$7 excise_end:$8\n";
          print PRECOUT "$precHash{$prec}{$mir_name}";
        }else{
          print SIGOUT "$1\t$2\t$3..$4\t$prec_name\t$6\t$7..$8\t$9\t$10\t$11\t$12\n";
        }
      }
    }
  }else{
    foreach my $query (keys(%{$subjectHash{$prec}})){
      my @queries = keys(%{$precHash{$prec}});
      my $sig = $subjectHash{$prec}{$query};
      $sig =~ /^(\S+)\s+(\S+)\s+(\d+)\.+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\.+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(.+)$/;
      if($query eq $queries[0]){
        print SIGOUT "$1\t$2\t$3..$4\t$prec_name\t$6\t$7..$8\t$9\t$10\t$11\t$12\tmature\n";
        my $strand = find_strand($sig);
        print PRECOUT ">$prec_name strand:$strand excise_beg:$7 excise_end:$8\n";
        print PRECOUT "$precHash{$prec}{$query}";
      }else{
        print SIGOUT "$1\t$2\t$3..$4\t$prec_name\t$6\t$7..$8\t$9\t$10\t$11\t$12\n";
      }
    }
  }
}

sub find_strand{

    #A subroutine to find the strand, parsing different blast formats
    #Sub coming from miRDeep script miRDeep.pl

    my($other)=@_;

    my $strand="+";

    if($other=~/-/){
	$strand="-";
    }

    if($other=~/minus/i){
	$strand="-";
    }
    return($strand);
}
