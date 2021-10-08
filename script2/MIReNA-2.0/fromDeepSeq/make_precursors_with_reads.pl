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
#################################################################################

use strict;
use warnings;
use Getopt::Std;

################################################################################
# Main
################################################################################
my $usage =
"
$0 -p <precursors file> -s <signatures file>

";

my %inputs = ();
getopts("s:p:", \%inputs);
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

## Parsing the precursors file
open(STREAM, "<$prec_file") or die "Cannot open $prec_file\n";
my %precHash;
while(my $line = <STREAM>){
  if($line =~ /^>(\S+)\s+.*$/){
    my $prec = $1;
    my $seq = <STREAM>;
    my $struct = <STREAM>;
    $precHash{$prec} = "$seq$struct";
  }
}
close(STREAM);

## Parsing the miRDeep signature file
open(STREAM, "<$sig_file") or die "Cannot open $sig_file\n";
while(my $line = <STREAM>){
  if($line =~ /^(\S+)\s+(\S+)\s+(\d+)\.+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\.+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(.+)$/){
    my $before = $7 - 1;
    my $after = $6 - $8;
    print ">$1"."_$5"."_$7:$2"."_$before"."_$after before:$before after:$after\n";
    print $precHash{$5};
  }
}
close(STREAM);
