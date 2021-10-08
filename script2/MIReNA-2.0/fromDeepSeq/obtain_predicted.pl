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


################################################################################
# Parse the miRDeep-like output file to extract the predictions made by MIReNA
################################################################################

use strict;

################################################################################
# Main
################################################################################

my $usage = "$0 <input file>\n";
open STREAM, $ARGV[0] or die $usage;

my $mature_beg;
my $mature_end;
my $query;
my $id;
my $seq;
my $struct;
my %precursors_hash;
while(my $line = <STREAM>){
  if($line =~ /^mature_beg\s+(\d+)/){
    $mature_beg = $1;
  }elsif($line =~ /^mature_end\s+(\d+)/){
    $mature_end = $1;
  }elsif($line =~ /^mature_query\s+(\S+)/){
    $query = $1;
  }elsif($line =~ /^pri_id\s+(\S+)/){
    $id = $1;
    $id =~ /(.*_.*)_.*$/;
    $id = $1;
  }elsif($line =~ /^pri_seq\s+(\S+)/){
    $seq = $1;
  }elsif($line =~ /^pri_struct\s+(\S+)/){
    $struct = $1;
    if(exists $precursors_hash{$id}){
      my @tuple = @{$precursors_hash{$id}};
      my $q = $tuple[2];
      $q =~ /\S+_x(\d+)/;
      my $freq_q = $1;
      $query =~ /\S+_x(\d+)/;
      my $freq_query = $1;
      if($freq_query > $freq_q){
        $precursors_hash{$id} = [$mature_beg, $mature_end, $query, $id, $seq, $struct];
      }
    }else{
      $precursors_hash{$id} = [$mature_beg, $mature_end, $query, $id, $seq, $struct];
    }
  }
}
close(STREAM);

foreach $id (keys(%precursors_hash)){
  my @pred = @{$precursors_hash{$id}};
  print ">$pred[2] $pred[3] begin:$pred[0] end:$pred[1]\n$pred[4]\n$pred[5]\n";
}
