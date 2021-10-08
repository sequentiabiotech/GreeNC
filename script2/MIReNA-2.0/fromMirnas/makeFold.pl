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
#
################################################################################
# A partir d'un fichier de séquences avec la position des miRNAs putatifs
# associés, donne la liste des pre-miRNAs possibles
################################################################################
use strict;

my $makefold = $0;
$makefold =~ s/makeFold\.pl/\.\.\/folding\/makefold/;


my $usage = "\n$0 <-1|-2> <input file> <max amfe> <max mfei> <max unpaired> ";
$usage .= "<ratiomin> <ratiomax> <stems> <min len> <output file>\n-1 from getPutative.pl\n-2 otherwise\n\n";
if(@ARGV != 10){
  die $usage;
}

my $option = shift or die $usage;
my $input;
$input = shift or $input = "";
open(INPUTFILE, $input) or die "\nCannot open $input\n\n";
my $amfe = shift;# or die $usage;
my $mfei = shift;# or die $usage;
my $max_unpaired = shift;# or die $usage;
my $ratiomin = shift;# or die $usage;
my $ratiomax = shift;# or die $usage;
my $stem = shift;# or die $usage;
my $minlen = shift;# or die $usage;
my $output = shift or die $usage;
  
if($option == -1){
  my $name = "";
  my $nom = "";
  my $nt_before = 0;
  my $nt_after = 0;
  while (<INPUTFILE>) {
    if (/>(.*)_(.*)/) {
      $name = "$1_$2_";
    }elsif (/(\d+:\d+) before:(\d+) after:(\d+)/){
      $nom = $name.$1;
      $nt_before = $2;
      $nt_after = $3;
    }elsif ((/^[ACGTU]+$/) || (/^[acgtu]+$/)){
          chomp;
          system "$makefold $_ $nt_before $nt_after $amfe $mfei $max_unpaired $ratiomin $ratiomax $output $stem $minlen \"$nom\"";
      }
  }
  close INPUTFILE;
}elsif($option == -2){
  my $name = "";
  my $nom = "";
  my $nt_before = 0;
  my $nt_after = 0;
  while (<INPUTFILE>) {
    if (/>(.*)\s+before:(\d+)\s+after:(\d+)/) {
      $name = $1;
      $nom = $name;
      $nt_before = $2;
      $nt_after = $3;
    }elsif ((/^[ACGTU]+$/) || (/^[acgtu]+$/)){
      chomp;
      system "$makefold $_ $nt_before $nt_after $amfe $mfei $max_unpaired $ratiomin $ratiomax $output $stem $minlen \"$nom\"";
    }
  }
  close INPUTFILE;
}
