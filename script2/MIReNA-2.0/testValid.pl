#!/usr/bin/perl -w
#################################################################################
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
use strict;
use Getopt::Std;

sub testValidity{
  my $amfe = shift (@_);
  my $minAMFE = shift (@_);
  my $mfei = shift (@_);
  my $minMFEI = shift (@_);
  my $unpairedblock = shift (@_);
  my $nbstem = shift (@_);
  my $mirlen = shift (@_);
  my $hairlen = shift (@_);
  my $mirstruct = shift (@_);
  my $hairstruct = shift (@_);
  my $bef = shift (@_);
  my $app = shift (@_);
  my $maxNbStem = shift (@_);
  my $ratiomin = shift(@_);
  my $ratiomax = shift(@_);
  my $bool = 0;
  if ($mfei >= $minMFEI){
    $bool = 1;
  }
  my $nonpaired = $unpairedblock / $mirlen;
  if ($unpairedblock / $mirlen > $app){
    $bool = 1;
  }
  if ($maxNbStem && $nbstem > $maxNbStem){
    $bool = 1;
  }
  if ($amfe >= $minAMFE){
    $bool = 1;
  }

  if (($mirstruct =~ /\(/ && $mirstruct =~ /\)/) || $mirstruct =~ /^\.*$/) {
    $bool = 1;
  }else{
    my $lencomplement = mircomplement($hairstruct, $mirstruct, $bef);
    my $rapport = $mirlen / $lencomplement;
    if ($rapport > $ratiomax || $rapport < $ratiomin){
      $bool = 1;
    }
  }
  return $bool;
}

################################################################################
sub mircomplement{ #pré-condition: le miRna ne se replie pas sur lui-même
  my $hairstruct = shift (@_); #structure secondaire de l'hairpin
  my $mirstruct = shift (@_); #structure secondaire du miRna
  my $mirbegin = shift (@_); #position de départ du miRna dans l'hairpin
  my $bool = 0;
  if ($mirstruct =~ /\)/){
    $bool = 1;
  }
  my %appariement;
  my @pile = ();
  for (my $i = 0; $i < length ($hairstruct); $i++){
    my $char = substr ($hairstruct, $i, 1);
    if ($char eq '('){
      unshift (@pile, $i);
    }elsif ($char eq ')'){
      if ($bool){
        $appariement{$i} = shift (@pile);
      }else{
        $appariement{shift (@pile)} = $i;
      }
    }
  }
  my $returnlen;
  my $mirend = $mirbegin + length ($mirstruct) - 1;
  while (substr ($hairstruct, $mirbegin, 1) eq '.'){
    $mirbegin++;
  }
  while (substr ($hairstruct, $mirend, 1) eq '.'){
    $mirend--;
  }
  $returnlen = abs ($appariement{$mirbegin} - $appariement{$mirend}) + 1;
  return $returnlen;
}

################################################################################
sub withMiRnaPos{
  my $input = shift (@_);
  my $minAMFE = shift (@_);
  my $minMFEI = shift (@_);
  my $app = shift (@_);
  my $maxNbStem = shift (@_);
  my $ratiomin = shift(@_);
  my $ratiomax = shift(@_);
  open INPUT, $input or die "Cannot open $input\n";
  my $nbCG;
  my $line;
  my $pl     = 0;
  my $name;
  my $nom;
  my $seq;
  my $struct;
  while ($line = <INPUT>){
    if ($line =~ />(.*?) before:(\d+) after:(\d+)/){
      chomp ($line);
      $nom = $line;
      $name = $1;
      my $bef = $2;
      my $aft = $3;
      $line = <INPUT>; #lecture de la séquence
      chomp ($line);
      $seq = $line;
      my $length = length($line);
      $nbCG = 0;
      $nbCG++ while ( $line =~ /[CG]/g );
      $line = <INPUT>; #lecture de la structure secondaire
      chomp $line;
      $struct = $line;
      $line =~ s/ \(( ?-\d+\.\d+)\)//;
      my $nrj = $1;
      my $l      = $length - $aft - $bef;
      my $mirna  = substr $line, $bef, $l;
      my $propCG = ( $nbCG / $length ) * 100;
      my $AMFE   = ( $nrj / $length ) * 100;
      my $MFEI;
      if ($propCG == 0){
        $MFEI = -1000.;
      }else{
        $MFEI = $AMFE / $propCG;
      }
      my $unpairedBlock = 0;
      $unpairedBlock++ while ($mirna =~ /\./g);
      $pl += $unpairedBlock / $l;
      my $nbStem = 0;
      my $hairpin = $line;
      $line =~ s/\.//g;
      $nbStem++ while ($line =~ /(\)+)/g);
      my $bool = testValidity($AMFE, $minAMFE, $MFEI, $minMFEI, $unpairedBlock,
        $nbStem, $l, $length, $mirna, $hairpin, $bef, $app, $maxNbStem,
        $ratiomin, $ratiomax);
      if (! $bool){
        print "$nom\n$seq\n$struct\n";
      }
    }
  }
  close INPUT;
}

################################################################################
sub testMiRna{
  my $hairpin = shift (@_);
  my $app = shift (@_);
  my $ratiomin = shift(@_);
  my $ratiomax = shift(@_);
  my $length = length ($hairpin);
  my $bool = 0;
  for (my $i = 0; $i < $length - 18; $i++){
    for (my $len = 18; $len < 26 && $len + $i - 1 < $length; $len++){
      my $mirna = substr ($hairpin, $i, $len);
      if ($mirna =~ /\(/ && $mirna =~ /\)/){
        next;
      }
      my $countPoint = 0;
      $countPoint++ while ($mirna =~ /\./g);
      if ($countPoint == length ($mirna)){
        next;
      }
      my $nbUnpaired = 0;
      $nbUnpaired++ while ($mirna =~ /\./g);
      if ($nbUnpaired / length ($mirna) > $app){
        next;
      }
      my $lencomplement = mircomplement ($hairpin, $mirna, $i);
      my $rapport = length ($mirna) / $lencomplement;
      if ($rapport > $ratiomax || $rapport < $ratiomin){
        next;
      }
      $bool = 1;
    }
    if ($bool){
      last;
    }
  }
  return $bool;
}

################################################################################
sub withoutMiRnaPos{
  my $input = shift (@_);
  my $amfe = shift (@_);
  my $mfei = shift (@_);
  my $app = shift (@_);
  my $maxNbStem = shift (@_);
  my $ratiomin = shift(@_);
  my $ratiomax = shift(@_);
  my $minLen = shift(@_);
  open INPUT, $input or die "Cannot open $input\n";
  my $nbCG;
  my $line;
  my $name;
  while ($line = <INPUT>){
    if ($line =~ />(.*)/){
      $name = $1;
      $line = <INPUT>; #lecture de la séquence
      my $seq = $line;
      chomp($line);
      my $length = length($line);
      $nbCG = 0;
      $nbCG++ while ( $line =~ /[CG]/g );
      $line = <INPUT>; #lecture de la structure secondaire
      chomp $line;
      my $secstruct = $line;
      $line =~ s/ \(( ?-\d+\.\d+)\)//;
      my $nrj = $1;
      my $propCG = ( $nbCG / $length ) * 100;
      my $AMFE   = ( $nrj / $length ) * 100;
      my $MFEI = 0;
      if ($propCG == 0){
        $MFEI = -1;
      }else{
        $MFEI = $AMFE / $propCG;
      }
      my $minen_prop = $nrj / $length;
      my $nbStem = 0;
      my $hairpin = $line;
      $line =~ s/\.//g;
      $nbStem++ while ($line =~ /(\)+)/g);
      if ($MFEI < $mfei && ($nbStem <= $maxNbStem || $maxNbStem == 0)
        && $AMFE < $amfe && $length >= $minLen){
        my $test = testMiRna($hairpin, $app, $ratiomin, $ratiomax);
        if ($test){
          print ">$name\n";
          print $seq;
          print "$secstruct\n";
        }
      }
    }
  }
}

################################################################################
# MAIN
################################################################################

my $usage = "$0 [-p|-x] -i <input file> -s <unpaired nt (<1)> ";
$usage .= "-r <ratiomin> -R <ratiomax> -a <amfe> -m <mfei> ";
$usage .= "-n <max stems> -l <min length>\n\n";
my %options=();
getopts("i:a:m:s:r:R:n:xpl:",\%options);
die $usage if (!defined $options{i});
my $input = $options{i};
die $usage if (!defined $options{a});
my $amfe = $options{a};
$amfe =~ tr/,/\./;
die $usage if (!defined $options{m});
my $mfei = $options{m};
$mfei =~ tr/,/\./;
die $usage if (!defined $options{s});
my $app = $options{s};
die $usage if (!defined $options{r});
my $ratiomin = $options{r};
die $usage if (!defined $options{R});
my $ratiomax = $options{R};
die $usage if (!defined $options{n});
my $maxNbStem = $options{n};
if (defined $options{p}){ # we know the miRna positions
  withMiRnaPos ($input, $amfe, $mfei, $app, $maxNbStem, $ratiomin, $ratiomax);
}elsif (defined $options{x}){ # we want to extract the putative miRna
  die $usage if (!defined $options{l});
  my $minLen = $options{l};
  withoutMiRnaPos ($input, $amfe, $mfei, $app, $maxNbStem, $ratiomin, $ratiomax, $minLen);
}else{
  die $usage;
}
