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
################################################################################use strict;


use Getopt::Std;

#Usage: countage (counttab, minen_prop, amfei, unpaired_block, nb_steam, mirlen,
#hairlen, mirstruct, hairstruct, nb nt before)
sub countage{
  my $countref = shift (@_);
  my $amfe = shift (@_);
  my $minAMFE = shift (@_);
  my $mfei = shift (@_);
  my $minMFEI = shift (@_);
  my $unpairedblock = shift (@_);
  my $nbsteam = shift (@_);
  my $mirlen = shift (@_);
  my $hairlen = shift (@_);
  my $mirstruct = shift (@_);
  my $hairstruct = shift (@_);
  my $bef = shift (@_);
  my $app = shift (@_);
  my $maxNbSteam = shift (@_);
  my $ratiomin = shift(@_);
  my $ratiomax = shift(@_);
  my $bool = 0;
  if ($mfei >= $minMFEI){
    @$countref[0]++;
    $bool = 1;
  }
  my $nonpaired = $unpairedblock / $mirlen;
  if ($unpairedblock / $mirlen > $app){
    @$countref[1]++;
    $bool = 1;
  }
  if ($nbsteam > $maxNbSteam){
    @$countref[2]++;
    $bool = 1;
  }
  if ($amfe >= $minAMFE){
    @$countref[5]++;
    $bool = 1;
  }

#  if ($hairlen > 225){
#    @$countref[3]++;
#    $bool = 1;
#  }

  if (($mirstruct =~ /\(/ && $mirstruct =~ /\)/) || $mirstruct =~ /^\.*$/) {
    @$countref[4]++;
    $bool = 1;
  }else{
    my $lencomplement = mircomplement ($hairstruct, $mirstruct, $bef);
    my $rapport = $mirlen / $lencomplement;
#    printf("rapport: %f\n", $rapport);
#    print "nonPaired: $nonpaired\n";
#    printf("absolute diff from 1 rapport: %f\n", abs(1-$rapport));
#    print "$mfei\n";
#    print "$amfe\n";
    if ($rapport > $ratiomax || $rapport < $ratiomin){
      @$countref[6]++;
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
#  print "hairpin = $hairstruct\n";
#  print "mirna = $mirstruct\n";
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
#  print "mirbegin = $mirbegin, mirend = $mirend\n";
#  print "compbegin = $appariement{$mirbegin}, compend = $appariement{$mirend}\n";
  $returnlen = abs ($appariement{$mirbegin} - $appariement{$mirend}) + 1;
#  print "complement length: $returnlen\n";
  return $returnlen;
}

################################################################################
sub withMiRnaPos{
  my $input = shift (@_);
  my $minAMFE = shift (@_);
  my $minMFEI = shift (@_);
  my $app = shift (@_);
  my $maxNbSteam = shift (@_);
  my $ratiomin = shift(@_);
  my $ratiomax = shift(@_);
  open INPUT, $input or die "Cannot open $input\n";
  my $nbCG;
  my $line;
  my $count  = 0;
  my @counts = (0,0,0,0,0,0,0);
  my $pl     = 0;
  my $name;
  my $nom;
  my $seq;
  my $struct;
  while ($line = <INPUT>){
    if ($line =~ />(.*?) before:(\d+) after:(\d+)/){
#    if ($line =~ />(.*)_(\d+)_(\d+)/){
#      print $line;
      chomp ($line);
      $nom = $line;
      $count++;
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
      $name =~ /.*_\d+:(\d+)_(\d+)_(\d+)/;
      if($length != $1 + $2 + $3){
        next;
      }
      $line =~ s/ \(( ?-\d+\.\d+)\)//;
#      my $length = length ($line);
      my $nrj = $1;
      my $l      = $length - $aft - $bef;
      my $mirna  = substr $line, $bef, $l;
#      print "line= $line, length= $length\n";
#      print "mirn= $mirna, length= $l\n";
      my $propCG = ( $nbCG / $length ) * 100;
      my $AMFE   = ( $nrj / $length ) * 100;
      my $MFEI;
      if ($propCG == 0){
        $MFEI = -1000.;
      }else{
        $MFEI = $AMFE / $propCG;
      }
#      print "amfe = $AMFE, MFEI= $MFEI\n";
      my $unpairedBlock = 0;
      $unpairedBlock++ while ($mirna =~ /\./g);
      $pl += $unpairedBlock / $l;
      my $nbSteam = 0;
      my $hairpin = $line;
      $line =~ s/\.//g;
      $nbSteam++ while ($line =~ /(\)+)/g);
      my $bool = countage (\@counts, $AMFE, $minAMFE, $MFEI, $minMFEI, $unpairedBlock,
        $nbSteam, $l, $length, $mirna, $hairpin, $bef, $app, $maxNbSteam,
        $ratiomin, $ratiomax);
      if (! $bool){
        print "$nom\n$seq\n$struct\n";
      }
    }
  }
  close INPUT;

  my $plaverage;
  if ($count != 0){
    $plaverage = $pl / $count;
  }else{
    $plaverage = -1;
  }
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
#      printf ("appUnpaired: %d\n", $nbUnpaired / length($mirna));
      if ($nbUnpaired / length ($mirna) > $app){
        next;
      }
      my $lencomplement = mircomplement ($hairpin, $mirna, $i);
      my $rapport = length ($mirna) / $lencomplement;
#      print "rapport: $rapport\n";
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
  my $maxNbSteam = shift (@_);
  my $ratiomin = shift(@_);
  my $ratiomax = shift(@_);
  open INPUT, $input or die "Cannot open $input\n";
  my $nbCG;
  my $line;
  my $count  = 0;
  my $nbOK  = 0;
  my $name;
  while ($line = <INPUT>){
    if ($line =~ />(.*)/){
#      print $line;
      $count++;
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
#      my $length = length($line);
      my $propCG = ( $nbCG / $length ) * 100;
      my $AMFE   = ( $nrj / $length ) * 100;
#      print "AMFE: $AMFE\n";
      my $MFEI = 0;
      if ($propCG == 0){
        $MFEI = -1;
      }else{
        $MFEI = $AMFE / $propCG;
      }
#      print "MFEI: $MFEI\n";
      my $minen_prop = $nrj / $length;
      my $nbSteam = 0;
      my $hairpin = $line;
      $line =~ s/\.//g;
      $nbSteam++ while ($line =~ /(\)+)/g);
#      print "nbSteam: $nbSteam\n";
      if ($MFEI < $mfei && $nbSteam <= $maxNbSteam# && $length <= 225 
        && $AMFE < $amfe){
        my $test = testMiRna ($hairpin, $app, $ratiomin, $ratiomax);
        if ($test){
          print ">$name\n";
          print $seq;
          print "$secstruct\n";
          $nbOK++;
        }
      }
    }
  }
#  print "nb total: $count, nbOK: $nbOK\n";
}

################################################################################
# MAIN
################################################################################

my $usage = "Use: $0 [-e|-p] -i <input file> -s <% non-appariés (<1) > ";
$usage .= "-r <ratiomin> -R <ratiomax> -a <amfe> -m <mfei> -n <nb steam max>\
#-p if we have miRna positions\n\n";
my %options=();
getopts("i:a:m:s:r:R:n:ep",\%options);
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
my $maxNbSteam = $options{n};
if($maxNbSteam == 0){
  $maxNbSteam = 1000000;
}
if (defined $options{p}){ # we know the miRna positions
  withMiRnaPos ($input, $amfe, $mfei, $app, $maxNbSteam, $ratiomin, $ratiomax);
}elsif (defined $options{e}){ # we want to extract the putative miRna
  withoutMiRnaPos ($input, $amfe, $mfei, $app, $maxNbSteam, $ratiomin, $ratiomax);
}else{
  die $usage;
}
