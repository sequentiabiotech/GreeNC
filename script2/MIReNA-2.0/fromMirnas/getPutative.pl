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
# permet d'extraire les possibles preMiRna a partir d'un fichier resultat de
# l'algorithme asml entre des miRNAs connus et un texte
################################################################################

sub openError{
  my $file = shift;
  return "\nCannot open '$file'\n\n";
}

my $usage = "\n$0 <nb_nt_before> <nb_nb_after> <input file> ";
$usage .= "<outputFileName>\n\n";
die $usage if (@ARGV != 4);
my $before;
$before = shift or $before = 0;
my $after;
$after = shift or $after = 0;
my $inputfilename;
$inputfilename = shift or $inputfilename = "";
my $outputfile;
$outputfile = shift or $outputfile = "";

open (INPUTFILE, $inputfilename) or die openError($inputfilename);
open (OUTPUTFILE, ">", $outputfile) or die openError($outputfile);

my $textfile;
my $descriptionLine;
while (<INPUTFILE>){
  if (/on/){
    $_ =~ s/on //;
    chomp;
    $_ =~ /(.*) (.*)/;
    my $name =$1;
    $name =~ tr/ /_/;
    $textfile = $2;
    $descriptionLine = "\n>$name"."_$2\n"
  }elsif (/Match/){
    print OUTPUTFILE "$descriptionLine";
    $_ =~ /Match at (\d+) with length (\d+)/;
    my $pos = $1; # la position est donnée en numérotant à partir de 1 la
    # séquence du texte
    my $len = $2;
    open (TEXT, $textfile) or die openError($textfile);
    my $texte = <TEXT>;
    chomp($texte);
    my $deb = $pos - $len - $before;
    $deb = 0 if ($deb < 0);
    $bef = $before;
    if ($deb == 0){
      $bef = $pos - $len;
    }
    my $substring = substr ($texte, $deb, $len + $bef + $after);
    my $aft = $after;
    if (length ($substring) != $len + $bef + $after){
#      if ($pos - $len - $before < 0){
#        $bef = $pos - $len;
#      }
      if ($pos + $after >= length ($texte)){
        $aft = length ($texte) - $pos;
      }
    }
    printf OUTPUTFILE "%d:$len before:$bef after:$aft\n", $pos - $len + 1; # la position de
    # début du miRNA dans le texte est donné à partir de l'indice 1 pour la
    # 1ère lettre du texte
    print OUTPUTFILE $substring;
    print OUTPUTFILE ();
    print OUTPUTFILE "\n";
    close TEXT;
  }
}
close OUTPUTFILE;
close INPUTFILE;
