#!/usr/bin/env bash

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

usage="\n$0 <OPTIONS>
  OPTIONS:
    -f <fasta file> : file containing the miRNA sequence in fasta format
    -t <text file>  : file containing the DNA sequence
    -e <max error>  : maximum percentage of errors to match a miRNA sequence
    -l <max len>    : maximum length of a match
    -o <output>     : output file
    -h              : print this help\n"

if [ "$#" != 10 ]
then
  echo -e "$usage"
  exit 1
fi
while getopts "f:t:e:o:l:h" options
do
  case $options in
    f) file=$OPTARG;;
    t) text=$OPTARG;;
    e) errors=$OPTARG;;
    o) output=$OPTARG;;
    l) length=$OPTARG;;
    h) echo -e "$usage";;
    \?) echo -e "$usage"
        exit 1;;
    *) echo -e "$usage"
       exit 1;;
  esac
done
if [ ! -e $text ]
then
  echo -e "\n$text does not exist\n"
  exit 1
fi
if [ ! -e $file ]
then
  echo -e "\n$file does not exist\n"
fi
src=${0%%/executeAsml.sh}
TMP1=`mktemp`
TMP2=`mktemp`
awk '/>/{print $$0}' $file | tr -d '>' > $TMP1
awk '!/>/{print $$0}' $file > $TMP2

echo "$errors percent of errors and $length as maximal length" > $output
declare -i count=0
count=0
for j in `cat $TMP2`; do
  count=$count+1
  l=`head -n $count $TMP1 | tail -n 1`
  tmp=`mktemp`
  $src/../asml/asml.c $j $errors $text $length > $tmp
  declare -i nbl
  nbl=`wc -l $tmp | cut -f 1 -d ' '`
  if (($nbl != 0))
  then
    echo $l on $text >> $output
    tmp2=`mktemp`
    cat $output $tmp > $tmp2
    mv $tmp2 $output
  fi
  rm $tmp
done
rm $TMP1 $TMP2
