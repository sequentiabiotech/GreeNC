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

################################################################################
# MIReNA implementation made by Anthony Mathelier
################################################################################

SUCCESS=0
FAILURE=-1

isdigit ()    # Tests whether *entire string* is numerical.
{             # In other words, tests for integer variable.
  [ $# -eq 1 ] || return $FAILURE

  case $1 in
    *[!0-9]*|"") return $FAILURE;;
              *) return $SUCCESS;;
  esac
}

isfloatordigit(){
  [ $# -eq 1 ] || return $FAILURE

  case $1 in
    *[!0-9]*|"")
      case $1 in
        [0-9]*.[0-9]*) return $SUCCESS;;
        -[0-9]*.[0-9]*) return $SUCCESS;;
        [0-9]*) return $SUCCESS;;
        -[0-9]*) return $SUCCESS;;
        *) return $FAILURE;;
      esac;;
    *) return $SUCCESS;;
  esac
}

function checkFile
{
  if [ ! -e $1 ]
  then
    echo -e "\n$1 does not exist\n"
    exit 1
  fi
}

function checkOutput
{
  if [ -e $1 ]
  then
    echo -e "\n$1 already exists\n"
    exit 1
  fi
}

function checkErrors
{
  if ! isfloatordigit $1
  then
    echo -e "\nYou must specify a valid errors percentage ($1 given)\n"
    exit 1
  fi

  if [[ $(echo "$1 < 0" | bc) -eq 1 || $(echo "$1 > 100" | bc) -eq 1 ]]
  then
    echo -e "\nYou must specify a valid errors percentage ($1 given)\n"
    exit 1
  fi
}

function checkLength
{
  if ! isdigit $1
  then
    echo -e "\nYou must specify a valid $2 value ($1 given)\n"
    exit 1
  fi
  if [[ $(echo "$1 < 0" | bc) -eq 1 ]]
  then
    echo -e "\nYou must specify a valid $2 value ($1 given)\n"
    exit 1
  fi
}

function checkNRJ
{
  if ! isfloatordigit $1
  then
    echo -e "\nYou must specify a valid $2 value ($1 given)"
    echo -e "$usage"
    exit 1
  fi
}

function checkVal
{
  if [[ "$2" == "non-app" || "$2" == "ratiomin" ]]
  then
    if [[ $(echo "$1 < 0" | bc) -eq 1 || $(echo "$1 > 1" | bc) -eq 1 ]]
    then
      echo -e "\nInvalid $2 value ($1 given)\n"
      exit 1
    fi
  else
    if [[ $(echo "$1 < 1" | bc) -eq 1 ]]
    then
      echo -e "\nInvalid $2 value ($1 given)\n"
      exit 1
    fi
  fi
}

function checkOption
{
  if [[ "$1" != "-p" && "$1" != "-x" ]]
  then
    echo -e "\nInvalid option: $1\n"
    exit 1
  fi
}

function checkMIR
{
  checkFile $1
  checkFile $2
  checkOutput $3
  checkErrors $4
  checkLength $5 "Length"
  checkNRJ $6 "MFEI"
  checkNRJ $7 "AMFE"
  checkVal $8 "non-app"
  checkVal $9 "ratiomin"
  checkVal ${10} "ratiomax"
  checkLength ${11} "stem"
  checkLength ${12} "left"
  checkLength ${13} "right"
  checkLength ${14} "length"
}

function checkDEEP
{
  checkFile $1
  checkFile $2
  checkOutput $3
  checkFile $4
  checkNRJ $5 "MFEI"
  checkNRJ $6 "AMFE"
  checkVal $7 "non-app"
  checkVal $8 "ratiomin"
  checkVal $9 "ratiomax"
  checkLength ${10} "stem"
  checkLength ${11} "nbreads"
}

function checkPredict
{
  checkFile $1
  checkOutput $2
  checkNRJ $3 "MFEI"
  checkNRJ $4 "AMFE"
  checkVal $5 "non-app"
  checkVal $6 "ratiomin"
  checkVal $7 "ratiomax"
  checkLength $8 "stem"
  checkLength $9 "length"
}

function checkValid
{
  checkFile $1
  checkOutput $2
  checkNRJ $3 "MFEI"
  checkNRJ $4 "AMFE"
  checkVal $5 "non-app"
  checkVal $6 "ratiomin"
  checkVal $7 "ratiomax"
  checkLength $8 "stem"
  checkLength $9 "length"
  checkOption ${10}
}

function checkRNAfold
{
  if ! which RNAfold > /dev/null
  then
    echo -e "\nYou must install RNAfold and put it executable directory in the PATH\n";
    exit 1;
  fi
}

usage="$0 <OPTIONS>

DO NOT FORGET TO READ THE README FILE

  OPTIONS:
**-----------------**
    --MIR  | -M : execute MIReNA algorithm starting from known miRNAs
        Required PARAMETERS:
        --file    | -f : a fasta file containing the miRNA sequences
                    (default: miRNAs.fa)
        --text    | -t : a text where to search for miRNA similar sequences and
                    pre-miRNAs (default: text.txt)
        --errors  | -e : maximum number of admited errors between miRNAs in
                    the fasta file and putative miRNAs in a text file;
                    the value is given in percentage (default: 15)
        --Length  | -L : maximum length of a putative (similar) miRNA
                    (default: 25)
        --left    | -g : left extension of putative miRNAs (default: 200)
        --right   | -d : right extension of putative miRNAs (default: 200)
        --length  | -l : minimum length for an accepting pre-miRNA (default: 60)
        --mfei    | -m : maximum MFEI (default: -0.85)
        --amfe    | -a : maximum AMFE (default: -32)
        --non-app | -n : maximum unmatched nucleotide threshold (<1)
                    (default: 0.26)
        --ratiomin| -r : minimum length ratio between miRNA and miRNA* (<1)
                    (default: 0.83)
        --ratiomax| -R : maximum length ratio between miRNA and miRNA* (>1)
                    (default: 1.17)
        --stem    | -s : maximum number of stem(s) (default: 0, i.e. unbound)
        --output  | -o : output file (default: out.txt)

**-----------------**
    --DEEP  | -D : execute MIReNA starting from deep sequencing reads
        Required PARAMETERS:
        --file    | -f : a fasta file containing the deep sequencing reads
                    (default: deep.fa)
        --genome  | -j : a genome (in fasta format) where to search for deep
                    sequencing reads (default: genome.fa)
        --mirnas  | -k : a fasta file containing already known miRNAs to be
                    checked for conservation (nt 2-8) (default: miRNAs.fa)
        --nbreads | -z : the minimum number of reads matching a pre-miRNA and
                    satisfying Dicer processing conditions. This number is used
                    in condition 3. discussed in the section \"Dicer processing 
                    condtion and mapping of reads to filter precursors\" of
                    Methods of the article (default: 3)
        --mfei    | -m : maximum MFEI (default: -0.62)
        --amfe    | -a : maximum AMFE (default: -26)
        --non-app | -n : maximum unmatched nucleotide threshold (<1)
                    (default: 0.26)
        --ratiomin| -r : minimum length ratio between miRNA and miRNA* (<1)
                    (default: 0.83)
        --ratiomax| -R : maximum length ratio between miRNA and miRNA* (>1)
                    (default: 1.17)
        --stem    | -s : maximum number of stems (default: 0, i.e. unbound)
        --output  | -o : output file (default: out.txt)
        
        Optional PARAMETER:
        --blast   | -b : a file that is the output of the script
                    miRDeep/blastoutparsed.pl. This file corresponds to a blast
                    search of deep sequencing reads on the genome parsed by the
                    miRDeep script. If this option is used, MIReNA will not redo
                    the blast search but use this file to create potential
                    precursors


**-----------------**
    --predict | -p : predict valid pre-miRNA sequences including putative
                    miRNAs occuring in the input file
        Required PARAMETERS:
        --file    | -f : an input file containing sequences, putative miRNAs and
                    their positions the corresponding sequences; look for
                    pre-miRNA sequences including the miRNAs (default: input.txt)
        --length  | -l : minimal length for an accepting pre-miRNA (default: 60)
        --mfei    | -m : maximum MFEI (default: -0.85)
        --amfe    | -a : maximum AMFE (default: -32)
        --non-app | -n : maximum unmatched nucleotide threshold (<1)
                    (default: 0.26)
        --ratiomin| -r : minimum length ratio between miRNA and miRNA* (<1)
                    (default: 0.83)
        --ratiomax| -R : maximum length ratio between miRNA and miRNA* (>1)
                    (default: 1.17)
        --stem    | -s : maximum number of stems (default: 0, i.e. unbound)
        --output  | -o : output file (default: out.txt)

**-----------------**
    --valid   | -v : print pre-miRNA sequences that validate criteria I-V;
                    sequences are contained in the fasta file
        Required PARAMETERS:
        --file    | -f : a fasta file containing the putative pre-miRNA sequences
                        (default: preMi.fa)
        --length  | -l : minimal length for an accepting pre-miRNA (default: 60)
        --mfei    | -m : maximum MFEI (default: -0.85)
        --amfe    | -a : maximum AMFE (default: -32)
        --non-app | -n : maximum unmatched nucleotide threshold (<1)
                    (default: 0.26)
        --ratiomin| -r : minimum length ratio between miRNA and miRNA* (<1)
                    (default: 0.87)
        --ratiomax| -R : maximum length ratio between miRNA and miRNA* (>1)
                    (default: 1.17)
        --stem    | -s : maximum number of stems (default: 0, i.e. unbound)
        --output  | -o : output file (default: out.txt)

        Parameter TO CHOOSE:
        [-y | -x]      : if -y option is given, the input file contains
                    information about the miRNA positions in putative
                    pre-miRNAs, otherwise (with -x option) no specific miRNAs
                    are considered

**-----------------**
    --help    | -h : print this help\n"

    TEMP=`getopt -o MDpvhk:j:b:f:z:i:t:o:e:l:L:m:a:n:r:R:s:g:d:xy \
    --long help,MIR,DEEP,genome,blast,nbreads,mirnas,predict,valid,file:,text:,output:,errors:,Length:,length:,mfei:,amfe:,non-app:,ratiomin:,ratiomax:,stem:,left:,right:\
    -- "$@"`

if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi

# Note the quotes around `$TEMP': they are essential!
eval set -- "$TEMP"

src=${0%%/MIReNA.sh}

valoption=""; file=""; text="text.txt"; genome="genome.fa";mirnas="miRNAs.fa";
errors="15"; length="60"; stem="0"; nonapp="0.26";
ratiomin="0.83"; ratiomax="1.17"; output="out.txt"; left="200"; right="200";
Length="25"; nbreads="3";

while true
do
  case "$1" in
    -M|--MIR)     option="MIR"; file="miRNAs.fa"; shift;;
    -D|--DEEP)    option="DEEP"; file="deep.fa"; shift;;
    -p|--predict) option="predict"; file="input.txt"; shift;;
    -v|--valid)   option="valid"; file="preMi.fa"; shift;;
    -i|--input)   input=$2; shift 2;;
    -f|--file)    file=$2; shift 2;;
    -t|--text)    text=$2; shift 2;;
    -b|--blast)   blast=$2; shift 2;;
    -z|--nbreads) nbreads=$2; shift 2;;
    -j|--genome)  genome=$2; shift 2;;
    -k|--mirnas)  mirnas=$2; shift 2;;
    -o|--output)  output=$2; shift 2;;
    -e|--errors)  errors=$2; shift 2;;
    -l|--length)  length=$2; shift 2;;
    -L|--Length)  Length=$2; shift 2;;
    -g|--left)    left=$2; shift 2;;
    -d|--right)   right=$2; shift 2;;
    -m|--mfei)    mfei=$2; shift 2;;
    -a|--amfe)    amfe=$2; shift 2;;
    -n|--non-app) nonapp=$2; shift 2;;
    -r|--ratiomin)ratiomin=$2; shift 2;;
    -R|--ratiomax)ratiomax=$2; shift 2;;
    -s|--stem)    stem=$2; shift 2;;
    -y       )    valoption="-p"; shift;;
    -x       )    valoption="-x"; shift;;
    -h|--help)    echo -e "$usage"; exit 1;;
    --)           shift; break;;
    *)            echo "Internal error!"; exit 1;;
  esac
done

checkRNAfold;

case "$option" in
  "MIR")
            echo "Launching MIReNA from miRNA sequences"
            if [ ! $mfei ]
            then
              mfei="-0.85";
            fi
            if [ ! $amfe ]
            then
              amfe="-32"
            fi
            checkMIR $file $text $output $errors $Length $mfei $amfe\
              $nonapp $ratiomin $ratiomax $stem $left $right $length;
            TMPFILE1=`mktemp`;
            echo -n "Searching for similar sequences ... "
            bash $src/fromMirnas/executeAsml.sh -f $file -t $text -e $errors -l\
              $Length -o $TMPFILE1;
            echo "Done"
            TMPFILE2=`mktemp`;
            echo -n "Extracting the extended sequences ... "
            perl $src/fromMirnas/getPutative.pl $left $right $TMPFILE1 $TMPFILE2;
            echo "Done"
            rm $TMPFILE1;
            TMPFILE1=`mktemp`;
            echo -n "Looking for pre-miRNAs validating the criteria ... "
            perl $src/fromMirnas/makeFold.pl -1 $TMPFILE2 $amfe $mfei $nonapp $ratiomin\
              $ratiomax $stem $length $TMPFILE1;
            echo "Done"
            echo -n "Selection between overlapping miRNAs ... "
            python $src/fromMirnas/elim_doubl.py -f $TMPFILE1 -o $output;
            echo "Done"
            rm $TMPFILE1 $TMPFILE2;
            echo "MIReNA has terminated successfully";;
  "DEEP")
            echo "Launching MIReNA from deep sequencing data"
            if [ ! $mfei ]
            then
              mfei="-0.62"
            fi
            if [ ! $amfe ]
            then
              amfe="-22"
            fi
            checkDEEP $file $genome $output $mirnas $mfei $amfe $nonapp $ratiomin\
              $ratiomax $stem $nbreads;
            if [ ! $blast ]
            then
              blastout=`mktemp`
              echo -n "Making the megablast search ... "
              formatdb -i $genome -p F -o T;
              megablast -d $genome -i $file -W 12 -D 2 -p 100 -o $blastout;
              echo "Done"
              echo -n "Parsing the megablast output ... "
              perl $src/miRDeep/blastoutparse.pl $blastout > blastoutparsed;
              echo "Done"
            else
              cp $blast blastoutparsed;
            fi
            echo -n "Excising the putative pre-miRNAs ... "
            perl $src/miRDeep/filter_alignments.pl blastoutparsed -c 5 >\
              blastoutparsed_excision;
            perl $src/miRDeep/excise_candidate.pl $genome blastoutparsed_excision >\
              precursors.fa;
            echo "Done"
            echo -n "Mapping the reads on the pre-miRNAs ... "
            perl $src/miRDeep/filter_alignments.pl blastoutparsed_excision -b $file\
              > $file"_aligned";
            perl $src/miRDeep/auto_blast.pl $file"_aligned" precursors.fa -b >\
              signatures;
            echo "Done"
            echo -n "Computing secondary structures ... "
            RNAfold --noPS < precursors.fa > precursors.fold;
            echo "Done"
            echo -n "Validation of pre-miRNAs through the criteria ... "
            perl $src/fromDeepSeq/make_precursors_with_reads.pl -p precursors.fold \
              -s signatures > precursors_reads.fold;
            perl $src/testValid.pl -p -i precursors_reads.fold -s $nonapp \
              -r $ratiomin -R $ratiomax -a $amfe -m $mfei -n $stem >\
              precursors_reads.valid;
            perl $src/fromDeepSeq/make_potential_precursors_signature.pl -s signatures \
              -p precursors_reads.valid -o sig -a prec.valid;
            echo "Done"
            echo -n "Looking for Dicer processing conditions ... "
            perl $src/fromDeepSeq/dicer_signature.pl sig prec.valid -s $mirnas \
              -z $nbreads > pred;
            echo "Done"
            echo -n "Extracting the predictions ... "
            perl $src/fromDeepSeq/obtain_predicted.pl pred > $output;
            echo "Done"
            rm -f blastoutparsed* $file"_aligned" sig signatures pred \
              precursors_reads* prec.valid;
            if [ ! $blast ]
            then
              rm $blastout $genome.n*;
            fi;
            echo "MIReNA has terminated successfully";;
  "predict") 
            echo "Launching MIReNA for pre-miRNA predictions"
            if [ ! $mfei ]
            then
              mfei="-0.69"
            fi
            if [ ! $amfe ]
            then
              amfe="-32"
            fi
            checkPredict $file $output $mfei $amfe $nonapp $ratiomin\
              $ratiomax $stem $length;
            echo -n "Looking for pre-miRNAs validating the criteria ... "
            perl $src/fromMirnas/makeFold.pl -2 $file $amfe $mfei $nonapp $ratiomin $ratiomax\
              $stem $length $output;
            echo "Done";
            echo "MIReNA has terminated successfully";;
  "valid") 
            echo "Launching MIReNA for pre-miRNAs validation"
            if [ ! $mfei ]
            then
              mfei="-0.69"
            fi
            if [ ! $amfe ]
            then
              amfe="-32"
            fi
            checkValid $file $output $mfei $amfe $nonapp $ratiomin\
            $ratiomax $stem $length $valoption;
            #TMPFILE=`mktemp`;
            echo -n "Computing the secondary structures ... "
            RNAfold --noPS < $file > RNAfold.txt;
            echo "Done"
            #rm *.ps;
            echo -n "Testing the validation of the criteria ... "
            perl $src/testValid.pl $valoption -i RNAfold.txt -s $nonapp -r $ratiomin -R $ratiomax -a $amfe -m $mfei \
              -n $stem -l $length > $output;
            echo "Done";
            echo "MIReNA has terminated successfully";;
  "")   echo -e "$usage"; exit 1;;
esac
