#!/usr/bin/perl

use warnings;
use strict;

use Getopt::Std;

my $usage ="
$0 blastoutputfile

This script parses the output of any blastn, blastx or tblastn file.

The resulting blastparsed format contains most of the information from ncbi blast output in tabulated fields.
One line corresponds to one local alignment:

The first field is the query id
The second field is the query length
The third field is the range of nucleotides in the query covered by the alignment
The fourth field is the subject id
The fifth field is the subject length
The sixth field is the range of nucleotides in the subject covered by the alignment
The seventh field is the blast e-value
The eigth field is the percentage of identical nucleotides in the alignment
The ninth field is the bit score
The tenth field (can contain whitespaces) are the strands of the query/subject in the alignment 

Example:

fire_9_x1       22      1..22   chrII   15279308        11537795..11537816      4e-05   1.00    44.1    Plus / Plus
";


my $query;
my $query_lng;
my $subject;
my $subject_lng;
my $bits;
my $expect;
my $identities;
my $gaps;
my $frame;
my $q_st;
my $q_ed;
my $s_st;
my $s_ed;




while (<>) {
    
  X:if (/^Query\s*=\s*(\S+)/) {
      #resetting variables to catch parsing errors
     $query_lng="";
     
     $query=$1;
 }
    
    if(/\s+\(([\d\,]+)\s+letters\)/) {
	$query_lng=$1;
    }
    

    #Requiring that query_lng is defined is to avoid parsing parts of the query descriptor
    #containing '>' as being the subject and going into the do while loop. Is currently a test.
    if (/^>(\S+)/ && $query_lng) {
	#resetting variables to catch parsing errors
	$subject_lng="";
	
	$subject=$1;
	
	
     	do {
     	    $_ = <>;
     	    tr/,//d;
	    
	    if (/\s+Length\s*=\s*(\d+)/) {$subject_lng=$1;}
	    
  	} while (!/\s+Length\s*=\s*(\d+)/);
    }
    
    
    if(/^\s{0,2}Score/){
	if (/^\s+Score\s*=\s*(\S+)\s*.*Expect\S*\s*=\s*(\S+)$/) {
	    $bits=$1;
	    $expect=$2;
	}else{
	    print STDERR "problem parsing score and expect line\n";
	    exit;
	}

	$_ = <>;
	if (!/Identities\s*=\s*(\S+)/) {
	    print STDERR "problem parsing identities line\n";
	    exit;
	}
	else { $identities=$1; }
	
	if (/Gaps\s*=\s*(\d+)\//) { $gaps = $1; }
	else                      { $gaps = 0;  }
	
	
	$_ = <>;
	if (/\s+Frame\s*=\s*(\S+)\s*\n/){
	    $frame=$1;
	}elsif(/\s*Strand\s*=\s*(\w+\s*\/\s*\w+)/){
	    $frame=$1;
	}else{
	  print STDERR "error reading the strand or frame\n";
      } 
	
	# read positions
	$q_st = "";
	$q_ed = "";
	$s_st = "";
	$s_ed = "";
	do {
	    $_ = <>;
	    if (/^Query:\s+(\d+).*\s(\d+)\s*$/)  {
		if ($q_st eq "")                 { $q_st = $1; }
		if ($1 < $q_st)                  { $q_st = $1; }
		if ($2 < $q_st)                  { $q_st = $2; }
		if ($q_ed eq "" or $1 > $q_ed)   { $q_ed = $1; }
		if ($2 > $q_ed)                  { $q_ed = $2; }
	    }
	    if (/^Sbjct:\s+(\d+).*\s(\d+)\s*$/)  {
		if ($s_st eq "")	         { $s_st = $1; }
		if ($1 < $s_st)		         { $s_st = $1; }
		if ($2 < $s_st)		         { $s_st = $2; }
		if ($s_ed eq "" or $1 > $s_ed)   { $s_ed = $1; }
		if ($2 > $s_ed)		         { $s_ed = $2; }
	    }
	} while ((/^\s*$/ or /^\s{8,}(\w|\s|\+|\||\*)+$/ or /^Query:/ or /^Sbjct:/) and !/^>/);
	
	my $strand=findstrand($frame);
	if($strand eq "-"){($s_st,$s_ed)=reversepositions($subject_lng,$s_st,$s_ed);}

	
	my @v = split(/\//, $identities);
	my $pid = $v[0]/$v[1];
	
	print "$query\t$query_lng\t$q_st..$q_ed\t$subject\t$subject_lng\t$s_st..$s_ed\t$expect\t";
	printf "%1.2f", $pid;
	print "\t$bits\t$frame";
	print "\n";
	
	unless($query and $query_lng and $q_st and $q_ed and $subject and $subject_lng and $s_st and $s_ed and $expect and $pid and $bits){
	    print STDERR "problem parsing entry\n";
	    exit;
	}
	
	#resetting variables to catch parsing errors
	$bits="";
	$expect="";
	$identities="";
	$gaps="";
	$frame="";
	$q_st="";
	$q_ed="";
	$s_st="";
	$s_ed="";

	goto X;
	
    }
}




sub findstrand{

    #A subroutine to find the strand, parsing different blast formats
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








sub reversepositions{

    #A subroutine to find positions relative to the minus strand
    my($length,$begin,$end)=@_;

    my $new_end=$length-$begin+1;
    my $new_beg=$length-$end+1;

    return($new_beg,$new_end);
}
