#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Std;



################################# Script obtained by modifying MIRDEEP ####################

################################## USAGE ##################################################


my $usage=
"$0 <signatures> <structures> -s <known miRNAs> -z <min reads for cond. 3.>";

############################################################################################

################################### INPUT ##################################################


#signature file in blast_parsed format
my $file_blast_parsed=shift or die $usage;

#structure file outputted from RNAfold
my $file_struct=shift or die $usage;

#options
my %options=();
getopts("hs:z:",\%options);






#############################################################################################

############################# GLOBAL VARIABLES ##############################################


#parameters
my $nucleus_lng=7;

#hashes
my %hash_desc;
my %hash_seq;
my %hash_struct;
my %hash_mfe;
my %hash_nuclei;
my %hash_mirs;
my %hash_query;
my %hash_comp;
my %hash_bp;

#other variables
my $subject_old;
my $message_filter;
my $message_score;
my $lines;
my $out_of_bound;
### Modif made by A. Mathelier ###
my $mature;
my $diff_reads;
### End modif by A. Mathelier ###




##############################################################################################

################################  MAIN  ###################################################### 


#print help if that option is used
if($options{h}){die $usage;}

#parse structure file outputted from RNAfold
parse_file_struct($file_struct);

#if conservation is scored, the fasta file of known miRNA sequences is parsed
if($options{s}){create_hash_nuclei($options{s})};

### Modif by A. Mathelier ###
if($options{z}){
  $diff_reads = $options{z};
}else{
  $diff_reads = 4;
}
### End modif by Mathelier ###

#parse signature file in blast_parsed format and resolve each potential precursor
parse_file_blast_parsed($file_blast_parsed);

exit;




##############################################################################################

############################## SUBROUTINES ###################################################



sub parse_file_blast_parsed{

#    read through the signature blastparsed file, fills up a hash with information on queries
#    (deep sequences) mapping to the current subject (potential precursor) and resolve each
#    potential precursor in turn
 
    my $file_blast_parsed=shift;
    
    open (FILE_BLAST_PARSED, "<$file_blast_parsed") or die "can not open $file_blast_parsed\n";
    while (my $line=<FILE_BLAST_PARSED>){

	if($line=~/^(\S+)\s+(\S+)\s+(\d+)\.+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\.+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(.+)$/){
          
       	    my $query=$1;
	    my $query_lng=$2;
	    my $query_beg=$3;
	    my $query_end=$4;
	    my $subject=$5;
	    my $subject_lng=$6;
	    my $subject_beg=$7;
	    my $subject_end=$8;
	    my $e_value=$9;
	    my $pid=$10;
	    my $bitscore=$11;
	    my $other=$12;
          ### Modif made by A. Mathelier ###
            my $is_mature = 0;
            if($line =~ /\s+mature/){
              $is_mature = 1;
              $line =~ s/\s+mature//;
            }
          ### End modif by A. Mathelier ###
	    
	    #if the new line concerns a new subject (potential precursor) then the old subject must be resolved
	    if($subject_old and $subject_old ne $subject){
		resolve_potential_precursor();
	    }

	    #resolve the strand
	    my $strand=find_strand($other);

	    #resolve the number of reads that the deep sequence represents
	    my $freq=find_freq($query);

	    #read information of the query (deep sequence) into hash
	    $hash_query{$query}{"subject_beg"}=$subject_beg;
	    $hash_query{$query}{"subject_end"}=$subject_end;
	    $hash_query{$query}{"strand"}=$strand;
	    $hash_query{$query}{"freq"}=$freq;

            ### Modif made by A. Mathelier ###
            if($is_mature){
              $mature = $query;
            }
            ### End modif by A. Mathelier ###

	    #save the signature information
	    $lines.=$line;

	    $subject_old=$subject;
	}
    }
    resolve_potential_precursor();
}

sub resolve_potential_precursor{
    
#    dissects the potential precursor in parts by filling hashes, and tests if it passes the
#    initial filter and the scoring filter

#    binary variable whether the potential precursor is still viable
    my $ret=1;

    fill_structure();
    
    fill_pri();

    fill_mature();
   
    fill_star();

    fill_loop();

    fill_lower_flanks();

### Modif made by A. Mathelier ###
    unless(pass_filtering_initial()){$ret=0;}
    if($ret){
      unless($hash_comp{"star_perfect"} || 
        ((not test_nucleus_conservation()) && 
          $hash_comp{"diff_consistent"} >= $diff_reads) || 
        (test_nucleus_conservation() && $hash_comp{"mature_fuzzy"} >= 2)){$ret=0;}

    }
### End modif by A. Mathelier ###

    print_results($ret);
    
    reset_variables();
    
    return;
    
}



sub print_results{

    my $ret=shift;

#    print out if the precursor is accepted and accepted precursors should be printed out
#    or if the potential precursor is discarded and discarded potential precursors should
#    be printed out
	
    if((!$options{t} and $ret) or ($options{t} and !$ret)){
	#full output
	unless($options{u}){
	    if($message_filter){print $message_filter;}
	    if($message_score){print $message_score;}
	    print_hash_comp();
	    print $lines."\n\n";
	    return;
	}
	#limited output (only ids)
	my $id=$hash_comp{"pri_id"};
	print "$id\n";
    }    
}



sub test_nucleus_conservation{

    #test if nucleus is identical to nucleus from known miRNA, return 1 or 0

    my $nucleus=substr($hash_comp{"mature_seq"},1,$nucleus_lng);
    $nucleus=~tr/[T]/[U]/;
    if($hash_nuclei{$nucleus}){return 1;}
    
    return 0;
}



sub pass_filtering_initial{

    #test if the structure forms a plausible hairpin
    unless(pass_filtering_structure()){filter_p("structure problem"); return 0;}

    #test if >90% of reads map to the hairpin in consistence with Dicer processing
    unless(pass_filtering_signature()){filter_p("signature problem");return 0;}

    return 1;

}


sub pass_filtering_signature{

    #number of reads that map in consistence with Dicer processing
    my $consistent=0;

    #number of reads that map inconsistent with Dicer processing
    my $inconsistent=0;
   
#   number of potential star reads map in good consistence with Drosha/Dicer processing
#   (3' overhangs relative to mature product)
    my $star_perfect=0;

#   number of potential star reads that do not map in good consistence with 3' overhang
    my $star_fuzzy=0;
 

    #sort queries (deep sequences) by their position on the hairpin
    my @queries=sort {$hash_query{$a}{"subject_beg"} <=> $hash_query{$b}{"subject_beg"}} keys %hash_query;

    ### Modif by A. Mathelier ###
    my %hash_diff;
    my %hash_diff_consistent;
    $hash_comp{"mature_fuzzy"} = 0;
    $hash_comp{"star_fuzzy"} = 0;
    ### End modif by A. Mathelier ###
    foreach my $query(@queries){

	#number of reads that the deep sequence represents
	unless(defined($hash_query{$query}{"freq"})){next;}
	my $query_freq=$hash_query{$query}{"freq"};

        ### Modif made by A. Mathelier ###
        my $query_beg = $hash_query{$query}{"subject_beg"};
        my $query_end = $hash_query{$query}{"subject_end"};
        unless(defined($hash_diff{$query_beg.$query_end})){
          $hash_diff{$query_beg.$query_end} = 1;
        }
        unless(defined($hash_diff_consistent{$query_beg.$query_end})){
          if(test_query($query)){
            $hash_diff_consistent{$query_beg.$query_end} = 1;
          }
        }
        ### End modif by A. Mathelier ###

	#test which Dicer product (if any) the deep sequence corresponds to
	my $product=test_query($query);

	#if the deep sequence corresponds to a Dicer product, add to the 'consistent' variable
	if($product){$consistent+=$query_freq;}

	#if the deep sequence do not correspond to a Dicer product, add to the 'inconsistent' variable
	else{$inconsistent+=$query_freq;}

	#test a potential star sequence has good 3' overhang
	if($product eq "star"){
        ### Modif made by A. Mathelier ###
          $hash_comp{"star_fuzzy"} += 1;
        ### End modif by A. Mathelier ###
	    if(test_star($query)){$star_perfect+=$query_freq;}
	    else{$star_fuzzy+=$query_freq;}
	}

        ### Modif made by A. Mathelier ###
        if($product eq "mature"){
          $hash_comp{"mature_fuzzy"} += 1;
        }
        ### End modif by A. Mathelier ###
    }
    ### Modif by A. Mathelier ###
    $hash_comp{"diff"} = scalar(keys(%hash_diff));
    $hash_comp{"diff_consistent"} = scalar(keys(%hash_diff_consistent));

    if($star_perfect){
      $hash_comp{"star_perfect"} = 1;
    }
    ### End modif by A. Mathelier ###

#   if the majority of potential star sequences map in good accordance with 3' overhang
#    score for the presence of star evidence
#    print "star_perfect: $star_perfect\n";
#    print "star_fuzzy: $star_fuzzy\n";
    if($star_perfect>$star_fuzzy){$hash_comp{"star_read"}=1;}

    #total number of reads mapping to the hairpin
    my $freq=$consistent+$inconsistent;
    $hash_comp{"freq"}=$freq;
    unless($freq>0){filter_s("read frequency too low"); return 0;}

    #unless >90% of the reads map in consistence with Dicer processing, the hairpin is discarded
    my $inconsistent_fraction=$inconsistent/($inconsistent+$consistent);
    unless($inconsistent_fraction<=0.1){filter_p("inconsistent\t$inconsistent\nconsistent\t$consistent"); return 0;}

    #the hairpin is retained
    return 1;
}

sub test_star{

    #test if a deep sequence maps in good consistence with 3' overhang

    my $query=shift;

    #5' begin and 3' end positions
    my $beg=$hash_query{$query}{"subject_beg"};
    my $end=$hash_query{$query}{"subject_end"};

    #the difference between observed and expected begin positions must be 0 or 1
    my $offset=$beg-$hash_comp{"star_beg"};
    if($offset==0 or $offset==1 or $offset==-1){return 1;}

    return 0;
}



sub test_query{

    #test if deep sequence maps in consistence with Dicer processing
    
    my $query=shift;

    #begin, end, strand and read count
    my $beg=$hash_query{$query}{"subject_beg"};
    my $end=$hash_query{$query}{"subject_end"};
    my $strand=$hash_query{$query}{"strand"};
    my $freq=$hash_query{$query}{"freq"};

    #should not be on the minus strand (although this has in fact anecdotally been observed for known miRNAs)
    if($strand eq '-'){return 0;}

    #the deep sequence is allowed to stretch 2 nt beyond the expected 5' end
    my $fuzz_beg=2;
    #the deep sequence is allowed to stretch 5 nt beyond the expected 3' end
    my $fuzz_end=5;

    #if in accordance with Dicer processing, return the type of Dicer product 
    if(contained($beg,$end,$hash_comp{"mature_beg"}-$fuzz_beg,$hash_comp{"mature_end"}+$fuzz_end)){return "mature";}
    if(contained($beg,$end,$hash_comp{"star_beg"}-$fuzz_beg,$hash_comp{"star_end"}+$fuzz_end)){return "star";}
    if(contained($beg,$end,$hash_comp{"loop_beg"}-$fuzz_beg,$hash_comp{"loop_end"}+$fuzz_end)){return "loop";}
  
    #if not in accordance, return 0
    return 0;
}


sub pass_filtering_structure{

    #The potential precursor must form a hairpin with miRNA precursor-like characteristics

    #return value
    my $ret=1;

    #potential mature, star, loop and lower flank parts must be identifiable
    unless(test_components()){return 0;}

    return $ret;
}






sub test_components{

    #tests whether potential mature, star, loop and lower flank parts are identifiable

    unless($hash_comp{"mature_struct"}){
	filter_s("no mature");
	return 0;
    }

    unless($hash_comp{"star_struct"}){
	filter_s("no star");
	return 0;
    }

    unless($hash_comp{"loop_struct"}){
	filter_s("no loop");
   	return 0;
    }

    unless($hash_comp{"flank_first_struct"}){
	filter_s("no flanks");
   	return 0;
    }
     
    unless($hash_comp{"flank_second_struct"}){
	filter_s("no flanks");
    	return 0;
    }
    return 1;
}



sub fill_structure{

    #reads the dot bracket structure into the 'bp' hash where each key and value are basepaired

    my $struct=$hash_struct{$subject_old};
    my $lng=length $struct;

    #local stack for keeping track of basepairings
    my @bps;

    for(my $pos=1;$pos<=$lng;$pos++){
	my $struct_pos=excise_struct($struct,$pos,$pos,"+");

	if($struct_pos eq "("){
	    push(@bps,$pos);
	}

	if($struct_pos eq ")"){
	    my $pos_prev=pop(@bps);
	    $hash_bp{$pos_prev}=$pos;
	    $hash_bp{$pos}=$pos_prev;
	}
    }
    return;
}



sub fill_star{

    #fills specifics on the expected star strand into 'comp' hash ('component' hash)
    
    #if the mature sequence is not plausible, don't look for the star arm
    my $mature_arm=$hash_comp{"mature_arm"};
    unless($mature_arm){$hash_comp{"star_arm"}=0; return;}
 
    #if the star sequence is not plausible, don't fill into the hash
    my($star_beg,$star_end)=find_star();
    my $star_arm=arm_star($star_beg,$star_end);
    unless($star_arm){return;}

    #excise expected star sequence and structure
    my $star_seq=excise_seq($hash_comp{"pri_seq"},$star_beg,$star_end,"+");
    my $star_struct=excise_seq($hash_comp{"pri_struct"},$star_beg,$star_end,"+");

    #fill into hash
    $hash_comp{"star_beg"}=$star_beg;
    $hash_comp{"star_end"}=$star_end;
    $hash_comp{"star_seq"}=$star_seq;
    $hash_comp{"star_struct"}=$star_struct;
    $hash_comp{"star_arm"}=$star_arm;

    return;
}


sub find_star{

    #uses the 'bp' hash to find the expected star begin and end positions from the mature positions

    #the -2 is for the overhang
    my $mature_beg=$hash_comp{"mature_beg"};
    my $mature_end=$hash_comp{"mature_end"}-2;
    my $mature_lng=$mature_end-$mature_beg+1;

    #in some cases, the last nucleotide of the mature sequence does not form a base pair,
    #and therefore does not basepair with the first nucleotide of the star sequence.
    #In this case, the algorithm searches for the last nucleotide of the mature sequence
    #to form a base pair. The offset is the number of nucleotides searched through.
    my $offset_star_beg=0;
    my $offset_beg=0;

    #the offset should not be longer than the length of the mature sequence, then it
    #means that the mature sequence does not form any base pairs
    while(!$offset_star_beg and $offset_beg<$mature_lng){
	if($hash_bp{$mature_end-$offset_beg}){
	    $offset_star_beg=$hash_bp{$mature_end-$offset_beg};
	}else{
	    $offset_beg++;
	}
    }
    #when defining the beginning of the star sequence, compensate for the offset
    my $star_beg=$offset_star_beg-$offset_beg;

    #same as above
    my $offset_star_end=0;
    my $offset_end=0;
    while(!$offset_star_end and $offset_end<$mature_lng){
	if($hash_bp{$mature_beg+$offset_end}){
	    $offset_star_end=$hash_bp{$mature_beg+$offset_end};
	}else{
	    $offset_end++;
	}
    }
    #the +2 is for the overhang
    my $star_end=$offset_star_end+$offset_end+2;

    return($star_beg,$star_end);
}


sub fill_pri{

    #fills basic specifics on the precursor into the 'comp' hash
    
    my $seq=$hash_seq{$subject_old};
    my $struct=$hash_struct{$subject_old};
    my $mfe=$hash_mfe{$subject_old};
    my $length=length $seq;
    
    $hash_comp{"pri_id"}=$subject_old;
    $hash_comp{"pri_seq"}=$seq;
    $hash_comp{"pri_struct"}=$struct;
    $hash_comp{"pri_mfe"}=$mfe;
    $hash_comp{"pri_beg"}=1;
    $hash_comp{"pri_end"}=$length;
    
    return;
}


sub fill_mature{

    #fills specifics on the mature sequence into the 'comp' hash

    ### Modif made by A. Mathelier ###
    #my $mature_query=find_mature_query();
    my $mature_query = $mature;
    ### End modif by A. Mathelier ###
    my($mature_beg,$mature_end)=find_positions_query($mature_query);
    my $mature_strand=find_strand_query($mature_query);
    my $mature_seq=excise_seq($hash_comp{"pri_seq"},$mature_beg,$mature_end,$mature_strand);
    my $mature_struct=excise_struct($hash_comp{"pri_struct"},$mature_beg,$mature_end,$mature_strand);
    my $mature_arm=arm_mature($mature_beg,$mature_end,$mature_strand);

    $hash_comp{"mature_query"}=$mature_query;
    $hash_comp{"mature_beg"}=$mature_beg;
    $hash_comp{"mature_end"}=$mature_end;
    $hash_comp{"mature_strand"}=$mature_strand;
    $hash_comp{"mature_struct"}=$mature_struct;
    $hash_comp{"mature_seq"}=$mature_seq;
    $hash_comp{"mature_arm"}=$mature_arm;

    return;
}



sub fill_loop{

    #fills specifics on the loop sequence into the 'comp' hash

    #unless both mature and star sequences are plausible, do not look for the loop
    unless($hash_comp{"mature_arm"} and $hash_comp{"star_arm"}){return;}

    my $loop_beg;
    my $loop_end;

    #defining the begin and end positions of the loop from the mature and star positions
    #excision depends on whether the mature or star sequence is 5' of the loop ('first')
    if($hash_comp{"mature_arm"} eq "first"){
	$loop_beg=$hash_comp{"mature_end"}+1;
    }else{
	$loop_end=$hash_comp{"mature_beg"}-1;
    }
    
    if($hash_comp{"star_arm"} eq "first"){
	$loop_beg=$hash_comp{"star_end"}+1;
    }else{
	$loop_end=$hash_comp{"star_beg"}-1;
    }

    #unless the positions are plausible, do not fill into hash
    unless(test_loop($loop_beg,$loop_end)){return;}

    my $loop_seq=excise_seq($hash_comp{"pri_seq"},$loop_beg,$loop_end,"+");
    my $loop_struct=excise_struct($hash_comp{"pri_struct"},$loop_beg,$loop_end,"+");

    $hash_comp{"loop_beg"}=$loop_beg;
    $hash_comp{"loop_end"}=$loop_end;
    $hash_comp{"loop_seq"}=$loop_seq;
    $hash_comp{"loop_struct"}=$loop_struct;

    return;
}


sub fill_lower_flanks{

    #fills specifics on the lower flanks and unpaired strands into the 'comp' hash

    #unless both mature and star sequences are plausible, do not look for the flanks
    unless($hash_comp{"mature_arm"} and $hash_comp{"star_arm"}){return;}

    my $flank_first_end;
    my $flank_second_beg;

    #defining the begin and end positions of the flanks from the mature and star positions
    #excision depends on whether the mature or star sequence is 5' in the potenitial precursor ('first')
    if($hash_comp{"mature_arm"} eq "first"){
	$flank_first_end=$hash_comp{"mature_beg"}-1;
    }else{
	$flank_second_beg=$hash_comp{"mature_end"}+1;
    }
    
    if($hash_comp{"star_arm"} eq "first"){
	$flank_first_end=$hash_comp{"star_beg"}-1;
    }else{
	$flank_second_beg=$hash_comp{"star_end"}+1;
    }   

    #unless the positions are plausible, do not fill into hash
    unless(test_flanks($flank_first_end,$flank_second_beg)){return;}

    $hash_comp{"flank_first_end"}=$flank_first_end;
    $hash_comp{"flank_second_beg"}=$flank_second_beg;
    $hash_comp{"flank_first_seq"}=excise_seq($hash_comp{"pri_seq"},$hash_comp{"pri_beg"},$hash_comp{"flank_first_end"},"+");
    $hash_comp{"flank_second_seq"}=excise_seq($hash_comp{"pri_seq"},$hash_comp{"flank_second_beg"},$hash_comp{"pri_end"},"+");
    $hash_comp{"flank_first_struct"}=excise_struct($hash_comp{"pri_struct"},$hash_comp{"pri_beg"},$hash_comp{"flank_first_end"},"+");
    $hash_comp{"flank_second_struct"}=excise_struct($hash_comp{"pri_struct"},$hash_comp{"flank_second_beg"},$hash_comp{"pri_end"},"+");

    return;
}


sub arm_mature{
 
    #tests whether the mature sequence is in the 5' ('first') or 3' ('second') arm of the potential precursor

    my ($beg,$end,$strand)=@_;
 
    #mature and star sequences should alway be on plus strand
    if($strand eq "-"){return 0;}

    #there should be no bifurcations and minimum one base pairing
    my $struct=excise_seq($hash_comp{"pri_struct"},$beg,$end,$strand);
    if(defined($struct) and $struct=~/^(\(|\.)+$/ and $struct=~/\(/){
	return "first";
    }elsif(defined($struct) and $struct=~/^(\)|\.)+$/ and $struct=~/\)/){
	return "second";
    }
    return 0;
}


sub arm_star{

    #tests whether the star sequence is in the 5' ('first') or 3' ('second') arm of the potential precursor

    my ($beg,$end)=@_;

    #unless the begin and end positions are plausible, test negative
    unless($beg>0 and $beg<=$hash_comp{"pri_end"} and $end>0 and $end<=$hash_comp{"pri_end"} and $beg<=$end){return 0;}

    #no overlap between the mature and the star sequence
    if($hash_comp{"mature_arm"} eq "first"){
	($hash_comp{"mature_end"}<$beg) or return 0;
    }elsif($hash_comp{"mature_arm"} eq "second"){
	($end<$hash_comp{"mature_beg"}) or return 0;
    }

    #there should be no bifurcations and minimum one base pairing
    my $struct=excise_seq($hash_comp{"pri_struct"},$beg,$end,"+");
    if($struct=~/^(\(|\.)+$/ and $struct=~/\(/){
	return "first";
    }elsif($struct=~/^(\)|\.)+$/ and $struct=~/\)/){
	return "second";
    }
    return 0;
}


sub test_loop{

    #tests the loop positions

    my ($beg,$end)=@_;

    #unless the begin and end positions are plausible, test negative
    unless($beg>0 and $beg<=$hash_comp{"pri_end"} and $end>0 and $end<=$hash_comp{"pri_end"} and $beg<=$end){return 0;}

    return 1;
}


sub test_flanks{

    #tests the positions of the lower flanks

    my ($beg,$end)=@_;

    #unless the begin and end positions are plausible, test negative
    unless($beg>0 and $beg<=$hash_comp{"pri_end"} and $end>0 and $end<=$hash_comp{"pri_end"} and $beg<=$end){return 0;}

    return 1;
}


sub comp{

    #subroutine to retrive from the 'comp' hash

    my $type=shift;
    my $component=$hash_comp{$type};
    return $component;
}


sub find_strand_query{

    #subroutine to find the strand for a given query

    my $query=shift;
    my $strand=$hash_query{$query}{"strand"};
    return $strand;
}


sub find_positions_query{

    #subroutine to find the begin and end positions for a given query

    my $query=shift;
    my $beg=$hash_query{$query}{"subject_beg"};
    my $end=$hash_query{$query}{"subject_end"};
    return ($beg,$end);
}



sub find_mature_query{

    #finds the query with the highest frequency of reads and returns it
    #is used to determine the positions of the potential mature sequence

    my @queries=sort {$hash_query{$b}{"freq"} <=> $hash_query{$a}{"freq"}} keys %hash_query;
    my $mature_query=$queries[0];
    return $mature_query;
}




sub reset_variables{

    #resets the hashes for the next potential precursor

    %hash_query=();
    %hash_comp=();
    %hash_bp=();

    $message_filter=();
    $message_score=();
    $lines=();

    return;
}



sub excise_seq{

    #excise sub sequence from the potential precursor

    my($seq,$beg,$end,$strand)=@_;

    #begin can be equal to end if only one nucleotide is excised
    unless($beg<=$end){print STDERR "begin can not be smaller than end for $subject_old\n";exit;}

    #rarely, permuted combinations of signature and structure cause out of bound excision errors.
    #this happens once appr. every two thousand combinations
    unless($beg<=length($seq)){$out_of_bound++;return 0;}

    #if on the minus strand, the reverse complement should be excised
    if($strand eq "-"){$seq=revcom($seq);}

    #the blast parsed format is 1-indexed, substr is 0-indexed
    my $sub_seq=substr($seq,$beg-1,$end-$beg+1);

    return $sub_seq;

}

sub excise_struct{

    #excise sub structure

    my($struct,$beg,$end,$strand)=@_;
    my $lng=length $struct;

    #begin can be equal to end if only one nucleotide is excised
    unless($beg<=$end){print STDERR "begin can not be smaller than end for $subject_old\n";exit;}

    #rarely, permuted combinations of signature and structure cause out of bound excision errors.
    #this happens once appr. every two thousand combinations
    unless($beg<=length($struct)){return 0;}

    #if excising relative to minus strand, positions are reversed
    if($strand eq "-"){($beg,$end)=rev_pos($beg,$end,$lng);}

    #the blast parsed format is 1-indexed, substr is 0-indexed
    my $sub_struct=substr($struct,$beg-1,$end-$beg+1);
 
    return $sub_struct;
}


sub create_hash_nuclei{
 
    #parses a fasta file with sequences of known miRNAs considered for conservation purposes
    #reads the nuclei into a hash

    my ($file) = @_;
    my ($id, $desc, $sequence, $nucleus) = ();

    open (FASTA, "<$file") or die "can not open $file\n";
    while (<FASTA>)
    {
        chomp;
        if (/^>(\S+)(.*)/)
	{
	    $id       = $1;
	    $desc     = $2;
	    $sequence = "";
	    $nucleus  = "";
	    while (<FASTA>){
                chomp;
                if (/^>(\S+)(.*)/){
		    $nucleus                = substr($sequence,1,$nucleus_lng);
		    $nucleus                =~ tr/[T]/[U]/;
		    $hash_mirs{$nucleus}   .="$id\t";
		    $hash_nuclei{$nucleus} += 1;

		    $id               = $1;
		    $desc             = $2;
		    $sequence         = "";
		    $nucleus          = "";
		    next;
                }
		$sequence .= $_;
            }
        }
    }
    $nucleus                = substr($sequence,1,$nucleus_lng);
    $nucleus                =~ tr/[T]/[U]/;
    $hash_mirs{$nucleus}   .="$id\t";
    $hash_nuclei{$nucleus} += 1;
    close FASTA;
}
    

sub parse_file_struct{
 
    #parses the output from RNAfoldand reads it into hashes

    my($file) = @_;
    my($id,$desc,$seq,$struct,$mfe) = ();

    open (FILE_STRUCT, "<$file") or die "can not open $file\n";
    while (<FILE_STRUCT>)
    {
        chomp;
        if (/^>(\S+)\s*(.*)/)
	{
	    $id          = $1;
	    $desc        = $2;
	    $seq         = "";
	    $struct      = "";
	    $mfe         = "";
	    while (<FILE_STRUCT>){
                chomp;
                if (/^>(\S+)\s*(.*)/){
		    $hash_desc{$id}   = $desc;
		    $hash_seq{$id}    = $seq;
		    $hash_struct{$id} = $struct;
		    $hash_mfe{$id}    = $mfe;

		    $id          = $1;
		    $desc        = $2;
		    $seq         = "";
		    $struct      = "";
		    $mfe         = "";

		    next;
                }
		if(/^\w/){
		    tr/uU/tT/;
		    $seq .= $_;
		}if(/((\.|\(|\))+)/){
		    $struct .=$1;
		}
		if(/\((\s*-\d+\.\d+)\)/){
		    $mfe = $1;
		}
	    
	    }
        }
    }

    $hash_desc{$id}        = $desc;
    $hash_seq{$id}         = $seq;
    $hash_struct{$id}      = $struct;
    $hash_mfe{$id}         = $mfe;

    close FILE_STRUCT;
    return;
}


sub filter_s{

    #this filtering message is appended to the end of the string of filtering messages outputted for the potential precursor

    my $message=shift;
    $message_filter.=$message."\n";
    return;
}


sub filter_p{

    #this filtering message is appended to the beginning of the string of filtering messages outputted for the potential precursor

    my $message=shift;
    if(defined $message_filter){$message_filter=$message."\n".$message_filter;}
    else{$message_filter=$message."\n";}
    return;
}

    
sub find_freq{

    #finds the frequency of a given read query from its id.

    my($query)=@_;

    if($query=~/x(\d+)/){
	my $freq=$1;
	return $freq;
    }else{
	print STDERR "Problem with read format\n";
	return 0;
    }
}


sub print_hash_comp{

    #prints the 'comp' hash

    my @keys=sort keys %hash_comp;
    foreach my $key(@keys){
	my $value=$hash_comp{$key};
	print "$key  \t$value\n";
    }
}



sub find_strand{

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


sub contained{

    #Is the stretch defined by the first positions contained in the stretch defined by the second?

    my($beg1,$end1,$beg2,$end2)=@_;

    testbeginend($beg1,$end1,$beg2,$end2);

    if($beg2<=$beg1 and $end1<=$end2){
	return 1;
    }else{
	return 0;
    }
}


sub testbeginend{

    #Are the beginposition numerically smaller than the endposition for each pair?

    my($begin1,$end1,$begin2,$end2)=@_;

    unless($begin1<=$end1 and $begin2<=$end2){
	print STDERR "beg can not be larger than end for $subject_old\n";
	exit;
    }
}


sub rev_pos{

#   The blast_parsed format always uses positions that are relative to the 5' of the given strand
#   This means that for a sequence of length n, the first nucleotide on the minus strand base pairs with
#   the n't nucleotide on the plus strand

#   This subroutine reverses the begin and end positions of positions of the minus strand so that they
#   are relative to the 5' end of the plus strand	
   
    my($beg,$end,$lng)=@_;
    
    my $new_end=$lng-$beg+1;
    my $new_beg=$lng-$end+1;
    
    return($new_beg,$new_end);
}

sub rev{

    #reverses the order of nucleotides in a sequence

    my($sequence)=@_;

    my $rev=reverse $sequence;   

    return $rev;
}

sub com{

    #the complementary of a sequence

    my($sequence)=@_;

    $sequence=~tr/acgtuACGTU/TGCAATGCAA/;   
 
    return $sequence;
}

sub revcom{
    
    #reverse complement

    my($sequence)=@_;

    my $revcom=rev(com($sequence));

    return $revcom;
}
