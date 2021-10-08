#!/bin/bash

# ./lncRNA_NoDB.sh --fasta ./dataset.fa --output outdir --cpc ./run_predict.sh

clear

#####-Default databases-#####
WORKDIR=$(dirname $(readlink -e $0))

# rRNAdb=/media/sequentia/NAS/src/006-ncRNA/datasets/Generic/rRNA_20131128.fa
#protdb=/media/sequentia/NAS/src/006-ncRNA/datasets/Generic/uniprot_sprot_20131128.fasta
protdb=/media/sequentia/sharing/src/006-ncRNA/datasets/Generic/Uniprot_SwissProt_0421.fasta

echo -e "\n  +--------------------------------+"
echo -e "  |    lncRNA Detection Analysis   |"
echo -e "  +--------------------------------+\n"

echo -e "******* Analysis started at $(date) *******\n"

###############################################FUNCTIONS###############################################

#####-Function to check if all dependencies are met-#####
function deps {
	echo " - Testing dependencies..."
	script=$(dirname $0)
	DEPENDENCIES="ugene blastn blastx fasta_formatter samtools"
	deps_ok=YES
	for dep in $DEPENDENCIES
	do
		if ! which $dep &>/dev/null;  then
			echo This script requires $dep to run but it is not installed.                    
			deps_ok=NO
		fi
	done
	if [[ $deps_ok == NO ]]; then
		echo Unmet dependencies... Aborting!
		exit 1
	else
		echo -e "     ===> Dependencies OK!\n"
		return 0
	fi
}


###############################################PARAMETERS###############################################

#####-Parameter capture-#####
if [ "$#" -eq 0 ] ; then echo -e "\nNo files or parameters found. See $0 --help for command options"
	exit 1
fi
if [ "$#" -gt 0 ]
	then
	for (( i=1; i<="$#"; i++ )) ; do
		case ${!i} in
			--fasta)
				((i+=1))
				fasta=$(readlink -e ${!i})
				fastaname=$(basename $fasta)
				;;
			--output)
				((i+=1))
				if [ ! -d ${!i} ] ; then mkdir -p ${!i} ; echo -e " - Output directory set: ${!i}\n" ; else echo -e " - Output directory set: ${!i}\n" ; fi
				out=$(readlink -e ${!i})
				;;
			--cpc)
				((i+=1))
				cpc_run=$(readlink -e ${!i})
				;;
			--verbose)
				set -exv
				;;
			--ref)
				((i+=1))
				ref=$(readlink -e ${!i})
				;;
			--nodep)
				dependencies=0
				;;
			--rRNA)
				((i+=1))
				rRNAdb=$(readlink -e ${!i})
				rRNA_analysis="1"
				;;
			--swissprot)
                                ((i+=1))
                                protdb=$(readlink -e ${!i})
                                ;;
			--help)	
				echo " 			006-rRNA Detection Script - MANUAL

				SYNOPSIS
					
					./lncRNA_NoDB.sh --fasta ./dataset.fa --output outdir --cpc /../run_predict.sh
						-- rRNA [OPTIONAL]

				DESCRIPTION

					lncRNA detection script allows an accurate and quick detection of lncRNAs of raw 
					  sequence data. First it performs three independent filtering processes in order
					  to distinguish protein-coding RNAs from noncoding RNAs.
						1) 1st Filter: Ugene (ORF searcher) and size (< 200bp).
						2) 2nd Filter: BLAST-UNIPROT.
						3) 3rd Filter: Coding Potential Calculator (CPC) tool.
					
					Finally, with --rRNA option, it discards rRNAs from lncRNAs blasting rRNA DB.
				
				OUTPUT
					- lncRNA_results directory with several files containing coding sequences and
					   noncoding sequences (if --rRNA, then also a rRNA coding sequences file is
					   provided).
					- ncRNA_tmp directory with different temporal files.
					- Locuslist.txt file with input sequences.

				COMMANDS:
					--fasta [path]		# Path to experiment fasta
					--output [name] 	# For folder and files, do not provide a path
					--cpc [path]		# Path to run_predict.sh script, included in CPC directory
					--swissprot [path]	# Path to SwissProt db if it is not the default
					--ref [GTF]		# Reference (if you want to run the script as strand-specific mode)
				  Optional commands:
					--rRNA	[path]		# Path to rRNA db fasta file
					--nodep			# Do not check dependencies					
					--help			# Inception?
					--verbose		# DON'T do that! You're being warned.
					"	
				exit 1
				;;
							
				*)
				echo "\"${!i}\" is not a valid parameter. See --help"
				exit 1
				;;
		esac
	done
fi

if [ -z $dependencies ] ; then deps ; fi 	# To check the dependencies
if [[ -z $fasta ]] || [[ -z $out ]] ; then echo -e "Missing arguments: See --help"  ; exit 1 ; fi
if [[ -z $cpc_run ]] ; then cpc_run=$(readlink -e "$WORKDIR"/cpc-0.9-r2-modified/bin/run_predict.sh) ; fi

###############################################SCRIPT###############################################

ORFpassed=0 ; ORFdump=0 ; count=0 ; lendump=0 ; filt=0 ; dump=0 ; filt07=0 ; dump07=0 ; cpc=1

##### Script ######

#------------Removing old files------------#

if [ -d $out/ncRNA-tmp ]; then rm -R $out/ncRNA-tmp ; fi
if [ -d $out/lncRNA_results ]; then rm -R $out/lncRNA_results ; fi
mkdir -p $out/ncRNA-tmp
mkdir -p $out/lncRNA_results
echo -ne "" > "$out/rRNA_precursors.txt"

cat $fasta | sed 's/>/>SEQUENTIA_/g' | sed 's/ .*$//g' > "$out/ncRNA-tmp/00-Reads.fa" 
samtools faidx "$out/ncRNA-tmp/00-Reads.fa"
fasta=$(readlink -e "$out/ncRNA-tmp/00-Reads.fa")

echo -e "  #-------- 00-Starting analysis --------#"

if [[ $ref ]]; then
	echo -e "      Looking for stranded transcripts..."
	strand="perl $WORKDIR/strandspecific2.pl --in $fasta --reference $ref --out $out/ncRNA-tmp/" ; $strand

	positive=$(readlink -e "$out/ncRNA-tmp/00-positive.fa")
	negative=$(readlink -e "$out/ncRNA-tmp/00-negative.fa")
	unknown=$(readlink -e "$out/ncRNA-tmp/00-unknown.fa")

	cat "$out/ncRNA-tmp/00-positive.fa" "$out/ncRNA-tmp/00-negative.fa" > "$out/ncRNA-tmp/00-stranded.fa"
	stranded=$(readlink -e "$out/ncRNA-tmp/00-stranded.fa")
fi

#------------Fasta formatting------------#

echo -e "      Formatting FASTA to TAB ...\n"
fasta_formatter -t -i "$fasta" -o "$out/ncRNA-tmp/$fastaname.tab"

#if [[ $ref ]]; then
#	fasta_formatter -t -i "$stranded" -o "$out/ncRNA-tmp/stranded.tab"
#	fasta_formatter -t -i "$unknown" -o "$out/ncRNA-tmp/unknown.tab"
#fi

#------------Data extraction------------# 

typeset -A array
while read line 
	do text=( $line ) ; primer=${text[0]}; segon=${#text[1]} ; array[${primer}]=${segon}; 
done < "$out/ncRNA-tmp/$fastaname.tab"

#------------01-Locus list and Blacklist creation------------#

echo "  #-------- 01-Locus list and Blacklist creation --------# "
echo -e "      Filtering with Ugene (ORF searcher) and filtering by size (< 200bp) ...\n"

grep ">" "$fasta" | sed s/">"// | sed s/" "/\\t/ | cut -f 1 | uniq > "$out/ncRNA-tmp/Locuslist.txt"

if [[ $ref ]]; then
	grep ">" "$stranded" | sed s/">"// | sed s/" "/\\t/ | cut -f 1 | uniq > "$out/ncRNA-tmp/Locuslist-stranded.txt"
	grep ">" "$unknown" | sed s/">"// | sed s/" "/\\t/ | cut -f 1 | uniq > "$out/ncRNA-tmp/Locuslist-unknown.txt"
fi

#####-Ugene-######

ugene find-orfs --in="$fasta" --out="$out/ncRNA-tmp/01-$fastaname"_ORFs.txt --min-length=360 --tmp-dir="$out"/ncRNA-tmp	# Find ORFs in both strands

if [[ $ref ]]; then
	## forward ORFs
	grep -v "/dna_len" "$out/ncRNA-tmp/01-$fastaname"_ORFs.txt | grep -v "/protein_len" | grep -v "/ugene_name=" | grep -v "misc_feature    complement" | grep -B 2 "misc_feature" | grep "SEQUENTIA" | tr -d " " > "$out/ncRNA-tmp/ORFforward.txt"
	grep -w -F -f "$out/ncRNA-tmp/ORFforward.txt" "$out/ncRNA-tmp/Locuslist-stranded.txt" > "$out/ncRNA-tmp/Blacklist.txt"
	## Both strands for unstranded
	grep "FEATURES             Location/Qualifiers" -B 1 "$out/ncRNA-tmp/01-$fastaname"_ORFs.txt | grep "SEQUENTIA" | tr -d " " > "$out/ncRNA-tmp/ORF.txt"
	grep -w -F -f "$out/ncRNA-tmp/ORF.txt" "$out/ncRNA-tmp/Locuslist-unknown.txt" >> "$out/ncRNA-tmp/Blacklist.txt"
else
	grep "FEATURES             Location/Qualifiers" -B 1 "$out/ncRNA-tmp/01-$fastaname"_ORFs.txt | grep "SEQUENTIA" | tr -d " " > "$out/ncRNA-tmp/Blacklist.txt"
fi

cat "$out/ncRNA-tmp/Blacklist.txt" | sed 's/SEQUENTIA_//g' > "$out/ORF-detected.txt"	# ORFs detected in the fasta input file

for i in $(cat "$out/ncRNA-tmp/Locuslist.txt"); do
	((count+=1))
	if ! grep -wilq "$i" "$out/ncRNA-tmp/Blacklist.txt" ; then		
		if [ ${array[${i}]} -gt 200 ] ; then	
			#echo "$i greater than 200bp"
			samtools faidx "$fasta" "$i" >> "$out/ncRNA-tmp/01-ORF-filtered.fa"
			echo "$i" >> "$out/ncRNA-tmp/01-ORF-Passed-Locuslist.txt"
			((ORFpassed+=1))
			if [[ $ref ]]; then
				if grep -wilq "$i" "$out/ncRNA-tmp/Locuslist-stranded.txt" ; then	
					samtools faidx "$fasta" "$i" >> "$out/ncRNA-tmp/01-ORF-filtered-stranded.fa"
					echo "$i" >> "$out/ncRNA-tmp/01-ORF-Passed-Locuslist-stranded.txt"
				elif grep -wilq "$i" "$out/ncRNA-tmp/Locuslist-unknown.txt" ; then	
					samtools faidx "$fasta" "$i" >> "$out/ncRNA-tmp/01-ORF-filtered-unknown.fa"
					echo "$i" >> "$out/ncRNA-tmp/01-ORF-Passed-Locuslist-unknown.txt"
				fi
			fi
		else
			#echo "$i smaller than 200bp"
			echo "$i" >> "$out/ncRNA-tmp/01-Length-Dumped.txt"
			echo "$i" >> "$out/ncRNA-tmp/Blacklist.txt"
			((lendump+=1))
		fi
	else
		((ORFdump+=1))	
	fi
done

#------------02-Blast------------# 

echo "  #-------- 02-Blast UNIPROT --------# "
echo -e "      Filtering RNA sequences with similar UNIPROT hits...\n"

if [[ $ref ]]; then
	cmd1="blastx -strand plus -evalue 0.001 -query $out/ncRNA-tmp/01-ORF-filtered-stranded.fa -db $protdb -out $out/ncRNA-tmp/02-Blast-out-stranded.txt -num_threads 16 -outfmt 6  -max_target_seqs 1"
	cmd3="blastx -strand both -evalue 0.001 -query $out/ncRNA-tmp/01-ORF-filtered-unknown.fa -db $protdb -out $out/ncRNA-tmp/02-Blast-out-both.txt -num_threads 16 -outfmt 6  -max_target_seqs 1"
	$cmd1 ; $cmd3
	cat "$out/ncRNA-tmp/02-Blast-out-stranded.txt" "$out/ncRNA-tmp/02-Blast-out-both.txt" > "$out/ncRNA-tmp/02-Blast-out.txt"
else
	cmd="blastx -evalue 0.001 -query $out/ncRNA-tmp/01-ORF-filtered.fa -db $protdb -out $out/ncRNA-tmp/02-Blast-out.txt -num_threads 16 -outfmt 6  -max_target_seqs 1"
	$cmd
fi

cut -f 1 "$out/ncRNA-tmp/02-Blast-out.txt" | sort -u | while read -ra line ; do 	
		if ! grep -wilq "${line[0]}" "$out/ncRNA-tmp/Blacklist.txt" ; then	
			echo ${line[0]} >> "$out/ncRNA-tmp/02-blast-coding.txt"		
		else
			echo ${line[0]} >> "$out/ncRNA-tmp/02-blast-noncoding.txt"	
		fi
done

for a in $(cat "$out/ncRNA-tmp/01-ORF-Passed-Locuslist.txt"); do
 	if ! grep -wlq "$a" "$out/ncRNA-tmp/02-blast-coding.txt" ; then 	
		echo "$a" >> "$out/ncRNA-tmp/02-Blast-Filtered.txt"				
 	fi
 done

dump02=$(cat "$out/ncRNA-tmp/02-blast-coding.txt" | wc -l)

#------------CPC------------#

if [ $cpc ] ; then
	echo -e "  #-------- Running CPC --------#\n"
	if [[ $ref ]]; then
		cmd1="bash $cpc_run $out/ncRNA-tmp/01-ORF-filtered-stranded.fa $out/ncRNA-tmp/CPC-out-stranded $out/ncRNA-tmp/ $out/ncRNA-tmp/CPC-out-stranded2"
		cmd2="bash $cpc_run $out/ncRNA-tmp/01-ORF-filtered-unknown.fa $out/ncRNA-tmp/CPC-out-nostranded $out/ncRNA-tmp/ $out/ncRNA-tmp/CPC-out-nostranded2 OK"
		$cmd1 ; $cmd2
		cat "$out/ncRNA-tmp/CPC-out-stranded" "$out/ncRNA-tmp/CPC-out-nostranded" > "$out/ncRNA-tmp/CPC-out"
		cat "$out/ncRNA-tmp/CPC-out-stranded2.homo" "$out/ncRNA-tmp/CPC-out-nostranded2.homo" > "$out/ncRNA-tmp/CPC-out2.homo"
		cat "$out/ncRNA-tmp/CPC-out-stranded2.orf" "$out/ncRNA-tmp/CPC-out-nostranded2.orf" > "$out/ncRNA-tmp/CPC-out2.orf"
	else
		cmd="bash $cpc_run $out/ncRNA-tmp/01-ORF-filtered.fa $out/ncRNA-tmp/CPC-out $out/ncRNA-tmp/ $out/ncRNA-tmp/CPC-out2"
		$cmd
	fi	
	grep "noncoding" "$out/ncRNA-tmp/CPC-out" | cut -f 1 | sort -u > "$out/ncRNA-tmp/CPC-filter.txt"
	grep -w "coding" "$out/ncRNA-tmp/CPC-out" | cut -f 1 | sort -u >> "$out/ncRNA-tmp/02-CPC-coding.txt" 
	count_cpc=$(cat "$out/ncRNA-tmp/02-CPC-coding.txt" | wc -l)
	sort -u "$out/ncRNA-tmp/02-Blast-Filtered.txt" "$out/ncRNA-tmp/CPC-filter.txt" >  "$out/ncRNA-tmp/02-dumb-list.txt"		
	for a in $(cat "$out/ncRNA-tmp/02-dumb-list.txt"); do
		if ! grep -wlq "$a" "$out/ncRNA-tmp/Blacklist.txt" ; then
			echo "$a" >> "$out/ncRNA-tmp/02-Blast-Passed-Locuslist.txt" 	# Deberían ser os mismos que 02-dumb-list.txt
		else
			echo "$a" >> "$out/ncRNA-tmp/02-Blacklist.txt"
		fi
	done
	count_apte=$(cat "$out/ncRNA-tmp/02-Blast-Passed-Locuslist.txt" | wc -l)
	count_discard=$(($count - $count_apte)) 
	cat "$out/ncRNA-tmp/02-Blast-Passed-Locuslist.txt" | while read -ra line ; do
		samtools faidx "$fasta" "$line" >> "$out/ncRNA-tmp/02-Blast-filtered.fa"	# Todos los RNAs de interés pasan a este archivo fasta
	done
fi

grep ">" "$out/ncRNA-tmp/02-Blast-filtered.fa" | sed s/">"// | sort -u > "$out/ncRNA-tmp/02-Blast-Filtered.txt"		

sort "$out/ncRNA-tmp/02-Blast-Filtered.txt" "$out/ncRNA-tmp/Locuslist.txt" | uniq -u >> "$out/ncRNA-tmp/Blacklist.txt"


#####-"blasting function" call-#####
last_fasta="$out/ncRNA-tmp/02-Blast-filtered.fa"
last_locuslist="$out/ncRNA-tmp/02-Blast-Passed-Locuslist.txt"
counter=2

cat "$out/ncRNA-tmp/02-Blast-Passed-Locuslist.txt" | sed 's/SEQUENTIA_//g' > "$out/02-Blast-Passed-Locuslist.txt"

grep -w "noncoding" "$out/ncRNA-tmp/CPC-out" | awk '{print $1}' | grep -w -v -f "$out/ncRNA-tmp/02-blast-coding.txt" > "$out/lncRNA_results/high_confidence.txt"

#-----------04-Blasting rRNA---------#
if [ $rRNA_analysis ] ; then

	((counter+=1))
	echo -e "  #--------0$counter-Blasting against rRNA db --------#\n"
	cmd="blastn -task blastn -dust no -perc_identity 100 -evalue 0.01 -word_size 7 -query $last_fasta -db $rRNAdb -out $out/ncRNA-tmp/0$counter-Blast-rRNA-out.txt -num_threads 16 -outfmt 6 -num_alignments 1"
	$cmd 
	cat "$out/ncRNA-tmp/0$counter-Blast-rRNA-out.txt" | sort -u | while read -ra line ; do
		if ! grep -wlq "${line[0]}" "$out/ncRNA-tmp/Blacklist.txt" ; then
			echo ${line[0]} >> "$out/ncRNA-tmp/Blacklist.txt" 		# Se meten los rRNAs en la blacklist para ir descartando hits y seguir con el análisis de los demás ncRNAs...
			echo ${line[0]} | sed 's/SEQUENTIA_//g' >> "$out/rRNA_precursors.txt"	# RNAs RIBOSOMALES
		fi
	done
	for y in $(cat "$last_locuslist"); do
		if ! grep -wlq "$y" "$out/ncRNA-tmp/Blacklist.txt" ; then
			samtools faidx "$fasta" "$y" >> "$out/ncRNA-tmp/0$counter-Blast-rRNA-filtered.fa"
			((filt07+=1))
			else
			samtools faidx "$fasta" "$y" >> "$out/ncRNA-tmp/0$counter-Blast-rRNA-detected.fa"
			((dump07+=1))
		fi
	done

	cat "$out/ncRNA-tmp/0$counter-Blast-rRNA-detected.fa" | sed 's/SEQUENTIA_//g' > "$out/rRNA_precursors.fa"
	numrRNA=$(cat "$out/rRNA_precursors.txt" | wc -l)
	lncRNA=$(($count_apte - $numrRNA))
	
	
	for i in `cat "$out/02-Blast-Passed-Locuslist.txt"`; do echo $i | grep -x -v -f "$out/rRNA_precursors.txt" ; done > "$out/lncRNA_precursors.txt"  
	cat "$out/lncRNA_precursors.txt" | sed 's/^/SEQUENTIA_/g' | sed 's/ .*$//g' > "$out/lncRNA_precursors1.txt"
	cat "$out/lncRNA_precursors1.txt" | while read -ra line ; do
		samtools faidx "$fasta" "$line" >> "$out/lncRNA_precursors.fa"

	done

	rm "$out/lncRNA_precursors1.txt"

fi

count=$(wc -l "$out/ncRNA-tmp/Locuslist.txt" | tr " " "\t" | cut -f 1)		# Number of input sequences
# numrRNA=$(cat "$out/rRNA_precursors.txt" | wc -l)
# lncRNA=$(($count_apte - $numrRNA))


# cat "$out/ncRNA-tmp/0$counter-Blast-rRNA-detected.fa" | sed 's/SEQUENTIA_//g' > "$out/rRNA_precursors.fa"
# cat "$out/ncRNA-tmp/02-Blast-Passed-Locuslist.txt" | sed 's/SEQUENTIA_//g' > "$out/02-Blast-Passed-Locuslist.txt"
cat "$out/ncRNA-tmp/Locuslist.txt" | sed 's/SEQUENTIA_//g' > "$out/Locuslist.txt"

# for i in `cat "$out/02-Blast-Passed-Locuslist.txt"`; do echo $i | grep -x -v -f "$out/rRNA_precursors.txt" ; done > "$out/lncRNA_precursors.txt"  
# cat "$out/lncRNA_precursors.txt" | sed 's/^/SEQUENTIA_/g' | sed 's/ .*$//g' > "$out/lncRNA_precursors1.txt"
# cat "$out/lncRNA_precursors1.txt" | while read -ra line ; do
# 	samtools faidx "$fasta" "$line" >> "$out/lncRNA_precursors.fa"
# done

grep -w -v -F -f "$out/02-Blast-Passed-Locuslist.txt" "$out/Locuslist.txt" > "$out/Coding_reads.txt"
cat "$out/Coding_reads.txt" | sed 's/^/SEQUENTIA_/g' | sed 's/ .*$//g' > "$out/Coding_reads1.txt"
cat "$out/Coding_reads1.txt" | while read -ra line ; do
	samtools faidx "$fasta" "$line" >> "$out/Coding_reads.fa"
done

rm "$out/Coding_reads1.txt"
# rm "$out/lncRNA_precursors1.txt"
rm "$out/02-Blast-Passed-Locuslist.txt"
mv -t "$out/ncRNA-tmp" "$out/ORF-detected.txt"

if [ $rRNA_analysis ] ; then

	echo -e "\n +-------------------------------------+"
	echo -e " |  lncRNA Detection Analysis RESULTS: |" | tee -a "$out/lncRNA_results/Final_log.txt"
	echo -e " +-------------------------------------+"
	echo -e "   1) FILTERING:" | tee -a "$out/lncRNA_results/Final_log.txt"
	echo -e "       - Total number of reads proccessed: $count" | tee -a "$out/lncRNA_results/Final_log.txt"
	echo -e "       - Reads shorter than 200bp: $lendump" | tee -a "$out/lncRNA_results/Final_log.txt"
	echo -e "       - Reads containing a long ORFs in Ugene: $ORFdump" | tee -a "$out/lncRNA_results/Final_log.txt"
	echo -e "       - Reads with similar UNIPROT hits: $dump02" | tee -a "$out/lncRNA_results/Final_log.txt"
	echo -e "       - Reads with coding potential in CPC: $count_cpc" | tee -a "$out/lncRNA_results/Final_log.txt"
	echo -e "           ====> Total reads discarded: $count_discard" | tee -a "$out/lncRNA_results/Final_log.txt"
	echo -e "           ====> Total potential RNA's: $count_apte" | tee -a "$out/lncRNA_results/Final_log.txt"
	echo -e "    2) lncRNA DETECTION: " | tee -a "$out/lncRNA_results/Final_log.txt"
	echo -e "       - rRNA precursors: $numrRNA" | tee -a "$out/lncRNA_results/Final_log.txt"
	echo -e "       - lncRNA precursors: $lncRNA" | tee -a "$out/lncRNA_results/Final_log.txt"
	echo -e "\n******* Analysis finalised at $(date) *******#\n"
	mv -t "$out/lncRNA_results" "$out/Coding_reads.txt" "$out/Coding_reads.fa" "$out/lncRNA_precursors.txt" "$out/lncRNA_precursors.fa" "$out/rRNA_precursors.txt" "$out/rRNA_precursors.fa"

else

	mv -t "$out/lncRNA_results" "$out/Coding_reads.txt" "$out/Coding_reads.fa"
	mv "$out/ncRNA-tmp/02-Blast-filtered.fa" "$out/lncRNA_results/NonCoding_reads.fa"
	mv "$out/ncRNA-tmp/02-Blast-Passed-Locuslist.txt" "$out/lncRNA_results/NonCoding_reads.txt"
	echo -e "\n +-------------------------------------+"
	echo -e " |  lncRNA Detection Analysis RESULTS: |" | tee -a "$out/lncRNA_results/Final_log.txt"
	echo -e " +-------------------------------------+"
	echo -e "   1) FILTERING:" | tee -a "$out/lncRNA_results/Final_log.txt"
	echo -e "       - Total number of reads proccessed: $count" | tee -a "$out/lncRNA_results/Final_log.txt"
	echo -e "       - Reads shorter than 200bp: $lendump" | tee -a "$out/lncRNA_results/Final_log.txt"
	echo -e "       - Reads containing a long ORFs in Ugene: $ORFdump" | tee -a "$out/lncRNA_results/Final_log.txt"
	echo -e "       - Reads with similar UNIPROT hits: $dump02" | tee -a "$out/lncRNA_results/Final_log.txt"
	echo -e "       - Reads with coding potential in CPC: $count_cpc" | tee -a "$out/lncRNA_results/Final_log.txt"
	echo -e "   2) CLASSIFICATION:" | tee -a "$out/lncRNA_results/Final_log.txt"
	echo -e "       - Total CODING reads: $count_discard" | tee -a "$out/lncRNA_results/Final_log.txt"
	echo -e "       - Total NON CODING RNA's: $count_apte" | tee -a "$out/lncRNA_results/Final_log.txt"
	echo -e "\n******* Analysis finalised at $(date) *******#\n"

fi

