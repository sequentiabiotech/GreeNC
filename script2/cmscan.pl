#!/usr/bin/perl
#Andreu Paytuvi
#Required software in your path: blastn, Infernal 1.0.
#Parameters to set in your computer: path to the blast db, path to the db of CMs, path to MIReNA.sh (program).
#Typical usage: call 'script.pl -h' to know.

use warnings;
use Getopt::Long;
use Bio::SeqIO;

$srcpath = $0;
$srcpath =~ s/cmscan.pl//g;

if (!@ARGV) {
	print "\nNo flags or arguments found.\n";
	print "See 'script-pl --help' 'script.pl -h' for further information.\n\n";
	die;
}

$threads = 2; $out = "./";
$cmpath = $srcpath."Rfam.cm";
$blastdb = $srcpath."mature.fasta";
$mirpath = $srcpath."MIReNA-2.0/MIReNA.sh";
$nomirna = 0;

GetOptions ("threads=i" => \$threads,
            "out=s" => \$out,
            "cmpath=s" => \$cmpath,
            "mirpath=s" => \$mirpath,
            "blastdb=s" => \$blastdb,
            "nomirna+" => \$nomirna,
            "in=s" => \$in,
            "help" => \$help
            )
	or die("Error in command line arguments.\nSee 'script-pl --help' 'script.pl -h' for further information.\n\n");

if ($help) {
	print "\nUSAGE:\t'script.pl\n\t --in [FILE]\n\t --nomirna\t\t#Optional. No pre-miRNAs in your output.\n\t --threads [NUM]\t#Optional. Default 2.\n\t --cmpath [FILE]\t#Optional. Set by default.\n\t --mirpath [FILE]\t#Optional. Set by default.\n\t --blastdb [FILE]\t#Optional. Set by default.\n\t --out [DIRECTORY]'\t#Optional. Default in current directory.\n\n";
	die;
}

mkdir $out;
mkdir $out."/ncRNA_prediction";
mkdir $out."/ncRNA_prediction/tmp";
mkdir $out."/ncRNA_prediction/outputs";

sub Chomp_fasta {
	open INPUT, $in or die "Could not open the file.";
	open $fh, '>', $in.".fasta" or die "Could not open the file.";
	while ($line = <INPUT>) {
        if ($line =~ /^>/) {
		if ($. < 2) {
			print $fh $line;
			next;
		}
		print $fh $sequence."\n";
		chomp $line;
           	print $fh $line."\n";
            	$sequence = "";
            	next;
        }
        else {
            chomp $line;
            $sequence .= $line;
        }
	}
	print $fh $sequence;
	close INPUT;
	close $fh;
	$in = $in.".fasta";
}

sub Rfam {
	print "\n###############################################\n";
	print "###Calling Infernal against CMs of Rfam...\n";
	system("cmscan -E 0.01 --noali --cpu ".$threads." -o ".$out."/ncRNA_prediction/tmp/cmscan.txt ".$cmpath." ".$in);
	mkdir $out."/ncRNA_prediction/outputs/Rfam";

	$lncrna_counter = 0;
	$mirna_counter = 0;
	$other_counter = 0;
	$rrna_counter = 0;
	$trna_counter = 0;

	###families of lncRNAs in Rfam release 11.
	@lncrnas = ("FMR1-AS1_1","FMR1-AS1_2","Yar_2","WT1-AS_3","WT1-AS_2","WT1-AS_1","FTX_1","FTX_2","FTX_3","FTX_4","FTX_5","Yar_1","Six3os1_7","adapt33_4","CDKN2B-AS","CDKN2B-AS_2","CDKN2B-AS_3","ZNRD1-AS1_3","Six3os1_6","UCA1","ZNFX1-AS1_3","ZFAT-AS1_1","JPX_1","JPX_2","GHRLOS","Six3os1_5","KCNQ1DN","KCNQ1OT1_1","KCNQ1OT1_2","KCNQ1OT1_3","KCNQ1OT1_5","ST7-OT4_4","Six3os1_4","Six3os1_3","CLRN1-AS1","ZEB2_AS1_4","GNAS-AS1_1","GNAS-AS1_2","GNAS-AS1_3","GNAS-AS1_4","GNAS-AS1_5","ST7-OT4_3","TUG1_4","TUG1_3","TUG1_2","TUG1_1","Six3os1_2","Six3os1_1","lincRNA-p21_1","lincRNA-p21_2","TTC28-AS1_4","ST7-OT4_2","LOC285194","ST7-OT4_1","TTC28-AS1_3","MALAT1","ST7-OT3_4","ST7-OT3_3","ST7-OT3_2","TTC28-AS1_2","H19_1","H19_2","H19_3","TTC28-AS1_1","MEG3_1","MEG3_2","MEG8_1","MEG8_2","MEG8_3","MESTIT1_1","MESTIT1_2","MESTIT1_3","ZEB2_AS1_3","ZEB2_AS1_2","BC040587","ST7-OT3_1","MIAT_exon1","MIAT_exon5_1","MIAT_exon5_2","MIAT_exon5_3","ST7-AS2_2","HAR1A","ZEB2_AS1_1","Mico1","ST7-AS2_1","DAOA-AS1_1","MIMT1_2","DAOA-AS1_2","ST7-AS1_2","ST7-AS1_1","TP73-AS1","TP53TG1_2","RMST_9","RMST_8","RMST_7","RMST_6","RMST_5","RMST_4","RMST_3","RMST_2","RMST_10","RMST_1","SPRY4-IT1_2","SPRY4-IT1_1","TP53TG1_1","Sphinx_2","Sphinx_1","ZNFX1-AS1_2","SOX2OT_exon4","SOX2OT_exon3","SOX2OT_exon2","SOX2OT_exon1","ZNRD1-AS1_2","ZNFX1-AS1_1","ZFAT-AS1_3","ZFAT-AS1_2","DGCR5","ZNRD1-AS1_1","SMCR2_2","SMCR2_1","HOTAIRM1_1","HOTAIRM1_2","HOTAIRM1_3","RFPL3-AS1_2","RFPL3-AS1_1","HOTAIRM1_4","HOTAIRM1_5","HOTAIR_1","HOTAIR_2","HOTAIR_3","HOTAIR_4","HOTAIR_5","HOTTIP_1","HOTTIP_2","HOTTIP_3","HOTTIP_4","PVT1_7","PVT1_6","PVT1_5","PVT1_4","PVT1_3","PVT1_2","PVT1_1","HOXA11-AS1_1","HOXA11-AS1_2","HOXA11-AS1_3","HOXA11-AS1_4","HOXA11-AS1_5","HOXA11-AS1_6","HOXB13-AS1_1","HOXB13-AS1_2","DISC2","DLEU1_1","PRINS","DLEU1_2","HSR-omega_1","HSR-omega_2","DLEU2_1","DLEU2_2","HTT-AS1_1","HTT-AS1_2","HTT-AS1_3","HULC","DLEU2_3","HYMAI","DLEU2_4","DLEU2_5","DLEU2_6","DLG2-AS1_1","DLG2-AS1_2","VIS1","adapt33_1","PISRT1","Pinc","Xist_exon4","Xist_exon1","bxd_1","SMAD5-AS1_4","bxd_2","SMAD5-AS1_3","PCGEM1","PCA3_2","PCA3_1","PART1_3","PART1_2","PART1_1","bxd_3","EGOT","bxd_4","bxd_5","bxd_6","bxd_7","adapt33_2","adapt33_3","WT1-AS_8","Evf1_1","Evf1_2","SMAD5-AS1_2","FAM13A-AS1_1","FAM13A-AS1_2","FAS-AS1","SMAD5-AS1_1","Vax2os1_3","NPPA-AS1_3","NPPA-AS1_2","NPPA-AS1_1","Vax2os1_2","Nkx2-2as","Vax2os1_1","NEAT1_3","NEAT1_2","NEAT1_1","NCRUPAR_2","NCRUPAR_1","NBR2","NAMA_2","NAMA_1","WT1-AS_7","TCL6_3","TCL6_2","TCL6_1","WT1-AS_6","WT1-AS_5","WT1-AS_4","Yar_3");

	open FILE, $out."/ncRNA_prediction/tmp/cmscan.txt" or die "Could not open the file.";
	while ($lines = <FILE>) {
	    if ($lines =~ /Query:/) {
	    	$lines =~ s/ +/\t/g;
	    	@chomp = split(/\t/, $lines);
	      	$seqname = $chomp[1];
	       	next;
	    }
	    if ($lines =~ /\(1\)/ ) {
	    	$lines =~ s/ +/\t/g;
	       	chop $lines;
	       	$hash{$seqname} = $lines;
	    }
	}
	close FILE;

	for $key (keys %hash) {
		@chompline = split(/\t/, $hash{$key});
		for $lncrna (@lncrnas) {
			if ($lncrna eq $chompline[6]) {
				$hash{$key} = "lncRNA";
				$lncrna_counter++;
				push @lncrna_names, $key;
				last;
			}
		}
		if ($chompline[6] =~ /mir/ or $chompline[6] =~ /MIR/ or $chompline[6] eq "let-7" or $chompline[6] eq "lin-4" or $chompline[6] eq "lsy-6") {
			$hash{$key} = "miRNA";
			$mirna_counter++;
			push @mirna_names, $key;
			$mirnarfam{$key}{family} = $chompline[6];
			$mirnarfam{$key}{evalue} = $chompline[3];
			next;
		}
		elsif ($chompline[6] =~ /tRNA/) {
			$hash{$key} = "tRNA";
			$trna_counter++;
			push @trna_names, $key;
			next;
		}
		elsif ($chompline[6] =~ /rRNA/) {
			$hash{$key} = "rRNA";
			$rrna_counter++;
			push @rrna_names, $key;
			next;
		}
		else {
			$hash{$key} = "Other ncRNA";
			$other_counter++;
			push @other_names, $key;
		}
	}

	$total_counter = $lncrna_counter + $mirna_counter + $rrna_counter + $trna_counter + $other_counter;
	print "\n###############################################\n";
	print "######## Number of ncRNA found in Rfam ########\n";
	print "###############################################\n\n";
	printf "%7s %7s %7s %7s %7s %7s\n", "lncRNAs", "miRNAs", "rRNAs", "tRNAs", "Others", "Total";
	printf "%7s %7s %7s %7s %7s %7s\n", $lncrna_counter, $mirna_counter, $rrna_counter, $trna_counter, $other_counter, $total_counter;

	if (@lncrna_names) {
		open $fh, '>', $out."/ncRNA_prediction/outputs/Rfam/lncRNAs.txt";
		foreach (@lncrna_names) {
			print $fh "$_\n";
		}
		close $fh;
	}
	if (@mirna_names) {
		open $fh, '>', $out."/ncRNA_prediction/outputs/Rfam/miRNAs.txt";
		foreach (@mirna_names) {
			push @discarded, $_;
			print $fh "$_\n";
		}
		close $fh;
	}
	if (@rrna_names) {
		open $fh, '>', $out."/ncRNA_prediction/outputs/Rfam/rRNAs.txt";
		foreach (@rrna_names) {
			push @discarded, $_;
			push @discarded2, $_;
			print $fh "$_\n";
		}
		close $fh;
	}
	if (@trna_names) {
		open $fh, '>', $out."/ncRNA_prediction/outputs/Rfam/tRNAs.txt";
		foreach (@trna_names) {
			push @discarded, $_;
			push @discarded2, $_;
			print $fh "$_\n";
		}
		close $fh;
	}
	if (@other_names) {
		open $fh, '>', $out."/ncRNA_prediction/outputs/Rfam/others.txt";
		foreach (@other_names) {
			push @discarded, $_;
			push @discarded2, $_;
			print $fh "$_\n";
		}
		close $fh;
	}
	open $fh, '>', $out."/ncRNA_prediction/outputs/Rfam/miRNA_families.txt";
	for $key (keys %mirnarfam) {
		print $fh $key."\t".$mirnarfam{$key}{family}."\t".$mirnarfam{$key}{evalue}."\n";
	}
	close $fh;
}

sub Blast {
	print "\n###############################################\n";
	print "###Calling blastn against databse of mature miRNAs of plants (homology search)...\n";
	print "blastn -dust no -evalue 0.05 -word_size 7 -query ".$in." -db ".$blastdb." -num_threads ".$threads." -outfmt 6 -out ".$out."/ncRNA_prediction/tmp/blast1.txt\n";
	system("blastn -dust no -evalue 0.05 -word_size 7 -query ".$in." -db ".$blastdb." -num_threads ".$threads." -outfmt 6 -out ".$out."/ncRNA_prediction/tmp/blast1.txt");
	open CHOMPED, $in;
	while ($line = <CHOMPED>) {
		chomp $line;
		if ($line =~ />/) {
			@chompline = split(/ /, $line);
			$name = $chompline[0];
		}
		else {
			$length{$name} = length($line);
		}
	}
	close CHOMPED;
	mkdir $out."/ncRNA_prediction/outputs/Blast";
	print "###Checking hits with MIReNA...\n";
	open BLAST, $out."/ncRNA_prediction/tmp/blast1.txt";
	open $fh, '>', $out."/ncRNA_prediction/tmp/seq.fa";
	while ($lines = <BLAST>) {
		@chompline = split(/\t/, $lines);
		if ($chompline[3] > 18 and $chompline[4] <= 2) {
			$start = $chompline[6] -1;
			$end = $chompline[7] - 1;
			$after = $length{">".$chompline[0]} - $end;
			$header = $chompline[0]." before:".$start." after:".$after;
			next if defined $hashofheaders{$header};
			$reg = 0; open FAST, $in;
			while ($line = <FAST>) {
				chomp $line;
				if ($line =~ /^>$chompline[0]$/ or $line =~ /^>$chompline[0] /) {
					print $fh ">".$chompline[0]." before:".$start." after:".$after."\n";
					$hashofheaders{$header} = 1;
					$reg = 1;
				}
				elsif ($reg == 1 and $line !~ />/) {
					print $fh $line."\n";
				}
				elsif ($reg == 1 and $line =~ />/) {
					last;
				}
			}
			close FAST;
		}
	}
	close BLAST;
	close $fh;
	system("bash ".$mirpath." -v -y --file ".$out."/ncRNA_prediction/tmp/seq.fa --output ".$out."/ncRNA_prediction/tmp/checkblast.txt");
	open BLASTCHECK, $out."/ncRNA_prediction/tmp/checkblast.txt";
	open $fh, '>', $out."/ncRNA_prediction/outputs/Blast/miRNAs.txt";
	while ($line = <BLASTCHECK>) {
		if ($line =~ />/) {
			$line =~ s/>//g; $line =~ s/ +/\t/g; $line =~ s/before://g;
			@splitlines = split(/\t/, $line);
			$startkey = $splitlines[1] + 1;
			if ($splitlines[0] ~~ @blastkeys) {
				#do nothing
			}
			else {
				print $fh $splitlines[0]."\n";
				push @blastkeys, $splitlines[0];
			}
			next if $startkey ~~ @{$homologyhash{$splitlines[0]}};
			push @{$homologyhash{$splitlines[0]}}, $startkey;
		}
	}
	close BLASTCHECK; close $fh;
	open BLAST4, $out."/ncRNA_prediction/tmp/blast1.txt";
	while ($blastlines = <BLAST4>) {
		@splitblastlines = split(/\t/, $blastlines);
		if ($splitblastlines[6] ~~ @{$homologyhash{$splitblastlines[0]}}) {
			next if $splitblastlines[1] ~~ @{$homologyhits{$splitblastlines[0]}};
			push @{$homologyhits{$splitblastlines[0]}}, $splitblastlines[1];
		}
	}
	close BLAST4;
	open $fi, '>', $out."/ncRNA_prediction/outputs/Blast/homology.txt";
	for $key (keys %homologyhits) {
		print $fi $key."\t"; print $fi join(", ", @{$homologyhits{$key}}); print $fi "\n";
	}
	close $fi;	
	@blasthits = `cat $out/ncRNA_prediction/outputs/Blast/miRNAs.txt`;
	$blastcount = 0;
	for $hit (@blasthits) {
		chomp $hit;
		$blastcount++;
		next if $hit ~~ @discarded;
		push @discarded, $hit;
	}
	print "\nBlast found ".$blastcount." miRNAs.\n";
	$outpath = $out."/ncRNA_prediction/tmp/seq2.fa";
	@blasthits = ();
	out($outpath);
}

sub MIReNA {
	print "\n###############################################\n";
	print "###Finding transcripts that meet miRNA features with MIReNA...\n\n";
	$inpath = $out."/ncRNA_prediction/tmp/seq2.fa";
	mkdir $out."/ncRNA_prediction/outputs/MIReNA";
	system("bash ".$mirpath." -v -x --file ".$inpath." --output ".$out."/ncRNA_prediction/tmp/mirena.txt");
	system("grep '>' ".$out."/ncRNA_prediction/tmp/mirena.txt |  cut -d ' ' -f1| sed 's/^>//g' > ".$out."/ncRNA_prediction/outputs/MIReNA/miRNAs.txt");
	$mirenahits = `cat $out/ncRNA_prediction/outputs/MIReNA/miRNAs.txt | wc -l`; chomp $mirenahits;
	print "MIReNA found ".$mirenahits." putative novel miRNAs.\n";
	$finaloutput = $out."/ncRNA_prediction/lncRNA.fa";
	if ($nomirna == 0) {
		@discarded = ();
		@discarded = @discarded2;
	}
	out($finaloutput);
}

sub out {
	$outpath = shift (@_);
	open $outfile, '>', $outpath;
	$seqin = Bio::SeqIO->new(-file => $in, -format => "Fasta");
	$seqout = Bio::SeqIO->new(-fh => $outfile, -format => "Fasta");
	while ($seq = $seqin->next_seq()) {
		if ($seq->id ~~ @discarded) {
			#do nothing
		}
		else {
			$seqout->width($seq->length);
			$seqout->write_seq($seq);
		}
	}
	close $outfile;
	$seqin = ""; $seqout = "";
}

Chomp_fasta;
Rfam;
Blast;
MIReNA;

$numberdiscarded = $#discarded + 1;

print "\n###############################################\n";
print "Number of transcripts discarded: ".$numberdiscarded."\n";
print "###############################################\n\n";

