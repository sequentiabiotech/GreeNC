#!/usr/bin/perl

use warnings;
use Getopt::Long;
use Bio::SeqIO;

if (!@ARGV) {
	print "\nNo flags or arguments found.\n";
	print "See 'script-pl --help' 'script.pl -h' for further information.\n\n";
	die;
}

$out = "./";
GetOptions ("in=s" => \$in,
			"reference=s" => \$ref,			#Final reference from cufflinks // merged
			"out=s" => \$out,
			"help" => \$help)
			or die("Error in command line arguments\n");

$helptext = "\n USAGE:\t'script.pl\n\t --in [FASTA]\n\t --reference [GTF]\n\t --out [DIRECTORY]'\n\n The script takes a GTF and a set of sequences in fasta format and outputs three files according to its strand (+/-/unknown).\n\n";

if ($help) {
	print $helptext;
	die;
}

open REFERENCE, "<", $ref or die "Unable to open $ref\n";
while ($line = <REFERENCE>) {
	chomp $line;
	if ($line =~ /transcript_id "([A-Za-z0-9:._\-]*)"/) {
		$name = "SEQUENTIA_".$1;
		@chompline = split(/\t/, $line);
		if ($chompline[6] eq "+") {
			next if $name ~~ @positive;
			push @positive, $name;
			next;
		}
		elsif ($chompline[6] eq "-") {
			next if $name ~~ @negative;
			push @negative, $name;
			next;
		}
		elsif ($chompline[6] eq ".") {
			next if $name ~~ @unknown;
			push @unknown, $name;
			next;
		}
		else {
			print "Could not recognize the strand for ".$1."\n";
			next;
		}
	}
}
close REFERENCE;

sub out {
	($outpath, $array) = @_;
	@array = @{$array};
	open $outfile, '>', $outpath;
	$seqin = Bio::SeqIO->new(-file => $in, -format => "Fasta");
	$seqout = Bio::SeqIO->new(-fh => $outfile, -format => "Fasta");
	while ($seq = $seqin->next_seq()) {
		if ($seq->id ~~ @array) {
			$seqout->width($seq->length);
			$seqout->write_seq($seq);
		}
	}
	close $outfile;
	$seqin = ""; $seqout = "";
}

$par1 = $out."/00-positive.fa";
$par2 = $out."/00-negative.fa";
$par3 = $out."/00-unknown.fa";

out($par1, \@positive);
out($par2, \@negative);
out($par3, \@unknown);


