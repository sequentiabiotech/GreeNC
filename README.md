# GreeNC

This repository includes the scripts used to identify the lncRNAs from GreeNC, a Wiki-database
of plant lncRNAs.

![Pipelines](Pipeline.png?raw=true "Pipelines")

## Script 1: identification of lncRNA

A bash script was written in order to identify lncRNAs among a pool of transcripts filtering them
by length below 200nt using the scripting language and filtering by ORF smaller than 120aa with 
Ugene. On one side, BLASTX is called against SwissProt with parameters -evalue 0.001, -outfmt 6, 
and -max_target_seqs 1. The strand is set to "plus" of "both" depending on the group of transcripts 
being analyzed. On the other side, a modification of CPC (0.9-r2) is called with FrameFinder 
parameter -r set to "True" or "False" and BLASTX parameter -S set to "3" or "1" depending on the 
group of transcripts being analyzed.

```
bash script1/lncRNA_strand.sh --fasta FASTA --output OUTDIR --swissprot BLASTDB --ref GTF [optional]
```

## Script 2: discrimination of other ncRNA

A perl script was written in order to discriminate other non-coding transcripts from long non-coding 
transcripts and identify miRNA precursor long non-coding transcripts. On one side, transcripts are 
analyzed by cmscan against Rfam database. On the other side, BLASTn is called with parameters -dust no,
-evalue 0.05, and -word_size 7 against a database of plant mature miRNA sequences from miRBase and 
results are then validated by MIReNA (2.0) with parameters --valid, -y, --mfei -0.69, --amfe -32, 
--ratiomin 0.83, and --ratiomax 1.17. Finally, MIReNA is called again with parameters --valid, -x, 
--mfei -0.69, --amfe -32, --ratiomin 0.83, and --ratiomax 1.17.

```
perl script2/cmscan.pl --in FASTA --out OUTDIR --cmpath RFAMDB --blastdb BLASTDB [optional]
```
