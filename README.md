# Axminster7A

## Using a IWGSC reference gene matching my candidate contig to to pull out additional matches from my RenSeq WT assembly

if you have a fragmented assembly and yout primary candidate is only a partial NLR, the rest of your missing sequence might still be floating around as a seperate contig that didn't get joined due to poor read coverage over the gap. If a homolog or ortholog exists in a reference assembly (such as IWGSC RefSeq v1.0), you may be able to use that to fish for additional contigs that may form part of your gene candidate.

1. using bedtools/2.29.2, the IWGSC v1.0 reference genome and accompanying RefSeq gene annotation file, extract all the annotated gene sequences to a fasta file (Make sure the gene identifier is in the 3rd column of the gff file):

`bedtools getfasta -fi RefSeqv1_chromosomes.fasta -bed RefSeqv1_genes.gff -nameOnly > RefSeqv1_genes.fasta`

2. using blast+/2.9.0, BLAST your candidate to the RefSeqv1 genes:
```
makeblastdb -in RefSeqv1_genes.fasta -dbtype nucl
blastn -query candidate1.fasta -db RefSeqv1_genes.fasta -outfmt 6 -out cand1_vs_RefSeq-genes.blastn.txt
```
>from the output table, see your best hit(s) and whether it could represent a homolog/ortholog (sufficient percent id and alignment length). in this example, the top hit was to gene TraesCS7D01G540500.

3. using samtools/1.9.0, retrieve TraesCS7D01G540500 sequence from RefSeqv1_genes.fasta and BLAST it back to the WT assembly:
```
samtools faidx RefSeqv1_genes.fasta TraesCS7D01G540500 > TraesCS7D01G540500.fasta
makeblastdb -in WT_contigs.fa -dbtype nucl
blastn -query TraesCS7D01G540500.fa -db WT_contigs.fa -outfmt 6 -out homolog_vs_WTcontigs.blastn.txt
```
>from the output table you will see that your initial candidate is probably the top hit. Any adjoining contigs should align adjacent to it on the TraesCS7D01G540500 sequence and not overlap. Columns 7 and 8 can indicate coordinates for this. Also, filtering the BLAST table can narrow down the matches:
`awk '$3 >= 85.0 && $4 >= 500' homolog_vs_WTcontigs.blastn.txt | sort -rnk12 > homolog_vs_WTcontigs.filtered.txt`
>here, we select only matches with a percent id (col 3) of 85 or higher, and an alignment length (col 4) of 500 or more, then sorted by bscore (col 12) so best hits appear at the top.

4. From the shortlist of matching contigs, visually inspect their corresponding RenSeq alignments to see evidence of additional mutations or NLR domains that may complement the initial candidate contig.

It should be noted that this approach compared nucleotides (blastn), which was good enough. However, comparing translated nucleotides (tblastn) may be a more robust approach should initial attempts return too many ambiguous matches.

## Normalising BLAST bit scores and ordering alignments of IWGSC RefSeq v1.0 chr7A genes vs Axminster chr7A de novo assembled contigs

1. Use bedtools/2.26.0 and chr7A gene annotation file to extract IWGSC RefSeq v1.0 chr7A gene sequences:

`bedtools getfasta -fi RefSeqv1_chromosomes.fasta -bed RefSeqv1_chr7A-genes.gff > RefSeqv1_chr7A-genes.fasta`
>the fasta headers will be in the form "chr7A:xxxxx-xxxxx" retaining coordinate information of gene on chr7A

2. BLAST gene multifasta to Axminster 7A assembly:
```
makeblastdb -in Axminster7A_contigs.fasta -dbtype nucl
blastn -num_threads 8 -query RefSeqv1_chr7A-genes.fasta -db Axminster7A_contigs.fasta -outfmt 6 -out 7A-genes_vs_Ax7A-contigs.blastn.txt
```

3. Filter BLAST table by alignment length >=1000, sort by bit score, get top hit per query sequence using blast_filterV2 (TC-Hewitt/Misc_NGS):

`python blast_filterV2.pyc -i 7A-genes_vs_Ax7A-contigs.blastn.txt -o 7A-genes_vs_Ax7A-contigs.filtered.txt -a 1000 -sort1 bscore -topq 1`

4. Get bscore/kb for filtered hits and reorder by query gene position along RefSeq v1.0 chr7A using hspnormalize.py:

`python hspnormalize.py -i 7A-genes_vs_Ax7A-contigs.filtered.txt -o 7A-genes_vs_Ax7A-contigs_scored.ordered.txt`
>coordinate information from fasta headers of query sequences in 1st column of BLAST table is used for ordering, and the tabulated output has the following fields: query gene start, query alignment start, bit score, bit score/kb, subject seq ID

5. columns 1 (x-axis) and 4 (y-axis) can be plotted to see change in BLAST strength along the reference chromosome 

## Filter and sort MashMap output

1. Extract alignments >= 91% identity:

`awk '$10 >= 95' mashmap.out > mashmap.pi91.out`

2. Append column to output giving the mapping length (difference of start and end coordinates):

`awk '{ $11 = $4 - $3 } 1' mashmap.pi91.out | sed 's/ /\t/g' > mashmap.pi91.addlen.out`

3. Sort by appended column and retrieve best alignment per query based on mapping length using sortget.py:

`python sortget.py -i mashmap.pi91.addlen.out -sf 10 -nf 0 -r 1 > mashmap.pi91.addlen.topq1.out`

## Retrieve contigs mapping distal to a putative breakpoint based on MashMap output

From the dotplot of previous MashMap output, an approximate breakpoint could be inferred around position 728000000 on chromosome 7A...

1. extract only those alignments to Chr7A distal to the breakpoint:

`awk '$8 >= 728000000' mashmap.out | grep Chr7A > mashmap_Chr7A.dist728mb.out`

2. retrieve those contigs mapping after the breakpoint using get_contigs.py (TC-Hewitt/Misc_NGS):

`python get_contigs.py -i all_contigs.fasta -o distal_mapping_contigs.fasta -t mashmap_Chr7A.dist728mb.out -s contig_`
