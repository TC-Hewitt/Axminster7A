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
