![alt text](https://github.com/SAMtoBAM/MUMandCo/blob/master/MUM%26Co.png)

MUM&Co is a simple bash script that uses Whole Genome Alignment information provided by MUMmer (v3 and v4) to detect variants. <br/>

MUM&Co is able to detect: <br/>
Deletions (unique and mobile), insertions (novel and mobile) and tandem duplications (>50bp) <br/>
Inversions and translocations (>1kb)

MUM&Co requires two steps prior to SV detection:

1. Alignment to the reference which will be used to call SVs and renaming based on paired chromosomes <br/>
Contigs corresponding to 2 reference chromosomes should have all reference chromsome names split by an underscore '_' <br/> 
If multiple Contigs correspond to the same chromosome, they can also be split by '_' then given a unique identifier such as '*_contig1' and *_contig2' <br/>
Tools such as RaGOO and Ragout can do this alongside scaffolding of contigs

2. Update bash file with: <br/> 
Paths to MUMmer scripts (can be found using '$ which  nucmer' etc), <br/>
Paths to both reference and query genomes, <br/>
Size of reference genome in base pairs, <br/>
Names of both reference and query <br/>
And desired output prefix <br/>
Optional: Add paths to BLAST scripts for distinguishing novel and mobile INDELS.
