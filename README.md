![alt text](https://github.com/SAMtoBAM/MUMandCo/blob/master/MUM%26Co.png)

MUM&Co is a simple bash script that uses the Whole Genome Alignment information provided by the MUMmer package to detect variants. <br/>
MUM&Co requires two steps prior to SV detection:

1. Alignment to the reference which will be used to call SVs and check renaming based on paired chromosomes <br/>
Contigs corresponding to 2 reference chromosomes should have all reference chromsome names split by an underscore '_' <br/> 
Tools such as RaGOO and Ragout can do this alongside scaffolding of contigs

2. Update bash file with: <br/> 
Paths to MUMmer scripts (can be found using '$ which  nucmer' etc), <br/>
Paths to both reference and query genomes, <br/>
Names of both reference and query <br/>
And desired output prefix

MUM&Co is able to detect: <br/>
Deletions (unique and mobile), insertions (novel and mobile) and tandem duplications (>50bp) <br/>
Inversions and translocations (>1kb)

Output is a tabulated file containing columns is this order: <br/>
Reference chromosome; Query chromosome; reference start; reference end; size; SV type; query start; query end.
