# MUMandCo

MUM&Co is a simple bash script that uses the Whole Genome Alignment information provided by the MUMmer package to detect variants.
MUM&Co requires two steps prior to SV detection.

1. Alignment to the reference which will be used to call SVs and check renaming based on paired chromosomes       
NOTE: Tools such as RaGOO and Ragout can do this alongside scaffolding of contigs

2. Within MUMmer3 or 4, reciprocal nucmer genome alignments, Global and Many-to-Many filtering, and coordinate parsing steps provide the raw data analysed by MUM&Co.

MUM&Co is able to detect:
>50bp: Deletions (unique and mobile), insertions (novel and mobile) and tandem duplications
>1kb: Inversions and translocations
