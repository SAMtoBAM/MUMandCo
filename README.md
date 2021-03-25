![alt text](https://github.com/SAMtoBAM/MUMandCo/blob/master/MUM%26Co2.png)

MUM&Co is a simple bash script that uses Whole Genome Alignment information provided by MUMmer (v3 and v4) to detect variants. <br/>

MUM&Co is able to detect: <br/>
Deletions (novel and mobile), insertions (novel and mobile) and tandem duplications (>50bp) <br/>
Inversions and translocations (>1kb)

MUM&Co requires installation of MUMmer3 or 4.<br/>
MUM&Co will look for the MUMmer toolkit's scripts path using 'which xxxxx'.<br/>
This path can be editted directly in the script if required.

In order to help with downstream analysis: <br/>
Renaming and re-orientation of the query genome contigs to correspond to their reference counterparts <br/>
Tools such as RaGOO and Ragout can do this alongside scaffolding of contigs (this is not currently recommended for short-read based assemblies).<br/>

Options: <br/>

         -r or --reference_genome          path to reference genome
         -q or --query_genome              path to query genome
         -g or --genome_size               size of genome
         -o or --output                    output prefix

Test run script (keep order of options):
         
         bash mumandco_v*.sh -r ./yeast.tidy.fa -q ./yeast_tidy_DEL100.fa -g 12500000 -o DEL100_test
         
Output folder contains:<br/>
Folder with alignments used for SV detection<br/>
Txt file with summary of SVs detected<br/>
and a tsv file with all the detected SVs.

UPDATE: <br/>
For those using MUMmer4 you can use the MUM&Co script '..MUMmer4_multithreads.sh' that contains an option for increasing the thread use of the nucmer alignments
This can only be run if MUMmer4 is being used as MUMmer3 does not offer this option <br/>
Options: <br/>

         -r or --reference_genome          path to reference genome
         -q or --query_genome              path to query genome
         -g or --genome_size               size of genome
         -o or --output                    output prefix
         -t or --threads                   thread number

Test run script (keep order of options): <br/>
         
         bash mumandco_v*_MUMmer4_multithreads.sh -r ./yeast.tidy.fa -q ./yeast_tidy_DEL100.fa -g 12500000 -o DEL100_test -t 2


Additionally there is an option to search for insertion and deletion events in the reference in order to label them as either mobile or novel events.<br/>
This requires BLAST and samtools installation.<br/>
In order to use this feature, edit the bash script to change 'blast_step="no"' to 'yes'


Reference:<br/>
Samuel O’Donnell, Gilles Fischer, MUM&Co: accurate detection of all SV types through whole-genome alignment, Bioinformatics, Volume 36, Issue 10, 15 May 2020, Pages 3242–3243, https://doi.org/10.1093/bioinformatics/btaa115
