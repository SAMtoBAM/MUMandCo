![alt text](https://github.com/SAMtoBAM/MUMandCo/blob/master/MUM%26Co2.png)

MUM&Co is a simple bash script that uses Whole Genome Alignment information provided by MUMmer (v3 and v4) to detect variants. <br/>

MUM&Co is able to detect: <br/>
Deletions (unique and mobile), insertions (novel and mobile) and tandem duplications (>50bp) <br/>
Inversions and translocations (>1kb)

MUM&Co requires installation of MUMmer3 or 4.<br/>
MUM&Co will look for the MUMmer toolkit's scripts path using 'which xxxxx'.<br/>
This path can be editted directly in the script if required.

MUM&Co requires 1 step prior to SV detection: <br/>
Alignment to the reference which will be used to call SVs and renaming based on paired chromosomes <br/>
Contigs corresponding to 2 reference chromosomes should have all reference chromsome names split by an underscore<br/>
e.g. 'chromosome2_3' <br/> 
If multiple Contigs correspond to the same chromosome, they should be split by an underscore and a unique identifier<br/>
e.g. 'chromosome2_contig1' and chromosome2_contig2' <br/>
Tools such as RaGOO and Ragout can do this alongside scaffolding of contigs (this is not currently recommended for short-read based assemblies).<br/>
Please verify renaming and scaffolding.

Options: <br/>

         -r or --reference_genome          path to reference genome
         -q or --query_genome              path to query genome
         -g or --genome_size               size of genome
         -o or --output                    output prefix

Test run script (keep order of options):
         
         bash mumandco.sh -r ./yeast.tidy.fa -q ./yeast_tidy_DEL100.fa -g 12500000 -o DEL100_test
         
Output folder contains:<br/>
Folder with alignments used for SV detection<br/>
Txt file with summary of SVs detected<br/>
and a tsv file with all the detected SVs.
