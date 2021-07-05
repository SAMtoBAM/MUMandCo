![alt text](https://github.com/SAMtoBAM/MUMandCo/blob/master/MUM%26Co2.png)

MUM&Co is a simple bash script that uses Whole Genome Alignment information provided by MUMmer (v4) to detect variants. <br/>
VERSION >= 3 UPDATE <br/>
Only uses MUMmer4 now and requires the thread count option <br/>
Contains a VCF output file with all calls currently being imprecise <br/>
Contains another output file containing the calls alongside the respective DNA impacted <br/>
This new step requires samtools installation <br/>
Now calls the reverse of tandem duplications, tandem contractions (>50bp) <br/>

MUM&Co is able to detect: <br/>
Deletions, insertions, tandem duplications and tandem contractions (>=50bp & <=150kb) <br/>
Inversions (>=1kb) and translocations (>=10kb)

MUM&Co requires installation of MUMmer4 and samtools.<br/>
MUM&Co will look for the MUMmer toolkit's and samtools scripts path using 'which xxxxx'.<br/>
This path can be editted directly in the script if required.

In order to help with downstream analysis: <br/>
Renaming and re-orientation of the query genome contigs to correspond to their reference counterparts <br/>
Tools such as RaGOO and Ragout can do this alongside scaffolding of contigs (this is not currently recommended for short-read based assemblies).<br/>

Options: <br/>

         -r or --reference_genome          path to reference genome
         -q or --query_genome              path to query genome
         -g or --genome_size               size of genome
         -o or --output                    output prefix
         -t or --threads                   thread number

Test run script (keep order of options): <br/>
         
         bash mumandco_v*.sh -r ./yeast.tidy.fa -q ./yeast_tidy_DEL100.fa -g 12500000 -o DEL100_test -t 2

OUTPUT FOLDER:<br/>
Folder with alignments used for SV detection<br/>
Txt file with summary of SVs detected<br/>
TSV file with all the detected SVs <br/>
TSV file with all detected SVs plus the DNA associated with the event (all from reference except insertions) <br/>
VCF file with all calls currently being imprecise <br/>

TSV NOTES: <br/>
The last column in the TSV file contains notes. <br/>
'complicated' : multiple calls within the same region; generally overlapping insertions and deletions <br/>
'double' : several calls at the same coordinates; generally tandem duplications or contractions with multiple copy changes <br/>
']chrX:xxxxxx]' : a VCF inspired notation for the association of the translocation fragments with the other fragments <br/>
e.g. for chr1 with its right border at 250000bp assocaited with chr2 at 100000bp; <br/>
the note would be as follows for chr 1: ']chr2:100000]'     and for chr2 : '[chr1:250000[' <br/>
As such, each translocation fragment as called as an event, is now a breakend-like call and will be duplicated if both borders are involved in translocations <br\>

VCF TRA EVENT: <br/>
The later notation for the TSV file is currently being added to the alt column in the VCF for 'TRA' events. <br/>
Currently it is not a called a breakend site (contains no nucleotide at edge) but can be interpreted similarly <br/>

Note: <br/>
MUMmer4 is now required due to the hard wired thread option not available during alignment with MUMmer3 <br/>
Additionally there is an option to search for insertion and deletion events in the reference in order to label them as either mobile or novel events.<br/>
This requires BLAST installation.<br/>
In order to use this feature, edit the bash script to change 'blast_step="no"' to 'yes'


Reference:<br/>
Samuel O’Donnell, Gilles Fischer, MUM&Co: accurate detection of all SV types through whole-genome alignment, Bioinformatics, Volume 36, Issue 10, 15 May 2020, Pages 3242–3243, https://doi.org/10.1093/bioinformatics/btaa115
