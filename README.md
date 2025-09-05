 <p align="center" >
    <img src="https://github.com/SAMtoBAM/MUMandCo/blob/master/logo.svg" width=70%>
</p>

v3.8 release : [![DOI](https://zenodo.org/badge/187255031.svg)](https://zenodo.org/badge/latestdoi/187255031)

MUM&Co is a simple bash script that uses Whole Genome Alignment information provided by MUMmer (v4) to detect variants. <br/>
VERSION >= 3 UPDATE <br/>
Only uses MUMmer4 now and contains a thread count option <br/>
Contains a VCF output file with all calls currently being imprecise <br/>
Contains another output file containing the calls alongside the respective DNA impacted <br/>
This new step requires samtools installation <br/>
Now calls the reverse of tandem duplications, tandem contractions (>50bp) <br/>

MUM&Co is able to detect: <br/>
Deletions, insertions, tandem duplications and tandem contractions (>=50bp & <=150kb) <br/>
Inversions (>=1kb) and translocations (>=10kb)

MUM&Co requires installation of MUMmer4 and samtools.<br/>
MUM&Co will look for the MUMmer toolkit's and samtools scripts path using 'which xxxxx'.<br/>
An error warning will print and the script will stop if these paths cannot be found <br/>
This path can be editted directly in the script if required.

In order to help with downstream analysis: <br/>
Renaming and re-orientation of the query genome contigs to correspond to their reference counterparts <br/>
Tools such as RaGOO and Ragout can do this alongside scaffolding of contigs (this is not currently recommended for short-read based assemblies) <br/>

Options: <br/>

    Required inputs:
	   -r | --reference_genome		Fasta file containing an assembly
	   -q | --query_genome		    Fasta file containing another assembly
	   -g | --genome_size		    Rough estimation of genome size for both reference and query to determine alignment parameters

	   Recommended inputs:
	   -t | --threads			    Number of threads for alignment (default: 1)
	   -ml | --minlen			    Minimum length of alignments in basepairs (Default: 50)

	   Optional parameters:
	   -p | --prefix			    Prefix for output files and name of output folder ('prefix'_output) (Default: mumandco)
	   -b | --blast			        Adds the blast option to identify is insertions or deletions look repetitive or novel (takes significantly longer)
	   -h | --help			        Print this help message

Test run script: <br/>
         
    mumandco.sh -r reference.fa -q query.fa -g 12500000

OUTPUT FOLDER:<br/>
Folder with alignments used for SV detection <br/>
Txt file with summary of SVs detected <br/>
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
As such, each translocation fragment as called as an event, is now a breakend-like call and will be duplicated if both borders are involved in translocations <br/>

VCF TRA EVENT: <br/>
The later notation for the TSV file is currently being added to the alt column in the VCF for 'TRA' events. <br/>
Currently it is not a called a breakend site (contains no nucleotide at edge) but can be interpreted similarly <br/>

Note: <br/>
MUMmer4 is now required due to the hard wired thread option not available during alignment with MUMmer3 <br/>
The blast option (-b /--blast) using BLAST to search for insertion and deletion events in the reference/query in order to label them as either mobile or novel events.<br/>


Reference:<br/>
Samuel O’Donnell, Gilles Fischer, MUM&Co: accurate detection of all SV types through whole-genome alignment, Bioinformatics, Volume 36, Issue 10, 15 May 2020, Pages 3242–3243, https://doi.org/10.1093/bioinformatics/btaa115
