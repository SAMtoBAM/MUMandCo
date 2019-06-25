![alt text](https://github.com/SAMtoBAM/MUMandCo/blob/master/MUM%26Co2.png)

MUM&Co is a simple bash script that uses Whole Genome Alignment information provided by MUMmer (v3 and v4) to detect variants. <br/>

MUM&Co is able to detect: <br/>
Deletions (unique and mobile), insertions (novel and mobile) and tandem duplications (>50bp) <br/>
Inversions and translocations (>1kb)

MUM&Co requires two steps prior to SV detection:

1. Alignment to the reference which will be used to call SVs and renaming based on paired chromosomes <br/>
Contigs corresponding to 2 reference chromosomes should have all reference chromsome names split by an underscore <br/> 
If multiple Contigs correspond to the same chromosome, they can also be split by an underscore then given a unique identifier such as '*_contig1' and *_contig2' <br/>
Tools such as RaGOO and Ragout can do this alongside scaffolding of contigs

2. Perform alignment and filtering with MUMmer3 or 4 and place all files in folder <br/>
e.g. For genome size less than 500Mb <br/>

         nucmer --maxmatch -l 100 -c 500 -p ref reference_assembly.fa query_assembly.fa
         delta-filter -m ref.delta > ref.delta_filter
         show-coords -T -r -c -l -d -g ref.delta_filter.coordsg
         show-coords -T -r -c -l -d ref.delta_filter.coords

         nucmer --maxmatch -l 100 -c 500 -p query query_assembly.fa reference_assembly.fa 
         delta-filter -m query.delta > query.delta_filter
         show-coords -T -r -c -l -d -g query.delta_filter.coordsg
         show-coords -T -r -c -l -d query.delta_filter.coords
    
         mkdir query_alignments
         mv ref.delta* query alignments/ 
         mv query.delta* query_alignments/ 

NOTE: <br/>
Alignment options (-l -c) will depend on genome size <br/>
Alignment and either RaGOO or Ragout steps should be verified by dotplot 

3. Edit mumandco.sh bash script. <br/> 
Provide chosen ref and query names, path to alignments folder (as used above) and desired output prefix in the space provided. <br/>
e.g. <br/>

         ##output file names for alignment and filtering## 
         ref="yeast_tidy" 
         query="DEL100" 
         alignments_folder="DEL100_alignments" 

         ##output file name for SV detection## 
         prefix="DEL100_test" 

4. Run mumandco.sh script. <br/>
Alignments, a summary file and a file with all SVs will be placed in a folder with the chosen prefix.
