		    OLIZ 

=============================================
+      OLIZ is free for academic use.       +
+ Please see the  License file for details +
============================================= 
 
Copyright(c)2001-2002 Hao Chen(hchen@utmem.edu)



What is Oliz
------------
Oliz is a suite of Perl scripts that assist the design of oligonucleotide microarrays. Four Perl scripts are provided. The script uni.pl extracts UniGene clusters, which are assembled into contig(s) by contig.pl using the CAP3 assembly program (Huang and Madan, 1999). Then, utr.pl parses the 3'UTRs of the contigs and converts the output into HTML format for visual inspection. Lastly, uniq.pl runs EMBOSS prima to obtain 50mer oligos and performs blast searches to ensure gene specificity of the oligos. The shell script onestep.sh ties the above four scripts together to provide a one step operation.

What is included
----------------
Four Perl scripts
One shell script
All the files produced during a sample run

Programs/Files that are not included and how to obtain them (all publicly available)
--------------------------------------------------------------------------------------
1. Linux/Unix Operating system (I use redhat 7.1 www.redhat.com)
2. Perl (I use v5.6.0. www.Perl.org)
3. bioinformatics tools:
	BioPerl (base install, no extra package necessary. www.bioperl.org)
	blast (ftp://ncbi.nlm.nih.gov/blast/executables/)
	cap3 (please contact Dr. Xiaoqiu Huang huang@mtu.edu)
	clustalw (ftp://ftp-igbmc.u-strasbg.fr/pub/ClustalX/)
	EMBOSS prima (ftp://ftp.uk.embnet.org/pub/EMBOSS)
4. Databases:
	* UniGene (ftp://ncbi.nlm.nih.gov/repository/UniGene/) (required files: Your_species.seq.all, Your_species.data)
	* Blastn database 
		This database determines the specificity of the final oligo sequences. Several options are available. I recommend that you retrieve all the known sequences of your species from NCBI. To do that, go to http://www.ncbi.nlm.nih.gov/ , select Genbank, type in:
		rattus [organism] 
		replace "rattus" with the name of the species you are working with. Click the "Go" button. Change the DISPLAY option to "FASTA" (default is "Summary"). Click DISPLAY (This is important!). Then, on the new page that has no side bar, click the "save" button. You'll be prompted for if you want to save xxx number of sequences, click "yes" and save the file. Format the file following instructions in the blast README files. Alternatively, you might be able to use the UniGene *.seq.all file as your database, or if you already have the blast nt database, you can obtain a gi list to limit the search to your species (change display option to "gi list" and save, see below). However, these files usually contain less sequences than a fresh download from NCBI. 
		
	* A gi list that does not contain EST 
		Sometimes all the oligos produced by the program fail to pass the specificity screen when blasted against all the known sequences of that particular species, including ESTs. Since EST sequences are not always accurate, Oliz runs blast a second time, while limiting the similarity search to none EST sequences. This search requires a no EST gi list. To obtain this list, go to NCBI website: http://www.ncbi.nlm.nih.gov/ select Genbank, then, type in:		
		rattus [organism] NOT gbdiv_est[PROP]
replace "rattus" with the species you are working with. Click on the "Go" button to run the search. Select from the drop down menu "Gi list" (the default is "summary"). Click on the "Save" button, rename the file and save it. Remember the directory because you'll need that information the first time you run oliz.		 
	
	
5. Recommended hardware: 
	Intel Pentium IV, 1G Hz CPU
	1Gb RAM (or close to the size of your blast database)


An example:
----------
Once you have all these ready, here is an example of a basic run:
	1. unpack the package and put them in a directory in your path, such as /home/your_username/bin, make sure all files are executable.
	2. Identify unigenes of your interest, save their ID in a file, one ID per line. For example, we save the following IDs in my_rn.list.

Rn.220
Rn.999
Rn.1345

	3. run:
	% perl onestep.pl my_rn.list
	You'll see a greeting screen and be asked where the required programs/files are located. Type in the path, pay attention to if file name is also required. You only need to do this once. If everything goes well, you'll have your oligo sequences in a few minutes.
	
	4. You can also execute one script at a time. That allows you to modify the intermediate results. Start by running:
	%perl uni.pl my_rn.list
	There will be a few lines to tell you what the next step is at the end of every step.
	
Brief explanation of files
-----------
Oliz produces many files along the way. Here is a brief explanation of the file and directories generated, still assuming the file name is my_rn:

my_rn_blast:	directory stores 50mer blast results
my_rn.brief:	abbreviated version of my_rn.search_log
my_rn_clustal:	directory stores 50mer clustalw results
my_rn.final:	qualified oligos.
my_rn.html:	the grand html file that summarizes the contig assembly and utr parsing. Many links are provided.
my_rn.list: 	a list containing UniGene IDs (this is the input file)
my_rn.missing: 	the file contains UniGene IDs that are not found in the UniGene database used. These IDs might have retired. Duplicate IDs are also provided in this file
my_rn.orient: 	provides fasta formated contigs. The orientation of these contigs cannot be figured out from within UniGene. You are recommended to use blast2 or other programs to compare this contig with its homologue in other species to figure out its orientation. 
my_rn.ref:	orientation information obtained from blast2 is stored in this file. '+' for 5'-3', '-' for 3'-5' (without quotes, of course).
my_rn.search_log:detailed log of the analysis of blast hits
my_rn: 		the main directory that hold files for my_rn.list




