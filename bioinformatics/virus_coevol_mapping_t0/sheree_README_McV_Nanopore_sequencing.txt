##################

#Data delivery for by nanopore of Micromonas commoda virus

##################

#02/07/2024
@ DATA DELIVERY OF McV 
Resequencing was done by Nanopore on McV by Aurelie Claes from DNA extracted by Thomas Anicet.


Bonjour à tous,

Le séquençage MinION est terminé et vous pourrez trouvez les résultats en cliquant sur le lien ci-dessous :

https://nuage.obs-banyuls.fr/s/59yHLotGs3Takk5

Bonne journée

Aurélie
barcode 	sample ID
06 	MPV
07 	N3
08 	N6
09 	N9
10 	P19
11 	S50-2
12 	S50-7
13 	S50-9
14 	T14-12
15 	T14-23
16 	T14-25

--> downloaded from Nuage the folder barcode06-MPV/
Containing 358 fastq.gz files.

#18/12/2024
I unzipped all of them 
% gzip -d FAZ18065_pass_barcode06_68193fbb_85fe0de7_*.gz   

I concatenated into one big fastq
% cat FAZ18065_pass_barcode06_68193fbb_85fe0de7_* >> all_McV.fastq

How many reads in the all_McV.fastq?
 % grep -c "@" all_McV.fastq                                  
---> 158477 fastq reads

Check that the sum of all the individual files is the same as the concatenated
% grep -c "@" FAZ18065_pass_barcode06_68193fbb_85fe0de7_* > checksum.txt
--> yes checksum correct.

Copied to my account on Genotoul
% scp all_McV.fastq syau@genobioinfo.toulouse.inrae.fr:/home/syau/work/project_Mcommoda_virus/data_sequencing_McV_nanopore/



#20/12/2024
@ASSEMBLE THE McV READS

Wrote a script based on the test file in genobioinfo and what is on the flye github
https://github.com/mikolmogorov/Flye/blob/flye/docs/USAGE.md



#!/bin/bash
#SBATCH -p workq
#SBATCH -t 01-00:00:00 
#SBATCH --cpus-per-task 4
#SBATCH --mem=40G

#Load modules
module load bioinfo/Flye/2.9.2

flye --nano-raw /home/syau/work/project_Mcommoda_virus/data_sequencing_McV_nanopore/all_McV.fastq --out-dir out_McV_assem1  --threads 4


--> assembly failed. Probably because there is too many other reads. Try on meta genome mode
-----------End assembly log------------
[2024-12-20 14:04:08] root: ERROR: No disjointigs were assembled - please check if the read type and genome size parameters are correct
[2024-12-20 14:04:08] root: ERROR: Pipeline aborted

#test2 in meta genome mode that takes into account uneven coverages.
#!/bin/bash
#SBATCH -p workq
#SBATCH -t 01-00:00:00 
#SBATCH --cpus-per-task 4
#SBATCH --mem=40G

#Load modules
module load bioinfo/Flye/2.9.2

flye --nano-raw /home/syau/work/project_Mcommoda_virus/data_sequencing_McV_nanopore/all_McV.fastq --out-dir out_McV_assem1 --genome-size 250k --meta  --threads 4

--> This time it worked.
Assembly info
#seq_name	length	cov.	circ.	repeat	mult.	alt_group	graph_path
contig_83	214624	768	N	N	4	*	83,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,6,-5,-5,-5,-5,-5,-5,-5,-5,-5
contig_22	29628	3	N	N	1	*	*,22,*
contig_13	27715	4	N	N	1	*	*,13,*
contig_61	26602	3	N	N	1	*	*,61,*
contig_29	24919	3	N	N	1	*	*,29,*
contig_23	21671	4	N	N	1	*	*,23,*
contig_44	20159	3	N	N	1	*	*,44,*
contig_45	20002	3	N	N	1	*	*,45,*
contig_72	19691	3	N	N	1	*	*,72,*
contig_11	16400	3	N	N	1	*	*,11,*
contig_57	15638	4	N	N	1	*	*,57,*
contig_80	15450	4	N	N	1	*	*,80,*
contig_26	12988	5	N	N	1	*	*,26,*
contig_73	12689	4	N	N	1	*	*,73,*
contig_63	12526	3	N	N	1	*	*,63,*
contig_10	12131	3	N	N	1	*	*,10,*
contig_47	10364	5	N	N	1	*	*,47,*
contig_1	10354	3	N	N	1	*	*,1,*
contig_33	10106	3	N	N	1	*	*,33,*
contig_4	9735	4	N	N	1	*	*,4,*
contig_71	9602	4	N	N	1	*	*,71,*
contig_28	9501	4	N	N	1	*	*,28,*
contig_77	9416	4	N	N	1	*	*,77,*
contig_82	9116	3	N	N	1	*	*,82,*
contig_34	9039	3	N	N	1	*	*,34,*
contig_25	8414	3	N	N	1	*	*,25,*
contig_38	8310	3	N	N	1	*	*,38,*
contig_35	7782	4	N	N	1	*	*,35,*
contig_18	7389	3	N	N	1	*	*,18,*
contig_53	7149	3	N	N	1	*	*,53,*
contig_21	7031	3	N	N	1	*	*,21,*
contig_49	7011	3	N	N	1	*	*,49,*
contig_20	6893	3	N	N	1	*	*,20,*
contig_3	6609	4	N	N	1	*	*,3,*
contig_15	6539	3	N	N	1	*	*,15,*
contig_32	6394	3	N	N	1	*	*,32,*
contig_19	6293	3	N	N	1	*	*,19,*
contig_8	6174	4	N	N	1	*	*,8,*
contig_14	5912	4	N	N	1	*	*,14,*
contig_70	5884	3	N	N	1	*	*,70,*
contig_78	5109	3	N	N	1	*	*,78,*
contig_54	5077	3	N	N	1	*	*,54,*
contig_81	4790	4	N	N	1	*	*,81,*
contig_75	4786	3	N	N	1	*	*,75,*
contig_42	4305	5	N	N	1	*	*,42,*
contig_9	4279	4	N	N	1	*	*,9,*
contig_39	4059	4	N	N	1	*	*,39,*
contig_27	3937	7	N	N	1	*	*,27,*
contig_55	3726	4	N	N	1	*	*,55,*
contig_43	3322	3	N	N	1	*	*,43,*
contig_74	3319	3	N	N	1	*	*,74,*
contig_7	3219	4	N	N	1	*	*,7,*
contig_66	3166	5	N	N	1	*	*,66,*
contig_12	3123	4	N	N	1	*	*,12,*
contig_52	3042	3	N	N	1	*	*,52,*
contig_31	3040	4	N	N	1	*	*,31,*
contig_41	3010	3	N	N	1	*	*,41,*
contig_50	2944	3	N	N	1	*	*,50,*
contig_17	2911	3	N	N	1	*	*,17,*
contig_67	2767	3	N	N	1	*	*,67,*
contig_60	2701	4	N	N	1	*	*,60,*
contig_37	2693	3	N	N	1	*	*,37,*
contig_79	2682	4	N	N	1	*	*,79,*
contig_16	2662	3	N	N	1	*	*,16,*
contig_62	2536	4	N	N	1	*	*,62,*
contig_36	2531	3	N	N	1	*	*,36,*
contig_76	2514	3	N	N	1	*	*,76,*
contig_58	2478	4	N	N	1	*	*,58,*
contig_2	2466	5	N	N	1	*	*,2,*
contig_51	2417	4	N	N	1	*	*,51,*
contig_30	2386	3	N	N	1	*	*,30,*
contig_59	2352	4	N	N	1	*	*,59,*
contig_46	2279	3	N	N	1	*	*,46,*
contig_24	2276	5	N	N	1	*	*,24,*
contig_69	2275	3	N	N	1	*	*,69,*
contig_56	2269	4	N	N	1	*	*,56,*
contig_68	2255	3	N	N	1	*	*,68,*
contig_40	2225	4	N	N	1	*	*,40,*
contig_48	2207	4	N	N	1	*	*,48,*
contig_65	2150	4	Y	Y	1	*	65
contig_64	1492	7	N	N	1	*	*,64,*

downloaded to my MacBook
 % scp syau@genobioinfo.toulouse.inrae.fr:/home/syau/work/project_Mcommoda_virus/analysis_flye_assembly_McV/out_McV_assem1/assembly.fasta .

I assume that the longest contig with the highest coverage contig_83 is the McV genome, as it has coverage >700x while the others have only 10s of coverage.

Extract the contig_83 with seqtk on your MacBook
% echo "contig_83" > contig_list.txt
% ./seqtk subseq /Users/shereeyau/Documents/work/Lucien_Weymiens_masters_project/analysis_McV_flye_assembly/assembly.fasta /Users/shereeyau/Documents/work/Lucien_Weymiens_masters_project/analysis_McV_flye_assembly/contig_list.txt > /Users/shereeyau/Documents/work/Lucien_Weymiens_masters_project/analysis_McV_flye_assembly/contig_83.fasta

Submit contig_83 to blast on NCBI
--> megablast got no hits.
Try with blast
--> still no hits

I took the first 5858 bp  of the contig and did a megablast 
--> 
Description	max_score	total_score	qcov	perc_id	acc_length	accession
Micromonas commoda virus 20T, complete genome	Micromonas commoda virus 20T	2327	4837	51%	0.0	96.77%	209762	PQ442261.1
Micromonas commoda virus strain McV-KB4, complete genome	Micromonas commoda virus	713	862	24%	0.0	76.58%	212418	PQ359806.1

I did a blast2seq of contig_83 against McV 20T
--> Matches well but patchily over the length in the same direction
Micromonas commoda virus 20T, complete genome	Micromonas commoda virus 20T	14657	1.586e+05	66%	0.0	86.52%	209762	PQ442261.1



# Also try assembling after filtering out the reads >10000 bp
On your MacBook use seqtk to remove reads less than 10000
% ./seqtk seq -L 10000 /Users/shereeyau/Documents/work/Lucien_Weymiens_masters_project/data_DNA_extraction_sequencing_McV_nanopore/barcode06-MPV/all_McV.fastq > /Users/shereeyau/Documents/work/Lucien_Weymiens_masters_project/data_DNA_extraction_sequencing_McV_nanopore/barcode06-MPV/all_McV_gt10000.fastq

Upload to genobioinfo.
% scp all_McV_gt10000.fastq syau@genobioinfo.toulouse.inrae.fr:/home/syau/work/project_Mcommoda_virus/data_sequencing_McV_nanopore/ 

work in progress...

#15/01/2025
@FIND THE MCP GENE OF McV
Go on NCBI And look up the genome of McV 20T PQ442261.1
Do a keyword search for "capsid"
--> There are 8 loci

Look up the 6th copy of the MCP, as this is the homolog in OtV5
-->
>XKM46595.1 major capsid protein [Micromonas commoda virus 20T]
MAGGLMQLVAYGAQDVYLTGNPEVTFFQAKYKRHTNFAMENIEQTVNGTAADSGRVSVTIARNGDLVGDM
YVELLSAAAASISSDATDDSCWVAERAISSVEISIGGQRVDKHYQKWWRLYSELYLDESKKLTYGKMTSA
TTGGAVYLPLVFFFNRNPGLYLPLIALQYHEVRIDFDLASDFSTYLNTGTFKVWANYVYLDTEERRRFAQ
KGHEYLIEQVQHTGQDTVTASGGTKQVRLSYNHPVKELVWCCDEGVARTKMWNFTHKAQVAEIVLEQDLT
MADSNCFIAPGAAGAPLLVCGTGGGTSKFTEEAVGTIDKFKLVLNGQDRFKEQSGKYFNQVQPHFHHSGA
PYAGVYAYSFALKPEEHQPTGTCNFSRIDNAQVSITTTSGNDAATNLNMFAVNYNVLRVQSGMGGLAFSN

Here is the gene sequence of the 6th MCP copy
atggccggtg gtcttatgca actcgtagct
   156361 tatggtgccc aggatgtcta ccttaccggt aaccctgagg taactttctt ccaggcgaaa
   156421 tacaagcgcc acactaactt cgcgatggag aacatcgagc agaccgtcaa cggtactgcc
   156481 gctgactccg gtcgcgtctc cgtcaccatt gcccgtaacg gtgatctcgt cggcgacatg
   156541 tacgtcgagc tcctctccgc cgctgcggcg tccatctctt ccgatgccac tgacgattct
   156601 tgctgggtcg ctgagcgtgc gatctcctcc gtcgagatat cgatcggtgg acaaagggtg
   156661 gacaagcact accagaagtg gtggcgtctc tactccgagc tttaccttga cgagtccaag
   156721 aagctcactt acggtaagat gacttccgcc actactggcg gtgctgtcta tttgccccta
   156781 gtctttttct ttaaccgcaa tcccggtctc tatctcccac taattgctct gcagtaccat
   156841 gaggtccgta tcgatttcga tttagcgtct gatttcagca cctatctcaa taccggtacc
   156901 ttcaaggtgt gggccaacta cgtctacctt gacactgagg agcgtaggcg ttttgcccag
   156961 aagggccacg agtacctgat cgagcaggtt caacacaccg gtcaggacac tgttaccgct
   157021 tccggtggta ccaagcaggt ccgcctctcg tacaatcacc ccgtcaagga gctcgtgtgg
   157081 tgctgcgacg agggtgtcgc ccgtaccaag atgtggaact ttacccacaa ggcccaggtt
   157141 gctgagattg ttctcgagca ggacctcacc atggcggact ccaactgttt catcgccccc
   157201 ggcgccgcgg gtgcccctct ccttgtgtgc ggcaccggtg gtggcacttc caagttcacc
   157261 gaggaggctg tcggtaccat cgacaagttc aagcttgtgc ttaacggcca ggaccgcttc
   157321 aaggagcagt ctggtaagta cttcaaccag gtgcagcccc acttccacca ctccggcgcc
   157381 ccctacgcgg gtgtctacgc gtactccttc gcgctcaagc ccgaggagca ccagcctacc
   157441 ggcacttgca acttctcccg tatcgataac gcgcaggttt ccatcaccac cacctccggt
   157501 aacgatgccg cgaccaacct caacatgttc gctgttaact acaacgtcct ccgtgtccag
   157561 tcgggtatgg gtggccttgc cttctccaac taa

On the blast webpage, do a blastp of the McV20T - MCP above to check it is the same as the OtV5 homologue
-->
The best blastp match of the McV20T gene above is to OtV5 locus OtV5_170, which is indeed the 6th copy and the likely "real" homologue.

Select seq gb|AET43572.1|	hypothetical protein MPWG_00083 [Micromonas pusilla virus PL1]	Micromonas pusilla virus PL1	677	677	100%	0.0	76.25%	418	AET43572.1
Select seq ref|YP_001648266.1|	putative major capsid protein [Ostreococcus tauri virus OtV5]	Ostreococcus tauri virus OtV5	661	661	100%	0.0	76.76%	432	YP_001648266.1
Select seq gb|XKM47335.1|

Now try and do a tblastn of the McV_20T MCP protein sequence against contig_83.fasta nucleotide on the NCBI webserver
--> This gave 8 matches!!! Probably the 8 homologues of the MCP gene
--> top match at 74% identity
--> The second match is at 31%, so the first match is right, locus is 65433 ...64180
Score	Expect	Method	Identities	Positives	Gaps
Frame
657 bits(1696)	0.0	Compositional matrix adjust.	313/422(74%)	351/422(83%)	4/422(0%)
-2
Query  1      MAGGLMQLVAYGAQDVYLTGNPEVTFFQAKYKRHTNFAMENIEQTVNGTAADSGRVSVTI  60
              MAGGLMQLVAYGAQDVYLTGNPEVTFFQAKYKRHTNFAMENIEQTVNGTAA+SGRVSVT+
Sbjct  65433  MAGGLMQLVAYGAQDVYLTGNPEVTFFQAKYKRHTNFAMENIEQTVNGTAANSGRVSVTV  65254

Query  61     ARNGDLVGDMYVELLsaaaasissdatddsC-WVAERAISSVEISIGGQRVDKHYQKWWR  119
              ARNGDLVGDMY+EL   +  + +       C WVAERA+++VE+SIGGQR+DKHYQKWWR
Sbjct  65253  ARNGDLVGDMYIEL--ESDEATTITTAAADCNWVAERAVNNVELSIGGQRIDKHYQKWWR  65080

Query  120    LYSELYLDESKKLTYGKMTSATTGGAVYLPLVFFFNRNPGLYLPLIALQYHEVRIDFDLA  179
              +YSELYLDESKK T+GKMT+A  G  VYLPL+FFFNRNPGL LPLIALQYHEVRIDFDLA
Sbjct  65079  MYSELYLDESKKATWGKMTTAGDGKTVYLPLIFFFNRNPGLALPLIALQYHEVRIDFDLA  64900

Query  180    SDFSTYLNTGTFKVWANYVYLDTEERRRFAQKGHEYLIEQVQHTGQDTVTASGGTKQVRL  239
              S+F+TYLN   FKVWANYVYLDTEERRRFAQKGHEYLIEQVQHTG DTVTA GGTKQVRL
Sbjct  64899  SNFTTYLNASVFKVWANYVYLDTEERRRFAQKGHEYLIEQVQHTGTDTVTADGGTKQVRL  64720

Query  240    SYNHPVKELVWCCDEGVARTKMWNFTHKAQVAEIVLEQDL-TMADSNCFIAPGAAGAPLL  298
              SYNHPVKELVWC         MWNFT  +  A I L+ +  ++  SNCF+    AG P++
Sbjct  64719  SYNHPVKELVWCFSNTQTNNGMWNFTTASTDANIKLDSNQNSLEGSNCFVTTATAGTPMV  64540

Query  299    VCGTGGGTSKFTEEAVGTIDKFKLVLNGQDRFKEQSGKYFNQVQPHFHHSGAPYAGVYAY  358
                G  GG+S FTEEAVG +  FKL+LNGQDRFKEQ GKYFNQVQP+ HH+G PY G+Y+Y
Sbjct  64539  KVGAIGGSSIFTEEAVGPLSTFKLILNGQDRFKEQKGKYFNQVQPYNHHTGCPYPGIYSY  64360

Query  359    SFALKPEEHQPTGTCNFSRIDNAQVSITTTSGNDAATNLNMFAVNYNVLRVQSGMGGLAF  418
              SFALKPEEHQPTGTCNFSRIDNAQV + T    + A +++MFA NYNVLR+QSGMGGLAF
Sbjct  64359  SFALKPEEHQPTGTCNFSRIDNAQVQVVTAGTTNNAISMHMFATNYNVLRIQSGMGGLAF  64180



Query  419    SN  420
              SN
Sbjct  64179  SN  64174

Range 3: 138488 to 139531GraphicsNext MatchPrevious MatchFirst Match
Alignment statistics for match #3
Score	Expect	Method	Identities	Positives	Gaps
Frame
180 bits(456)	1e-51	Compositional matrix adjust.	139/432(32%)	206/432(47%)	96/432(22%)



I did a blastn of the above McV 20T MCP gene DNA sequence against contig_83.fasta
--> Top match to one locus corresponding to the right locus 65433...64171
contig_83
Sequence ID: Query_6795639Length: 214624Number of Matches: 1
Range 1: 64171 to 65433GraphicsNext MatchPrevious Match
Alignment statistics for match #1
Score	Expect	Identities	Gaps	Strand
898 bits(486)	0.0	1017/1276(80%)	26/1276(2%)	Plus/Minus
Query  1      ATGGCCGGTGGTCTTATGCAACTCGTAGCTTATGGTGCCCAGGATGTCTACCTTACCGGT  60
              ||||| || ||||||||||||||||||||||| |||||||||||||||||||||||||||
Sbjct  65433  ATGGCTGGCGGTCTTATGCAACTCGTAGCTTACGGTGCCCAGGATGTCTACCTTACCGGT  65374

Query  61     AACCCTGAGGTAACTTTCTTCCAGGCGAAATACAAGCGCCACACTAACTTCGCGATGGAG  120
              ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  65373  AACCCTGAGGTAACTTTCTTCCAGGCGAAATACAAGCGCCACACTAACTTCGCGATGGAG  65314

Query  121    AACATCGAGCAGACCGTCAACGGTACTGCCGCTGACTCCGGTCGCGTCTCCGTCACCATT  180
              ||||||||||||||||||||||||||||||||  ||||||| ||||||||||||||| ||
Sbjct  65313  AACATCGAGCAGACCGTCAACGGTACTGCCGCCAACTCCGGCCGCGTCTCCGTCACCGTT  65254

Query  181    GCCCGTAACGGTGATCTCGTCGGCGACATGTACGTCGAGCTCCTCTCCGCCGCTGCGGCG  240
              ||||||||||||||||||||||| ||||||||| |||||||    ||||  |  ||| | 
Sbjct  65253  GCCCGTAACGGTGATCTCGTCGGTGACATGTACATCGAGCTTGAGTCCGATGAGGCGAC-  65195

Query  241    TCCATCTCTTCCGATGCCACTGACGATTCTTGCTGGGTCGCTGAGCGTGCGATCTCCTCC  300
              | | | ||  ||  |||| | |  ||||    |||||| || ||||||||| |   |  |
Sbjct  65194  TAC-TATCA-CCACTGCCGCGGCTGATTGCAACTGGGTTGCCGAGCGTGCGGTTAACAAC  65137

Query  301    GTCGAGATATCGATCGGTGGACAAAGGGTGGACAAGCACTACCAGAAGTGGTGGCGTCTC  360
              || ||  |||| || ||||| ||  |  | ||||||||||||||||||||||||||  | 
Sbjct  65136  GTAGAATTATCAATTGGTGGCCAGCGTATTGACAAGCACTACCAGAAGTGGTGGCGCATG  65077

Query  361    TACTCCGAGCTTTACCTTGACGAGTCCAAGAAGCTCACTTACGGTAAGATGACTTCCGCC  420
              ||||||||||| ||||| || ||||||||||||  |||||  |||||||||||  | || 
Sbjct  65076  TACTCCGAGCTCTACCTCGATGAGTCCAAGAAGGCCACTTGGGGTAAGATGACCACTGCG  65017

Query  421    ACT-ACTGGCGGTGCTGTCTATTTGCCCCTAGTCTTTTTCTTTAACCGCAATCCCGGTCT  479
                | || |||    |||||||  | |||||  | || ||||||||| | ||||| |||||
Sbjct  65016  GGTGAC-GGCAAGACTGTCTACCTCCCCCTTATTTTCTTCTTTAACAGGAATCCTGGTCT  64958

Query  480    CTATCTCCCACTAATTGCTCTGCAGTACCATGAGGTCCGTATCGATTTCGATTTAGCGTC  539
              |   |||||||||||||| ||||||||||||||||| || ||||||||||||||||||||
Sbjct  64957  CGCGCTCCCACTAATTGCCCTGCAGTACCATGAGGTGCGCATCGATTTCGATTTAGCGTC  64898

Query  540    TGATTTCAGCACCTATCTCAATACCGGTACCTTCAAGGTGTGGGCCAACTACGTCTACCT  599
              | | |||| |||||| |||||  |   |  ||||||||| |||||||||||||| |||||
Sbjct  64897  TAACTTCACCACCTACCTCAACGCGTCTGTCTTCAAGGTCTGGGCCAACTACGTGTACCT  64838

Query  600    TGACACTGAGGAGCGTAGGCGTTTTGCCCAGAAGGGCCACGAGTACCTGATCGAGCAGGT  659
              |||||| ||||||||| | || |||||||||||||| ||||||||||| || ||||||||
Sbjct  64837  TGACACCGAGGAGCGTCGCCGATTTGCCCAGAAGGGTCACGAGTACCTCATTGAGCAGGT  64778

Query  660    TCAACACACCGGTCAG-GACACTGTTACCGCTTCCGGTGGTACCAAGCAGGTCCGCCTCT  718
              ||| |||||||| ||  |||||||||||||||   |||||||||||||||||||||||||
Sbjct  64777  TCAGCACACCGG-CACTGACACTGTTACCGCTGATGGTGGTACCAAGCAGGTCCGCCTCT  64719

Query  719    CGTACAATCACCCCGTCAAGGAGCTCGTGTGGTGCTGCGACGAGGGTGTCGCCCGTACCA  778
              | ||||| |||||||| ||||||||||| ||||||| |  | |       | ||  || |
Sbjct  64718  CCTACAACCACCCCGTTAAGGAGCTCGTATGGTGCTTCTCCAACACCCA-GACCA-ACAA  64661

Query  779    AG--ATGTGGAACTTTACCCACAAGGCCCAG-GTTGCTGAGATT--GTTC---TCGAGCA  830
               |  ||||||||||| ||| ||   | | |  | |||  | ||   | ||   || | | 
Sbjct  64660  CGGTATGTGGAACTTCACC-ACCGCGTCTACCGATGCCAACATCAAGCTCGACTCCAACC  64602

Query  831    GGACCTCACCATGGCGGACTCCAACTGTTTCATCGCCCCCGGCGCCGCGGGTGCCCCTCT  890
               || ||| || | | || ||| ||||||||| |  || ||| | |||| ||| ||||| |
Sbjct  64601  AGAACTC-CC-TCGAGGGCTCTAACTGTTTCGTTACCACCGCCACCGCTGGTACCCCTAT  64544

Query  891    CCTTGTGTGCGGCACCGGTGGTGGCACTTCCAAGTTCACCGAGGAGGCTGTCGGTACCAT  950
                ||  |   ||  ||  ||| ||  | ||||  ||||| |||||||| |||||| || |
Sbjct  64543  GGTTAAGGTTGGTGCCATTGGCGGTTCGTCCATCTTCACTGAGGAGGCCGTCGGTCCCCT  64484

Query  951    CGACAAGTTCAAGCTTGTGCTTAACGGCCAGGACCGCTTCAAGGAGCAGTCTGGTAAGTA  1010
              |  ||  ||||||||  | || |||||||||||||| ||||||||||||   || |||||
Sbjct  64483  CTCCACCTTCAAGCTCATCCTCAACGGCCAGGACCGTTTCAAGGAGCAGAAGGGCAAGTA  64424

Query  1011   CTTCAACCAGGTGCAGCCCCACTTCCACCACTCCGGCGCCCCCTACGCGGGTGTCTACGC  1070
              |||||||||||| |||||| ||  ||||||| |||||  ||||||| | ||| ||||| |
Sbjct  64423  CTTCAACCAGGTCCAGCCCTACAACCACCACACCGGCTGCCCCTACCCCGGTATCTACTC  64364

Query  1071   GTACTCCTTCGCGCTCAAGCCCGAGGAGCACCAGCCTACCGGCACTTGCAACTTCTCCCG  1130
              |||||| |||||||||||||||||||||||||||||||||||||| ||||||||||| ||
Sbjct  64363  GTACTCTTTCGCGCTCAAGCCCGAGGAGCACCAGCCTACCGGCACCTGCAACTTCTCGCG  64304

Query  1131   TATCGATAACGCGCAGGTTTCCATCACCACCACCTCCGGTA--ACGATGC-CGCGACCAA  1187
               ||||| |||||||||||  |||      | ||| | ||||  || |    ||||| |  
Sbjct  64303  CATCGACAACGCGCAGGT--CCAGGTTGTC-ACCGCGGGTACCACCAACAACGCGATCTC  64247

Query  1188   CCTCAACATGTTCGCTGTTAACTACAACGTCCTCCGTGTCCAGTCGGGTATGGGTGGCCT  1247
              | |  ||||||||||   ||||||||||||||||||  ||||||||||||||||||||||
Sbjct  64246  CATGCACATGTTCGCCACTAACTACAACGTCCTCCGCATCCAGTCGGGTATGGGTGGCCT  64187

Query  1248   TGCCTTCTCCAACTAA  1263
              ||||||||||||||||
Sbjct  64186  TGCCTTCTCCAACTAA  64171


#16/01/2025
@Extract the MCP gene sequence
Went to the website https://www.reverse-complement.com/ and made the reverse complement of contig_83
--> Saved as contig_83revcomp.fasta
Did word search for the starting and ending ~20 bp of the putative MCP and copied to a new file 
--> contig_83_MCP_cp6.fasta 
