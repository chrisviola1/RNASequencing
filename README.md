RNA Sequencing

These python scripts use various sequencing software to locate error types and combine mutiple sequencing techniques for increased performance.

**First File**: MatchingMusclePbdagcon.py<br/>
The source script. This uses both the minimap2 error detection and matching script explained below. Takes two files with similar sequences and finds errors within the seqeunces as defined by minimap2. Information is parsed from minimap2 output and manipulated to determine where errors match. Each file is intended to be sequenced by different techniques (ex. Fasta file 1 is muscle and fasta file 2 is pbdagacon).
  
  Script Commands:
  
    FastaFile1: Sequences that need to be aligned in fasta format. Muscle Sequences.
    FastaFile2: Sequences that need to be aligned in fasta format. Pbdagcon Sequences.
    Reference: A reference sequence in fasta file format.
    OutputCsv: The name of the csv output.
    Minimap2: The location of the minimap2 application locally or on a server. Necessary only if minimap2 is not in the current directory.
  Example command:
  ```console
  foo@bar:~$ python MatchingMusclePbdagcon.py fastafile1.fasta fastafile2.fasta reference.fasta output.csv
  ```
 Csv file output:

  Name: ```The name of the strand```<br/>
  Similar Mismatches Intersection: ```Number of errors occurring in both files```<br/>
  Mismatches [fasta file 1 name]: ```Number of errors file 1```<br/>
  Mismatches [fasta file 2 name]: ```Number of errors file 2```<br/>
  NonHP Insertion Similar from Reference: ```number of non-homopolymers insertion errors```<br/> 
  NonHP Deletion Similar from Reference: ```number of non-homopolymers deletion errors```<br/> 
  Substitution Similar from Reference: ```number of substitution errors```<br/> 
  HP Insertion Similar from Reference: ```number of homopolymers insertion errors```<br/> 
  HP Deletion Similar from Reference: ```number of homopolymers deletion errors```<br/> 
  Long Insertion Similar from Reference: ```number of insertion errors longer than a length of 1```<br/>
  Long Deletion Similar from Reference: ```number of deletion errors longer than a length of 1```<br/>
  Reference NonHP Insertion Positions: ```the actual position numbers of non-homopolymers insertion errors```<br/> 
  Reference NonHP Deletion Positions: ```the actual position numbers of non-homopolymers deletion errors```<br/> 
  Reference Substitution Positions: ```the actual position numbers of substitution errors```<br/> 
  Reference HP Insertion Positions: ```the actual position numbers homopolymers insertion errors```<br/> 
  Reference HP Deletion Positions: ```the actual position numbers homopolymers deletion errors```<br/>
  Reference Long Insertion Positions: ```the actual position numbers of insertion errors longer than a length of 1```<br/>
  Reference Long Deletion Positions: ```the actual position numbers of deletion errors longer than a length of 1```<br/>

The minimap2 output explained below for each fasta file is concatenated to these results. 
<br/>
<br/>
<br/>
<br/>

**Second File**: Minimap2 Error Detection<br/>
Parses through minimap2 data and tabulates error type information into a clean and neat csv file.
  
  Script Commands:
  
    FastaFile: A fasta file format with numerous sequences.
    Reference: A reference sequence in fasta file format.
    OutputCsv: The name of the csv output.
    Minimap2: The location of the minimap2 application locally or on a server. Necessary only if minimap2 is not in the current directory.
  Example command:
  ```console
  foo@bar:~$ python Minimap2_Error_Detection.py Pbdagcon_Optimizers/Consensus_result_pbdagcon_pNLnefSBR2_ddAmp368Kb_RCA.fasta_penalty0no1.0nan.fasta NL4.3_reference.txt Pbdagcon_Optimizers/fasta_penalty0.csv
  ```
 Csv file output:
  ```console
  Name:
  Length:
  Reference Start Position:
  Total Matching Range:	
  Misaligned Front:
  Aligned:
  Misaligned Back:
  Subreads:
  Mismatches:
  Reference NonHP Insertion:
  Reference NonHP Deletion:
  Reference Substitution:
  Reference HP Insertion:
  Reference HP Deletion:
  Reference Long Insertion:
  Reference Long Deletion:
  Consensus NonHP Insertion Positions:
  Consensus NonHP Deletion Positions:
  Consensus Substitution Positions:
  Consensus HP Insertion Positions:
  Consensus HP Deletion Positions:
  Reference NonHP Insertion Positions:
  Reference NonHP Deletion Positions:
  Reference Substitution Positions:
  Reference HP Insertion Positions:
  Reference HP Deletion Positions:
  Reference Long Insertion Sites:
  Reference Long Deletion Sites:
  ```
  
  
