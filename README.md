RNA Sequencing

These python scripts use various sequencing software to locate error types (subsitution, insertion, deletion) and combine mutiple sequencing techniques. This increases detection and validates modification of similar consensus strands. 
<br/>
<br/>
Scripts:<br/>
Script 1: MatchingMusclePbdagcon.py; combines the two scripts below to output the matching data <br/>
Script 2: Minimap2ErrorDetection.py; parses the minimap2 results and outputs them into a csv file <br/>
Script 3: MatchingIdentification.py; finds matching errors from two different Minimap2ErrorDetection.py output csv files <br/>

**First File**: MatchingMusclePbdagcon.zip<br/>
The source script. This uses both the minimap2 error detection and matching script explained below. Takes two files with similar sequences and finds errors within the seqeunces as defined by minimap2. Information is parsed from minimap2 output and manipulated to determine where errors match. Each file is intended to be sequenced by different techniques (ex. Fasta file 1 is muscle and fasta file 2 is pbdagacon). Place all the files in the zip file into the same directory.
  
  Script Commands:
  
    FastaFile1: Sequences that need to be aligned in fasta format. Example: Muscle Consensus Sequence.
    FastaFile2: Sequences that need to be aligned in fasta format. Example: Pbdagcon Consensus Sequences.
    Reference: A reference sequence in fasta file format.
    OutputCsv: The name of the csv output.
    Minimap2: The location of the minimap2 application locally or on a server. Necessary only if minimap2 is not in the current directory.
  Example command:
  ```console
  foo@bar:~$ python MatchingMusclePbdagcon.py fastafile1.fasta fastafile2.fasta reference.fasta output.csv
  ```
 CSV file output:

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

The minimap2 output explained below for each fasta file is concatenated to these results.<br/> 
A Bed file is also outputted. The format for this file is explained here https://genome.ucsc.edu/FAQ/FAQformat.html#format1.
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
 CSV file output:
  Name: ```The name of the strand```<br/>
  Length: ```Total number of nucleotides```<br/>
  Reference Start Position: ```Actual start in a strand```<br/>
  Total Matching Range:	```Length minimap2 used in it's analysis (can be longer or shorter than length depending on reference)```<br/>
  Misaligned Front: ```Length minimap2 could not align with reference in the front (upstream)```<br/>
  Aligned: ```Number aligned correctly with reference for analysis```<br/>
  Misaligned Back: ```Length minimap2 could not align with reference in the back (downstream)```<br/>
  Subreads: ``` ********************* ```<br/>
  Mismatches: ```Number of errors in to5tal```<br/>
  Reference NonHP Insertion: ```number of non-homopolymers insertion errors```<br/>
  Reference NonHP Deletion: ```number of non-homopolymers deletion errors```<br/>
  Reference Substitution: ```number of substitution errors```<br/>
  Reference HP Insertion: ```number of homopolymers insertion errors```<br/>
  Reference HP Deletion: ```number of homopolymers deletion errors```<br/>
  Reference Long Insertion:  ```number of insertion errors longer than a length of 1```<br/>
  Reference Long Deletion: ```number of deletion errors longer than a length of 1```<br/>
  Consensus NonHP Insertion Positions: ```position of non-homopolymers insertion errors (start position is 1) (errors length of 1 and >1 are not seperated)```<br/>
  Consensus NonHP Deletion Positions:  ```position of non-homopolymers deletion errors```<br/>
  Consensus Substitution Positions:  ```position of substitution errors```<br/>
  Consensus HP Insertion Positions: ```position of homopolymers insertion errors```<br/>
  Consensus HP Deletion Positions: ```position of homopolymers deletion errors```<br/>
  Reference NonHP Insertion Positions:  ```position of non-homopolymers insertion errors (start position is based on reference)```<br/>
  Reference NonHP Deletion Positions: ```position of non-homopolymers deletion errors```<br/>
  Reference Substitution Positions: ```position of substitution errors```<br/>
  Reference HP Insertion Positions: ```position of homopolymers insertion errors```<br/>
  Reference HP Deletion Positions: ```position of homopolymers deletion errors```<br/>
  Reference Long Insertion Sites: ```position numbers of insertion errors longer than a length of 1```<br/>
  Reference Long Deletion Sites: ```position numbers of deletion errors longer than a length of 1```<br/>
<br/>
<br/>
<br/>
<br/>
  
  **Third File**: Matching Identification<br/>
Parses through two Minimap2 Error Detection Script's csv files explained above and determines where errors in a strand match. Validates and confirms stable modification sites.
  
  Script Commands:
  
    Minimap2ResultFile1: 1st file path for CSV with the minimap2 Results using Minimap2_Error_Detection.py.
    Minimap2ResultFile2: 2nd CSV file path using Minimap2_Error_Detection.py.
    OutputCsv: The matching sites from two differnt consensus generators compiled into a csv file.
  Example command:
  ```console
  foo@bar:~$ python Matching_Identification.py [csv file path 1] [csv file path 2] output.csv
  ```
 CSV file output: The output is the same as MatchingMusclePbdagcon.py (first script).
  ```
  
 
  
  
