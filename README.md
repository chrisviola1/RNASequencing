RNA Sequencing

These python scripts use various sequencing software to locate error types and combine mutiple sequencing techniques for increased performance.

**First File**: MatchingMusclePbdagcon.py<br/>
The source script. This uses both the minimap2 error detection and matching script explained below. Takes two files with similar sequences and finds errors within the seqeunces as defined by minimap2. Information is parsed from minimap2 output and manipulated to determine where errors match. Each file is intended to be sequenced by different techniques (ex. muscle or pbdagacon).
  
  Script Commands:
  
    FastaFile1: Sequences that need to be aligned in fasta format. Muscle Sequences.
    FastaFile2: Sequences that need to be aligned in fasta format. Pbdagcon Sequences.
    Reference: A reference sequence in fasta file format.
    OutputCsv: The name of the csv output.
    Minimap2: The location of the minimap2 application locally or on a server. Necessary only if minimap2 is not in the current directory.
  Example command:
  ```console
  foo@bar:~$ python MatchingMusclePbdagcon.py Sparc_HIV_ddAmp368Kb_RCA_reduced_no1.0nan_20.fasta Consensus_result_pbdagcon_pNLnefSBR2_ddAmp368Kb_RCA.fasta_penalty0
no1.0nan.fasta NL4.3_reference.txt trial_output.csv
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
  
  
