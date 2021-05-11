#! /usr/bin/env python
#packages that come with python
import pandas as pd
import string
import argparse
import os
import sys
import subprocess
import re

#Methods for Results

def cs_tags_str(mini):
    """Grabs the cs tag from the minimap2 output file

    Args:
        mini(sam file as a string): the minimap2 sam file

    Returns:
        Tags(string): identify summary inforamtion and are before the cigar string.
        cs(string): the cigar string identfying error types

    """
    tags = mini
    cs = mini
    
    tags = tags[tags.find('*')+1:]
    tagsCS = tags[tags.find('*'):]
    tags = tags[tags.find('*'):tags.find('cs')]
    
    cs = tagsCS[tagsCS.find('cs'):]
    cs = cs[cs.find('cs'):cs.find('rl:i:0')+6]
    return tags, cs

def primary_strand(mini, name):
    """Grabs the primary strand (first one listed) from the minimap2 sam file
        
    Args:
        mini(sam file as a string): the minimap2 sam file

    Returns:
        strand(string): Primary strand analyzed. This can be reverse complementary.

    """
    primary = False; info = mini; strand = ''; name = name[:-1]
    info = info[info.find('CL'):info.find('*')]
    num = re.split(r'\t+', info)[1]
    if int(num) == 0 or int(num) == 16:
        primary = True
    if primary == True:
        strand = mini[0:mini.find('rl:i:0')+6]
    else:
        strand = mini[0:mini.find('rl:i:0')+6]
        print("The strand aligned " + "(" + name + ")" + " is not a primary strand.")
    return strand

#For consective errors in a row within in the same error type
def consec_errors(mini,index):
    """The number of consective errors within a block of errors. A block of 
    errors is the number of errors of the same type in a row. For example 
    if one is -ctaaggttgc, the number returned is 10. 
        
    Args:
        mini(sam file as a string): the minimap2 sam file
        index(int): the position of the start of the block. Actual position not 
            based on the reference

    Returns:
        consec(int): The number of consective errors in a row.

    """
    consec = 0;index+=1
    let = mini[index]
    indentifers = ['-','+','=','*','\t']
    while let not in indentifers:
        consec+=1;index+=1
        let = mini[index]
    return consec

#Gets Substitution error information
#1.Number of Substitution errors
#2.Positions in Consensus and Ref
def sub(mini):
    """Subsitution data info when there is a subsitution error.
        
    Args:
        mini(sam file as a string): the minimap2 sam file

    Returns:
        0: the number of errors (int)
        1: the positions (actual) (string)
        2: the positions (based on reference) (string)

    """
    data,sub_pos,ref_sub_pos = [],[],[]
    sub_indexs = findOccurrences(mini,'*')
    for index in sub_indexs:
        shift = pos_shift(mini,index)
        sub_pos.append(index - shift)
        ref_shift = pos_shift_ref(mini,index)
        ref_sub_pos.append(index - ref_shift)
    data.append(len(sub_indexs)); data.append(sub_pos); data.append(ref_sub_pos)
    return data

def indels(mini, ref_start_pos, del_flag):
    """Outputs deletion information from minimap2 sam file depeding on the flag.

    Args:
        mini(sam file as a string): the minimap2 sam file
        ref_start_pos(integer): start position of testing strand in the reference strand
        
    Returns:
        data list with varying levels of information. Each spot has a desingated data point.
            0: number of non-homopolymers errors
            1: number of homopolymer errors 
            2: number of long deletion errors
            3: non_homopolymer consensus positions
            4: homopolymer consensus positions
                Example for formatting of positions below: 12*2^3 starts at reference position 12 followed by 
                2 non-hp errors followed by 3 hp errors.
            5: non_homopolymer reference positions
            6: homopolymer refererence positions
            7: long deletion positions
            
    """
    mini = mini.upper() 
    data, nonHP, HP, ref_nonHP, ref_HP, long_indel = [],[],[],[],[],[]
    indexs = []
    if del_flag == True:
        indexs = findOccurrences(mini,'-')
    else:
        indexs = findOccurrences(mini,'+')
    for index in indexs:
        consec = consec_errors(mini,index)
        indel_str = ""
        nonHp_str = ""; hp_str = ""
        for letter in range(consec):
            if letter == 0:
                first_letter_index = index
                first_indel_ref_pos = first_letter_index - pos_shift_ref(mini,index) + 1 #last actual char in ref 
                first_indel_ref_pos = first_indel_ref_pos + ref_start_pos #updating to actual position in alignment
                indel_str = str(first_indel_ref_pos)
            actual_letter_index = index+letter+1
            homopolymer = check_HP(mini, actual_letter_index, letter+1, consec, index)
            if homopolymer:
                shift = pos_shift(mini,index)
                HP.append(actual_letter_index - shift)
                indel_str = indel_str + '^'
            else:
                shift = pos_shift(mini,index)
                nonHP.append(actual_letter_index-shift)
                indel_str = indel_str + '*'
        indel_str, long_indel = indel_format(indel_str, long_indel)
        if indel_str.find('*') != -1 and indel_str not in long_indel:
            ref_nonHP.append(indel_str)
        if indel_str.find('^') != -1 and indel_str not in long_indel:
            ref_HP.append(indel_str)
    data.append(len(ref_nonHP)) #+ len([elt for elt in long_dels if "*" in elt]))
    data.append(len(ref_HP))  #+len([elt for elt in long_dels if "^" in elt]))
    data.append(len(long_indel))
    data.append(nonHP)
    data.append(HP)
    data.append(ref_nonHP)
    data.append(ref_HP)
    data.append(long_indel)
    return data

def indel_format(indel_str, long_indels):
    """Helper method to expand the number of error and hp/
        non-hp types visually. ^ is hp and * is non-hp.
        
    Args:
        indel_str(string): the block of same error type
        long_indels(list): list of error blocks with more than one error

    Returns:
        new_indel_str(string): new indel str with formatted errors
        long_indels(list): long indels errors and formatted appropriately

    """
    new_indel_str = ""; hp_str = ""; nonHp_str = ""; largeCount = 0
    for char in indel_str:
        if char.isdigit():
            new_indel_str = new_indel_str + char
        if char == "^":
            largeCount+=1
            if len(nonHp_str) > 0:
                new_indel_str = new_indel_str + "*" + str(len(nonHp_str))
                nonHp_str = ""
            hp_str = hp_str + "^"
        if char == "*":
            largeCount+=1
            if len(hp_str) > 0:
                new_indel_str = new_indel_str + "^" + str(len(hp_str))
                hp_str = ""
            nonHp_str = nonHp_str + "*"
    if len(nonHp_str) > 0:
        new_indel_str = new_indel_str + "*" + str(len(nonHp_str))
    if len(hp_str) > 0:
        new_indel_str = new_indel_str + "^" + str(len(hp_str))
    if largeCount > 1:
        long_indels.append(new_indel_str)
    return new_indel_str, long_indels

def check_HP(mini, actual_letter_index, pos_block, consec, block_start):
    """Check the tested nucleotide to see if it is homopolymer or not.
        
    Args:
        mini(sam file as a string): the minimap2 sam file
        actual_letter_index(int): index of letter based on position
        pos_block(int): the position in the block
        consec(int): the number in the block
        block_start(int): actual start of block of errors

    Returns:
        HP(boolean): shows whether the nucleotide is hp or not hp.

    """
    HP = False; HP_next = False; HP_prev = False
    
    next_letter, next_index  = next_let(mini,actual_letter_index) 
    pr_letter, pr_index = prev_let(mini,block_start,actual_letter_index)
    test_let = mini[actual_letter_index]
    
    if test_let == next_letter:
        HP_next = check_HPNext(mini, actual_letter_index, pos_block, consec)
    if test_let == pr_letter and HP_next == False:
        HP_prev = check_HPPrev(mini, actual_letter_index, pos_block, consec, block_start)
    if HP_next == True or HP_prev == True:
        HP = True
    return HP

def check_HPNext(mini, actual_letter_index, pos_block, consec):
    """Helper method to check the tested nucleotide of homopolymerness or not.
        Looks at nucleotides downstream of tested letters.
        
    Args:
        mini(sam file as a string): the minimap2 sam file
        actual_letter_index(int): index of letter based on position
        pos_block(int): the position in the block
        consec(int): the number in the block

    Returns:
        HP(boolean): shows whether the nucleotides downstream is hp or not hp.

    """
    HP = False
    
    next_letter, next_index  = next_let(mini,actual_letter_index) 
    test_let = mini[actual_letter_index]

    if (test_let == next_letter):
        if pos_block >= consec:
            if mini[actual_letter_index + 1] == '=':
                HP = True
            if mini[actual_letter_index + 1] == '*':
                HP = False
            #if mini [actual_letter_index + 1] == '+'
            #if mini [actual_letter_index + 1] == '-'
        else:
            HP = check_HPNext(mini, next_index, pos_block+1, consec)
    return HP

def check_HPPrev(mini, actual_letter_index, pos_block, consec, block_start):
    """Helper method to check the tested nucleotide of homopolymerness or not.
    Looks at nucleotides upstream of tested letters.
        
    Args:
        mini(sam file as a string): the minimap2 sam file
        actual_letter_index(int): index of letter based on position
        pos_block(int): the position in the block
        consec(int): the number in the block
        block_start(int): actual start of block of errors

    Returns:
        HP(boolean): shows whether the nucleotides upstream is hp or not hp.

    """
    HP = False
    
    pr_letter, pr_index = prev_let(mini,block_start,actual_letter_index)
    test_let = mini[actual_letter_index]
    
    if (test_let == pr_letter):
        if pos_block <= 1:
            if mini[last_type(mini, actual_letter_index, block_start)] == '=':
                HP = True
            if mini[last_type(mini, actual_letter_index,block_start)] == '*':
                HP = False
            #if mini [actual_letter_index + 1] == '+'
            #if mini [actual_letter_index + 1] == '-'
        else:
            HP = check_HPPrev(mini, pr_index, pos_block-1, consec, block_start)
    return HP
    
#gets mismatch tag
def mismatch(mini):
    """Identifies number of errors based on the NM tag
        
    Args:
        mini(sam file as a string): the minimap2 sam file

    Returns:
        num(int): integer of the number of errors in the aligned strand

    """
    NM = findOccurrences(mini,'N')
    if len(NM) == 0:
        num = 0
    else:
        NM = NM[0]
        num = ''
        index = NM + 5
        while mini[index] != '\t':
            if mini[index].isnumeric():
                num+=mini[index]
            index+=1
        num = int(num)
    return num

#used for identifying error positons; excludes extra symbols in cs tag
def pos_shift(mini,index):
    """Adjustment needed to transform the index to the actual position in the 
    aligned strand.
        
    Args:
        mini(sam file as a string): the minimap2 sam file
        index(int): start of a block for a error type

    Returns:
        pos_shift(int): integer of the amount to subtract to get the reference position.

    """
    smaller = mini[:index]
    smaller = smaller.replace('CS:Z:','')
    pos_shift = 5
    types = ['*','-','+','=']
    for char in smaller:
        if char in types:
            if char == '*':
                pos_shift+=2
            else:
                pos_shift+=1
    return pos_shift

def pos_shift_ref(mini,index):
    """Adjustment needed to transform the index to the actual position in the 
    aligned strand according to the reference.
        
    Args:
        mini(sam file as a string): the minimap2 sam file
        index(int): start of a block for a error type

    Returns:
        pos_shift(int): integer of the amount to subtract to get the reference position.

    """
    smaller = mini[:index]
    smaller = smaller.replace('CS:Z:','')
    smaller = smaller + '#' #indicate end of string from start to specific identifer
    pos_shift = 5
    for char in smaller:
        if char == '*':
            pos_shift = pos_shift + 2
        if char == '-' or char == '=':
            pos_shift = pos_shift + 1
    in_indexs = findOccurrences(smaller,'+')
    types = ['*','-','=','RL:I:0','#']
    for idx in in_indexs:
        block = idx
        while smaller[block] not in types:
            pos_shift+=1
            block+=1
    return pos_shift

#Gets alignment start position from reference fasta
def ref_start(mini):
    """Indentifies the reference start position using the CL tag.
    Args:
        mini(sam file as a string): the minimap2 sam file
    Returns:
        start_pos(int): an integer of the reference start position
    """
    info = mini
    occs = findOccurrences(info, '*')
    info = info[info.find('CL'):occs[1]]
    start_pos = int(re.split(r'\t+', info)[3])
    return start_pos

#finds occurences of a certain character
def findOccurrences(s, ch):
    """From a string, indentifies all index positions of a character.
    Args:
        s(string): some string to look at
        ch(string): the character to look for in the string, s.
    Returns:
        a list of index positions where ch shows up
    """
    return [i for i, letter in enumerate(s) if letter == ch]

#gathers the actual aligned sequence from the strand
def mini_seq(mini):
    """transforms back to the actual strand being aligned from the cs string.
    Args:
        mini(sam file as a string): the minimap2 sam file
    Returns:
        strand[0](string): the actual strand being aligned
    """
    chars = set('AGCT')
    strings = re.split(r'\t+', mini)
    strings = strings[1:]
    strand = [ele for ele in strings if all((c in chars) for c in ele)]
    return strand[0]

def lastSub(notMapped):
    """Helper function to get the number misaligned in the back
    Args:
        notMapped(list): list of mapped information split by misaligned areas
    Returns:
        subCount(int): number of misaligned in the back
    """
    revNotMapped = ''
    if len(notMapped) == 3:
        revNotMapped = notMapped[1][::-1]
    else:
        revNotMapped = notMapped[0][::-1]
    i = 0; subCount = ''
    while revNotMapped[i].isdigit():
        subCount+=revNotMapped[i];i+=1
    subCount = subCount[::-1]
    return(int(subCount))

def matching_area(mini):
    """Gets total matching area from "MDI" data string
    Args:
        mini(sam file as a string): the minimap2 sam file
    Returns:
        nums(list): misalignment information
            0: total matching area
            1: number of misaligned in the front of the strand
            2: number of aligned
            3: number of misaligned in the back of the strand
    """
    chars = set('MDIS')
    strings = re.split(r'\t+', mini) #splits data by tab
    #search for string with matching area data
    #gets strs that only contains M, D, or I
    potentialMatchData = [ele for ele in strings if any((c in chars) for c in ele)] 
    matchData = ""
    #gets last str that only contains digits and M, D, or I
    for dataPiece in potentialMatchData:
        if dataPiece[0].isdigit():
              matchData = dataPiece
    nums = [0,0,0,0]
    nums[0] = sum([int(s) for s in re.findall(r'\d+', matchData)]) #only digits
    if 'S' in matchData:
        notMapped = re.split(r'S', matchData) 
        if len(notMapped) == 3: #misaligned front and back
            subBack = lastSub(notMapped)
            nums[1] = int(notMapped[0])
            nums[2] = sum([int(s) for s in re.findall(r'\d+', notMapped[1])]) - subBack
            nums[3] = subBack
        else:
            if notMapped[1] == '': #misaligned back
                subBack = lastSub(notMapped)
                nums[2] = sum([int(s) for s in re.findall(r'\d+', notMapped[0])]) - subBack
                nums[3] = subBack
            else:
                nums[1] = int(notMapped[0])
                nums[2] = sum([int(s) for s in re.findall(r'\d+', notMapped[1])])
    else:
        nums[2] = sum([int(s) for s in re.findall(r'\d+', matchData)])
    return nums
    
def next_let(mini,next_let):
    """Finds the next letter depending on the error type

    Args:
        mini(string): the cs string from the minimap2 sam file
        next_let(int): the position of the letter being tested to get the next letter
        del_flag(binary): determine if looking for insertion or deletion errors

    Returns:
        let(string): the next letter

    """
    let = ''
    
    #hard cases when end of an error block
    if mini[next_let + 1] == '=':
        next_let+=2
        let = mini[next_let]
   
    elif mini[next_let + 1] == '-' or mini[next_let + 1] == '+':
        next_let+=2
        let = mini[next_let]

    elif mini[next_let + 1] == '*':
        next_let+=2
        let = mini[next_let]
    #when not end of a block
    else:
        #for del grab the next letter
        if mini[next_let + 1] != '\t':
            next_let+=1
            let = mini[next_let]
        #for in grab the next letter in next error block
        else:
            #block is number after the testing letter
            #testing letter index + number after in block + 2 for next letter
            block = consec_errors(mini,next_let)
            next_let = next_let+block+2
            let = mini[next_let+block+2]
    return let, next_let

def prev_let(mini,block_start,index):
    #mini: cs long str
    #index: location of error being compared
    """Finds the next letter depending on the error type

    Args:
        mini(string): the cs string from the minimap2 sam file
        block_start(int): the start position of error types. Ex: -ttggaac, - is the block start position.
        index(int): the position of the letter being tested
        del_flag(binary): determine if looking for insertion or deletion errors

    Returns:
        let(string): the previous letter

    """
    types = ['-','+','*','=']
    lastType = last_type(mini,index,block_start)
    if mini[lastType] == '=':
        if mini[index - 1] not in types:
            index = index - 1
            let = mini[index]
        else:
            index = block_start-1
            let = mini[index]

    elif mini[lastType] == '-' or mini[lastType] == '+': #goes to last indel error (+ga) -> prev letter is a
        #manage conidition eventually
        index = index - 1
        let = mini[index]

    elif mini[lastType] == '*':
        index = index - 2
        let = mini[index]
    else:
        index = index - 1
        let = mini[index]
    return let, index

def last_type(mini, index, block_start):
    """Indentifies the last error type in the block previous 
        to the one being tested.
    Args:
        mini(sam file as a string): the minimap2 sam file
        index(int): the position of the letter being tested
        block_start(int): the start position of error types
    Returns:
        int(last_type)(int): the position of the last error type in the cs string
    """
    last_type = 0
    types = ['-','+','*','=']
    smaller = mini[:block_start]
    last_type = 0
    for t in types:
        certain_type = 0
        occ = findOccurrences(smaller,t)
        if len(occ) == 0:
            certain_type = 0
        else:
            certain_type = max(occ)
        if certain_type > last_type:
            last_type = certain_type
    return int(last_type)

parser = argparse.ArgumentParser()
parser.add_argument("FastaFile", help="Sequences that need to be aligned in fasta format")
parser.add_argument("Reference", help="Reference to align sequences against")
parser.add_argument("OutputCsv", nargs='?', const = 1, default="Minimap2_Results.csv", help="The minimap2 error types compiled into a csv file.")
parser.add_argument("Minimap2", nargs='?', const = 1, default=os.getcwd() + "/minimap2", help="Minimap2 application. Assumes minimap2 is in the same directory. Otherwise, enter the directory where minimap2 is located.")

args = parser.parse_args()
fasta = args.FastaFile
ref = args.Reference
minimap2_app = args.Minimap2
result_name = args.OutputCsv

#Parsing Sequences
seqs = []
names = []
fasta = open(fasta)
for line in enumerate(fasta.readlines()):
    if line[0] % 2 == 0:#even line (id)
        seq = line[1]
        names.append(line[1])
    else:
        seq = seq + line[1]
        seqs.append(seq)
fasta.close()

erorr_sequences = []
name_i = 0
seq_count = 0; SEQS_LENGTH = len(seqs)
for seq in seqs:
    #Minimap2 alignment
    write = open("fasta_seq.fasta",'w')
    write.writelines(seq)
    write.close()
    minimap2 = subprocess.run(minimap2_app + " -a -L --cs=long --secondary=no" + " " + ref + " " + "fasta_seq.fasta", shell = True, capture_output = True)
    if minimap2.returncode != 0: #when theres an error
        print(minimap2.stderr)
    minimap2 = str(minimap2.stdout.decode())
    srs = len(minimap2.split("rl:i:0")) - 1
    #stream = os.popen(minimap2_app + " -a -L --cs=long" + " " + ref + " " + "fasta_seq.fasta", "r")
    #minimap2 = stream.read()
    #stream.close()

    #Gathering Results
    tags, cs = cs_tags_str(minimap2)
    minimap2 = primary_strand(minimap2, names[name_i])
    ref_start_pos = ref_start(minimap2)
    matchingInfo = matching_area(minimap2)

    #list of data -> 1.nonHp Count, 2. HP Count, 3. Con NonHP pos, 4. Con HP pos, 5. Ref nonHP pos, 6. Ref HP pos, 7. Long del sites
    ins_data = indels(cs, ref_start_pos - 1, False)
    dels_data = indels(cs, ref_start_pos - 1, True)
    subs_data = sub(cs)
    
    #updating ref positions to actual position for subs
    #dels_data[4] = [x + ref_start_pos for x in dels_data[4]] 
    #dels_data[5] = [x + ref_start_pos for x in dels_data[5]] 
    subs_data[2] = [x + ref_start_pos for x in subs_data[2]] 

    mini_info = {}
    mini_info["Name"] = names[name_i].replace("\n",'')
    mini_info["Length"] = len(mini_seq(minimap2))
    mini_info["Reference Start Position"] = ref_start_pos
    mini_info["Total Matching Range"] = matchingInfo[0]
    mini_info["Misaligned Front"] = matchingInfo[1]
    mini_info["Aligned"] = matchingInfo[2]
    mini_info["Misaligned Back"] = matchingInfo[3]
    mini_info["Subreads"] = srs
    mini_info["Mismatches"] = ins_data[0] + ins_data[1] + ins_data[2] + dels_data[0] + dels_data[1] + dels_data[2] + subs_data[0] 
    mini_info["Reference NonHP Insertion"] = ins_data[0] 
    mini_info["Reference NonHP Deletion"] = dels_data[0]
    mini_info["Reference Substitution"] = subs_data[0]
    mini_info["Reference HP Insertion"] = ins_data[1]
    mini_info["Reference HP Deletion"] = dels_data[1]
    mini_info["Reference Long Insertion"] = ins_data[2]
    mini_info["Reference Long Deletion"] = dels_data[2]
    mini_info["Consensus NonHP Insertion Positions"] = ins_data[3]
    mini_info["Consensus NonHP Deletion Positions"] = dels_data[3]
    mini_info["Consensus Substitution Positions"] = subs_data[1]
    mini_info["Consensus HP Insertion Positions"] = ins_data[4]
    mini_info["Consensus HP Deletion Positions"] = dels_data[4]
    mini_info["Reference NonHP Insertion Positions"] = ins_data[5]
    mini_info["Reference NonHP Deletion Positions"] = dels_data[5]
    mini_info["Reference Substitution Positions"] = subs_data[2]
    mini_info["Reference HP Insertion Positions"] = ins_data[6]
    mini_info["Reference HP Deletion Positions"] = dels_data[6]
    mini_info["Reference Long Insertion Sites"] = ins_data[7]
    mini_info["Reference Long Deletion Sites"] = dels_data[7]
    erorr_sequences.append(mini_info)
    name_i+=1
    
    seq_count+=1
    progress = (seq_count/SEQS_LENGTH)*100
    print('Progress [%.5f%%]\r'%progress, end="")
      
    
#Place into df
minimap2_results = pd.DataFrame.from_dict(erorr_sequences)
minimap2_results = minimap2_results[["Name", "Length", "Reference Start Position", "Total Matching Range",
                                     "Misaligned Front","Aligned","Misaligned Back", "Subreads",
                                     "Mismatches", "Reference NonHP Insertion","Reference NonHP Deletion","Reference Substitution",
                                     "Reference HP Insertion","Reference HP Deletion","Reference Long Insertion","Reference Long Deletion",
                                     "Consensus NonHP Insertion Positions","Consensus NonHP Deletion Positions",
                                     "Consensus Substitution Positions","Consensus HP Insertion Positions",
                                     "Consensus HP Deletion Positions",
                                    "Reference NonHP Insertion Positions", "Reference NonHP Deletion Positions",
                                    "Reference Substitution Positions", "Reference HP Insertion Positions", 
                                     "Reference HP Deletion Positions", "Reference Long Insertion Sites", "Reference Long Deletion Sites"]]
minimap2_results.set_index('Name')
minimap2_results.to_csv(result_name)
print(minimap2_results)