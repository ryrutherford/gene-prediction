import argparse
import json
from Bio import SeqIO
import numpy as np
from collections import deque
import math

states = ["Intergenic", "Genic Start", "Genic Middle", "Genic End"]
intergenic_index = 0
genic_start_index = 1
genic_middle_index = 2
genic_end_index = 3

#since log of 0 is -infinity we define our own log of 0 to be -1milli
def custom_log(number):
    if(number == 0):
        return -1000000
    else:
        return math.log(number)

def create_probability_tables(stats_dict):
    probabilities = {
        "Initial": {
            "Intergenic": 1,
            "Genic Middle": 0,
            "Genic Start": 0,
            "Genic End": 0,
        },
        "Transition": {
            "Intergenic": {"Intergenic": 0, "Genic Middle": 0, "Genic Start": 0, "Genic End": 0},
            "Genic Middle": {"Intergenic": 0, "Genic Middle": 0, "Genic Start": 0, "Genic End": 0},
            "Genic Start": {"Intergenic": 0, "Genic Middle": 0, "Genic Start": 0, "Genic End": 0},
            "Genic End": {"Intergenic": 0, "Genic Middle": 0, "Genic Start": 0, "Genic End": 0}
        },
        "Emission": {
            "Intergenic": {},
            "Genic Middle": {},
            "Genic Start": {},
            "Genic End": {}
        },
    }

    #emission probabilities
    probabilities["Emission"]["Intergenic"] = stats_dict["ig_nucleotides_freq"]
    probabilities["Emission"]["Genic Middle"] = stats_dict["g_codons_freq"]
    probabilities["Emission"]["Genic Start"] = stats_dict["codon_start_p"]
    probabilities["Emission"]["Genic End"] = stats_dict["codon_end_p"]

    #probabilities for transitioning from intergenic to next state
    average_length_ig = stats_dict["average_length_ig"]
    probabilities["Transition"]["Intergenic"]["Intergenic"] = (average_length_ig-1)/average_length_ig
    probabilities["Transition"]["Intergenic"]["Genic Start"] = 1/average_length_ig

    #probabilities for transitioning from genic start state to next state
    probabilities["Transition"]["Genic Start"]["Genic Middle"] = 1

    #probabilities for transitioning from genic end state to next state
    probabilities["Transition"]["Genic End"]["Intergenic"] = 1

    #probabilities for transitioning from genic middlestate to next state
    #we must divide the average_length_g by 3 because we emit codons not nucleotides
    average_length_g = int(stats_dict["average_length_g"]/3)
    probabilities["Transition"]["Genic Middle"]["Genic End"] = 1/average_length_g
    probabilities["Transition"]["Genic Middle"]["Genic Middle"] = (average_length_g - 1)/average_length_g
    
    return probabilities

def initial_state(arr, probabilities, first_nucleotide):
    for i in range(4):
        state = states[i]
        #the first entry is the value, the second is the traceback
        if first_nucleotide in probabilities["Emission"][state]:
            arr[i][0] = (custom_log(probabilities["Initial"][state])+custom_log(probabilities["Emission"][state][first_nucleotide])), (None, None)
        else:
            arr[i][0] = (custom_log(probabilities["Initial"][state]) + custom_log(0), (None, None))

def compute_trellis(trellis, probabilities, sequence, num_nucs):
    print("IN COMPUTE TRELLIS")
    #j tracks at which column we are in the trellis (index of the nucleotide in the sequence)
    j = 0
    for nucleotide in sequence:
        #we skip the first nucleotide since we've already computed it's probabilities
        if(j == 0):
            j += 1
            continue

        #computing the values for each of the 4 states
        for i in range(4):
            state = states[i]
            if(state == "Intergenic"):
                emissionProb = probabilities["Emission"][state][nucleotide]
            else:
                codon = None if j + 2 >= num_nucs else str(nucleotide + sequence[j+1] + sequence[j+2])
                if(codon in probabilities["Emission"][state]):
                    emissionProb = 0 if j + 2 >= num_nucs else probabilities["Emission"][state][codon]
                else:
                    emissionProb = 0

            #if i=0 then we're at the entry for the Intergenic state
            if(i == intergenic_index):
                #the transition from genic middle and genic start to intergenic is 0 so we only compare
                #the transitions from genic stop and intergenic to the intergenic state
                int_to_int = custom_log(emissionProb)+custom_log(probabilities["Transition"][state][state])+trellis[intergenic_index][j-1][0]
                #if j-3<0 then we haven't passed a codon yet so we ignore the transition from genic stop to intergenic
                stop_to_int = custom_log(0)+custom_log(0) if j - 3 < 0 else custom_log(emissionProb)+custom_log(probabilities["Transition"]["Genic End"][state])+trellis[genic_end_index][j-3][0]
                if(int_to_int >= stop_to_int):
                    trellis[i][j] = (int_to_int, (intergenic_index, j-1))
                else:
                    trellis[i][j] = (stop_to_int, (genic_end_index, j-3))
            elif(i == genic_start_index):
                #there is only one transition probability for the start state that is nonzero so this is the only one we need to check
                int_to_start = custom_log(emissionProb)+custom_log(probabilities["Transition"]["Intergenic"][state])+trellis[intergenic_index][j-1][0]
                trellis[i][j] = (int_to_start, (intergenic_index, j-1))
            elif(i == genic_middle_index):
                if(j - 3 < 0):
                    trellis[i][j] = (custom_log(0) + custom_log(0), (None, None))
                else:
                    start_to_genic = custom_log(emissionProb)+custom_log(probabilities["Transition"]["Genic Start"][state])+trellis[genic_start_index][j-3][0]
                    middle_to_middle = custom_log(emissionProb)+custom_log(probabilities["Transition"][state][state])+trellis[genic_middle_index][j-3][0]
                    if(start_to_genic >= middle_to_middle):
                        trellis[i][j] = (start_to_genic, (genic_start_index, j-3))
                    else:
                        trellis[i][j] = (middle_to_middle, (genic_middle_index, j-3))
            elif(i == genic_end_index):
                if(j - 3 < 0):
                    trellis[i][j] = (custom_log(0) + custom_log(0), (None, None))
                else:
                    #there is only one transition probability for the stop state that is nonzero so this is the only one we need to check
                    middle_to_end = custom_log(emissionProb)+custom_log(probabilities["Transition"]["Genic Middle"][state])+trellis[genic_middle_index][j-3][0]
                    trellis[i][j] = (middle_to_end, (genic_middle_index, j-3))      
        j += 1

def traceback(trellis, num_nucs):
    print("IN TRACEBACK")
    #we use a deque because it's O(1) for adding to the front
    seq = deque()

    #iterating from the end of the trellis
    column = num_nucs - 1
    #tracks which row, column gave the current value
    val_comes_from = (None, None)
    while column >= 0:
        if(column == num_nucs - 1):
            #getting the value of each of the 4 entries in column j
            values_in_col = [trellis[0][column][0], trellis[1][column][0], trellis[2][column][0], trellis[3][column][0]]
            #finding the row which has the max value in this column
            row = values_in_col.index(max(values_in_col))
            column = trellis[row][column][1][1]
        else:
            #storing the previous row because it's used at the end
            prev_row = row
            #getting the row which gave the current entry it's value
            row = val_comes_from[0]
            column = val_comes_from[1]
            if(row == None and column == None):
                row = prev_row
                if(row == intergenic_index):
                    seq.appendleft("I")
                else:
                    #print("State is G")
                    seq.appendleft("G")
                    seq.appendleft("G")
                    seq.appendleft("G")
                break

        val_comes_from = trellis[row][column][1]

        #if the row is an intergenic index we add 1 intergenic state, otherwise 3 genic states
        if(row == intergenic_index):
            seq.appendleft("I")
        else:
            seq.appendleft("G")
            seq.appendleft("G")
            seq.appendleft("G")
    return seq

#helper method to convert the genic_intergenic_sequence data into a gff3 file based on the info in the fasta_file
def output_to_file(genic_intergenic_sequence, fasta_file):
    print("OUTPUTTING TO FILE")
    #cds_regions stores all the cds_regions based on the seq_record from the fasta file
    cds_regions = {}
    #prev_state tracks the previous state value
    prev_state = None
    
    genic_intergenic_sequence = list(genic_intergenic_sequence)

    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        seq = seq_record.seq
        #getting the part of the genic_intergenic_sequence that matches this sequence record
        matching_seq = genic_intergenic_sequence[0:len(seq)]
        #deleting the matching part from the genic_intergenic_sequence list since we don't need it anymore
        del genic_intergenic_sequence[:len(seq)]
        index = 0
        for state in matching_seq:
            #if the current state is G and the previous state is I or none then we track the start index of this genic region
            if(state == "G" and (prev_state == "I" or prev_state == None)):
                start_of_state = 1 if index == 0 else index
                prev_state = "G"
            #if the current state is I and the previous state is G we have reached the end of the genic region so we add this cds to the dict
            if(state == "I" and prev_state == "G"):
                end_of_state = index - 1
                cds_regions[seq_record.id] = [(start_of_state, end_of_state)] if seq_record.id not in cds_regions else cds_regions[seq_record.id] + [(start_of_state, end_of_state)]
                start_of_state = None
                prev_state = "I"
            index += 1
        #if after iterating over the states for this seq_record we are still in a genic region
        #we say this genic region has ended at this point and add this region to the dict
        if(start_of_state != None):
            end_of_state = len(seq)
            cds_regions[seq_record.id] = [(start_of_state, end_of_state)] if seq_record.id not in cds_regions else cds_regions[seq_record.id] + [(start_of_state, end_of_state)]
            prev_state = None

    with open("viterbi_output.gff3", "w") as fp:
        fp.write("##gff-version 3\n")
        for key, value in cds_regions.items():
            for cds_region in value:
                start = cds_region[0]
                stop = cds_region[1]
                fp.write(f"{key}\tena\tCDS\t{start}\t{stop}\t.\t+\t0\t.\n")

#viterbi algorithm for 4 state HMM (where 3 states output 3 symbols and 1 outputs 1 symbol)
def viterbi(stats_dict, fasta_file):
    probabilities = create_probability_tables(stats_dict)
    sequence = []
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        sequence += list(seq_record.seq)
    num_nucs = len(sequence)
    trellis = np.zeros((4, num_nucs), dtype=tuple)
    #setting the initial values based on the initial probabilities
    initial_state(trellis, probabilities, sequence[0])
    compute_trellis(trellis, probabilities, sequence, num_nucs)
    genic_intergenic_sequence = traceback(trellis, num_nucs)
    print(f"Length of genic_intergenic sequence: {len(genic_intergenic_sequence)}, Number of nucleotides: {num_nucs}")
    output_to_file(genic_intergenic_sequence, fasta_file)
    
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("stats", help="The json file that gives information about the input fast sequence")
    parser.add_argument("fasta", help="The Fasta file that contains gene info")

    args = parser.parse_args()

    #path to the gff3 and fasta files
    stats = args.stats
    fasta = args.fasta

    with open(stats, "r") as fp:
        stats = json.load(fp)

    viterbi(stats, fasta)

if __name__ == '__main__':
    main()