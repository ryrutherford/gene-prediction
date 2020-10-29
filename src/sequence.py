import argparse
import json
from Bio import SeqIO
import numpy as np

#the list of possible symbols that can be observed
symbols = ["A", "C", "G", "T", "AAT", "GCG", "TTC", "AAA", "TTA", "GTA", "CAT", "TCT", "AAG", "ATT", "AAC", "TTT",\
        "TTG", "TCG", "TCC", "ATC","GGT", "TGG", "GAC", "CAA", "GCC", "GTC", "ATG", "CCC", "AGT", "GGC", "GCT","GAA",\
            "CGC","GAG","CTT","GTG","CTA", "ACC","CAG","CCA","CTG","GAT","GCA","ACG","CTC","GTT","CCT","TGT","CAC","CGT",\
                "AGC","ACA","ACT","GGG","TAC","TAT","CCG","CGA","TCA","CGG","ATA","TAA","GGA","TGC","AGA","AGG","TGA","TAG"]
states = ["Intergenic", "Genic Start", "Genic Middle", "Genic End"]

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
            arr[i][0] = (probabilities["Initial"][state]*probabilities["Emission"][state][first_nucleotide], None)
        else:
            arr[i][0] = (0, None)

def viterbi(stats_dict, fasta_file):
    probabilities = create_probability_tables(stats_dict)
    sequence = []
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        sequence += list(seq_record.seq)
    arr = np.zeros((4, len(sequence)), dtype=tuple)
    initial_state(arr, probabilities, sequence[0])
    for nucleotide in sequence:
        pass
    
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


