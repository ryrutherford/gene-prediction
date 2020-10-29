import os.path as osp
from Bio import SeqIO
import pandas as pd
import argparse
import json

directory = osp.normpath(osp.dirname(osp.dirname(__file__)))

gff3_file = osp.join(directory, "Vibrio_cholerae.GFC_11.37.gff3")
fasta_file = osp.join(directory, "Vibrio_cholerae.GFC_11.dna.nonchromosomal.fa")

#helper function to extract the cds regions from the gff3 file
def extract_cds(gff3_file):

    with open(gff3_file, "r") as gff3:
        lines = gff3.readlines()

        #if the gff3 file doesn't have the right header we throw an exception
        if(lines[0].startswith("##gff-version 3") == False):
            raise Exception("gff3 file missing header ##gff-version 3")

        gff3_dict = {
            "seqid": [],
            "source": [],
            "type": [],
            "start": [],
            "end": [], 
            "strand": [],
            "phase": [],
        }
        for line in lines:
            if("CDS" in line):
                split_lines = line.split("\t")

                #we don't care about negative strands
                #if(split_lines[6] == "-"):
                #    continue

                #adding the cds value to the dictionary
                index = 0
                for key in gff3_dict.keys():
                    if(key == "strand"):
                        index += 1
                    gff3_dict.setdefault(key, []).append(split_lines[index])
                    index += 1
    df = pd.DataFrame(data=gff3_dict)
    return df

def sum_dict(d):
    d_sum = 0
    for value in d.values():
        d_sum += value
    return d_sum

def average_dict(d):
    total = sum_dict(d)
    for key in d.keys():
        d[key] = d[key]/total
    return d

def average_dict2(d):
    for key in d.keys():
        d[key] = average_dict(d[key])
    return d
        

def stats_calc(fasta_file, df):
    #store the nucleotide frequency in intergenic regions
    nuc_freq_ig = {
        "A": 0,
        "C": 0,
        "T": 0,
        "G": 0,
    }
    #total number of nucleotides in intergenic regions
    #total_nuc_ig = 0

    #store the codon frequency in genic regions
    codon_freq = {}
    #total number of codons
    #total_codons = 0

    #average length of intergenic regions
    av_length_ig = 0
    #total number of intergenic regions
    total_ig = 0

    #average length of genic regions
    av_length_g = 0
    #total number of genic regions
    total_g = 0

    #variable to check if we've entered a new intergenic or genic region (if we have we will update some values)
    prev_region = None
    
    #variable to check if we are at a start codon
    is_start_codon = True
    start_codon_freq = {}

    #variable to track previous codon (used to check end codon probabilties)
    prev_codon = ""
    end_codon_freq = {}

    """
    #variable to track start nucleotides for intergenic regions
    start_nucleotide_freq = {
        "A": 0,
        "C": 0,
        "T": 0,
        "G": 0,
    }
    
    #variable to tract transition frequencies for intergenic regions
    ig_transitions = {
        "A": {},
        "C": {},
        "T": {},
        "G": {},
    }
    prev_nuc = ""

    #variable to track transition frequencies for genic regions
    g_transitions = {}
    prev_codon = ""
    """
    """
    g_test_len = 0
    g_test_av = 0
    """

    #iterating over every sequence in the fasta file
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        position = 1

        #getting the subset of the dataframe that corresponds to this section of the fasta sequence
        df_seq = df.loc[df["seqid"] == seq_record.id]

        row_iterator = df_seq.iterrows()
        all_nucs_are_ig = False
        
        try:
            #getting the first cds region for this sequence record
            _, row = next(row_iterator)
            """
            g_test_av += int(row["end"]) + 1 - int(row["start"])
            g_test_len += 1
            """
        except StopIteration:
            #if there are no cds's in this record then every nucleotide is in the intergenic region
            all_nucs_are_ig = True

        
        codon = ""
        done = False
        phase = int(row["phase"])
        for nucleotide in seq_record.seq:
            if(done == False and all_nucs_are_ig == False):
                #moving to the next cds region since we've passed the previous one
                if(position > int(row["end"])):
                    try:
                        done, row = next(row_iterator)
                        """
                        if(prev_region == "Genic" and position - 1 == row["start"]):
                            #do this if we want to count back to back genic regions as two genic regions
                            prev_region = "Intergenic"
                        """
                        """
                        if(row["strand"] != "-"):
                            g_test_av += int(row["end"]) + 1 - int(row["start"])
                            g_test_len += 1
                            #do this if we want to count back to back genic regions as two genic regions
                            prev_region = "Intergenic"
                        """
                        done = False
                        phase = int(row["phase"])
                    except StopIteration:
                        #when done is true we have read all the cds's for this seq record
                        done = True

            #if we haven't entered the next cds yet we update the stats for the ig regions
            if(all_nucs_are_ig == True or done == True or position < int(row["start"]) or row["strand"] == "-"):
                #checking if we've entered a new intergenic region
                if(prev_region == None or prev_region == "Genic"):
                    total_ig += 1
                    prev_region = "Intergenic"
                    #start_nucleotide_freq[nucleotide] += 1
                    #prev_nuc = ""
                """
                if(prev_nuc != ""):
                    #updating the transition probabilities
                    if(nucleotide in ig_transitions[prev_nuc]):
                        ig_transitions[prev_nuc][nucleotide] += 1
                    else:
                        ig_transitions[prev_nuc][nucleotide] = 1
                
                #updating the value for the previous nucleotide
                prev_nuc = nucleotide
                """
                #updating the nucleotide frequency table for intergenic regions
                nuc_freq_ig[nucleotide] += 1
                
                #total_nuc_ig += 1
            #otherwise we are in a cds region
            elif(position >= int(row["start"])):
                #checking if we've entered a new genic region
                if(prev_region == None or prev_region == "Intergenic"):
                    total_g += 1
                    #we would need to account for overlapping regions here
                    prev_region = "Genic"
                    is_start_codon = True

                    #updating the end frequency values
                    if(prev_codon != ""):
                        if(prev_codon not in end_codon_freq):
                            end_codon_freq[prev_codon] = 1
                        else:
                            end_codon_freq[prev_codon] += 1
                        #updating the general codon frequency since this codon shouldn't have counted for it
                        codon_freq[prev_codon] -= 1

                    #resetting the prev_codon value
                    prev_codon = ""

                #if the phase isn't 0 we need to skip phase nucleotides until it's 0
                if(phase != 0):
                    phase -= 1
                else:
                    codon += nucleotide
                    if(len(codon) == 3):
                    
                        """
                        if(prev_codon != ""):
                            #checking if the previous codon is already in our dictionary
                            if(prev_codon not in g_transitions):
                                g_transitions[prev_codon] = {}
                            #updating the value in our codon transitions table
                            if(codon in g_transitions[prev_codon]):
                                g_transitions[prev_codon][codon] += 1
                            else:
                                g_transitions[prev_codon][codon] = 1
                        """
                        #updating the prev_codon value for the next iteration
                        prev_codon = codon

                        #total_codons += 1
                        #if we are at a start codon we update the frequencies for this start codon
                        if(is_start_codon == True):
                            is_start_codon = False
                            if(codon in start_codon_freq):
                                start_codon_freq[codon] += 1
                            else:
                                start_codon_freq[codon] = 1
                        #otherwise we update the regular codon frequencies
                        else:
                            #update the codon frequency dictionary
                            if(codon in codon_freq):
                                codon_freq[codon] += 1
                            else:
                                codon_freq[codon] = 1
                        #reset the codon var to the empty string
                        codon = ""
                #for each nucleotide in the genic region we update the average length variable
                av_length_g += 1
            position += 1

    #after exiting the loop we need to add the last codon seen since it won't have been counted yet
    if(prev_codon != ""):
        if(prev_codon not in end_codon_freq):
            end_codon_freq[prev_codon] = 1
        else:
            end_codon_freq[prev_codon] += 1
        codon_freq[prev_codon] -= 1
    #Normalize averages
    total = sum_dict(nuc_freq_ig)
    av_length_ig = total/total_ig#total_nuc_ig/
    """
    for key in nuc_freq_ig.keys():
        nuc_freq_ig[key] = nuc_freq_ig[key]/total#round(nuc_freq_ig[key]/total_nuc_ig, 4)
    """
    
    nuc_freq_ig = average_dict(nuc_freq_ig)
    codon_freq = average_dict(codon_freq)
    #start_nucleotide_freq = average_dict(start_nucleotide_freq)
    start_codon_freq = average_dict(start_codon_freq)
    end_codon_freq = average_dict(end_codon_freq)
    #ig_transitions = average_dict2(ig_transitions)
    #g_transitions = average_dict2(g_transitions)
    
    #l = av_length_g + total_nuc_ig
    #num_g = av_length_g
    av_length_g = av_length_g/total_g

    final_json = {
        "ig_nucleotides_freq": nuc_freq_ig,
        "g_codons_freq": codon_freq,
        "average_length_ig": int(round(av_length_ig,0)),
        "average_length_g": int(round(av_length_g,0)),
        "codon_start_p": start_codon_freq,
        "codon_end_p": end_codon_freq
        #nucleotide_initial_p": start_nucleotide_freq,
        #"codon_transition_p": g_transitions,
        #"nucleotide_transition_p": ig_transitions
    }

    with open("stats.json", "w") as fp:
        json.dump(final_json, fp, indent=2)

    #print(json.dumps(nuc_freq_ig, indent=2))
    #print(json.dumps(codon_freq, indent=2))
    #print(av_length_ig)
    #rint(av_length_g)
    #print(g_test_av/g_test_len)
    #print(g_test_av)
    #print(g_test_len)
    #print(total_nuc_ig)
    #print(total_ig)
    #print(num_g)
    #print(total_g)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("gff3", help="The GFF3 file that accompanies the Fasta file")
    parser.add_argument("fasta", help="The Fasta file that contains gene info")

    args = parser.parse_args()

    #path to the gff3 and fasta files
    gff3 = args.gff3
    fasta = args.fasta

    df = extract_cds(gff3)
    stats_calc(fasta, df)        

if __name__ == '__main__':
    main()