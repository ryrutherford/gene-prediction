import os.path as osp
from Bio import SeqIO
import pandas as pd
import argparse
import json


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
                #WAS COMMENTED BEFORE CHANGES
                if(split_lines[6] == "-"):
                    continue

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

def get_genic_intergenic_regions(fasta_file, df):
    genic_regions = []
    intergenic_regions = []

    #iterating over every sequence in the fasta file
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        position = 0

        #getting the subset of the dataframe that corresponds to this section of the fasta sequence
        df_seq = df.loc[df["seqid"] == seq_record.id]
        for index, row in df_seq.iterrows():
            start = int(row["start"])
            end = int(row["end"])
            genic_regions.append(seq_record.seq[start - 1 : end])
            if(position < start - 1):
                intergenic_regions.append(seq_record.seq[position : start - 1])
            else:
                print(f"Position: {position} was greater than or equal to start - 1: {start -1}")
            position = end
        intergenic_regions.append(seq_record.seq[position:])
    return genic_regions, intergenic_regions

#helper function that gets the average length of a list of regions
def get_average_length_of_region(regions):
    lengths = [len(i) for i in regions]
    total_length = sum(lengths)
    total_regions = len(lengths)
    return total_length/total_regions

def calc_codon_stats(genic_regions):
    #variable to track start codons
    start_codon_freq = {}

    #store the codon frequency in genic regions (that aren't the end or start region)
    middle_codon_freq = {}

    #variable to track previous codon (used to check end codon probabilties)
    end_codon_freq = {}

    for region in genic_regions:
        #updating the frequencies for the middle_codon region 
        #(we will overcount by including start and end codons but account for this after iteration)
        for i in range(int(len(region)/3)):
            codon = str(region[i*3:(i+1)*3])
            middle_codon_freq[codon] = 1 if codon not in middle_codon_freq else middle_codon_freq[codon] + 1

        #adding the start codon to the start_codon_freq
        codon = str(region[0:3])
        start_codon_freq[codon] = 1 if codon not in start_codon_freq else start_codon_freq[codon] + 1
        #we need to remove a count from the middle codon_freq since we would have counted the start codon as a middle codon
        middle_codon_freq[codon] -= 1

        #adding the end codon to the end_codon_freq
        codon = str(region[-3:])
        end_codon_freq[codon] = 1 if codon not in end_codon_freq else end_codon_freq[codon] + 1
        #we need to remove a count from the middle codon_freq since we would have counted the end codon as a middle codon
        middle_codon_freq[codon] -= 1
    
    #getting the frequencies
    start_codon_freq = average_dict(start_codon_freq)
    middle_codon_freq = average_dict(middle_codon_freq)
    end_codon_freq = average_dict(end_codon_freq)
    return start_codon_freq, middle_codon_freq, end_codon_freq

def calc_nucleotide_stats(intergenic_regions):
    nucleotide_freq = {"A": 0, "C": 0, "G": 0, "T": 0}
    for region in intergenic_regions:
        for nucleotide in region:
            nucleotide_freq[str(nucleotide)] += 1
    nucleotide_freq = average_dict(nucleotide_freq)
    return nucleotide_freq

def compute_stats(fasta_file, df):
    #getting genic and intergenic regions based on gff3 file (in dataframe) and fasta file
    genic_regions, intergenic_regions = get_genic_intergenic_regions(fasta_file, df)

    #getting average region length for both genic and intergenic regions
    genic_region_av_length = int(round(get_average_length_of_region(genic_regions), 0))
    intergenic_region_av_length = int(round(get_average_length_of_region(intergenic_regions), 0))
    start_codon_freq, middle_codon_freq, end_codon_freq = calc_codon_stats(genic_regions)
    nucleotide_freq = calc_nucleotide_stats(intergenic_regions)

    final_json = {
        "average_length_ig": intergenic_region_av_length,
        "average_length_g": genic_region_av_length,
        "ig_nucleotides_freq": nucleotide_freq,
        "g_codons_freq": middle_codon_freq,
        "codon_start_p": start_codon_freq,
        "codon_end_p": end_codon_freq
    }

    with open("stats.json", "w") as fp:
        json.dump(final_json, fp, indent=2)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("gff3", help="The GFF3 file that accompanies the Fasta file")
    parser.add_argument("fasta", help="The Fasta file that contains gene info")

    args = parser.parse_args()

    #path to the gff3 and fasta files
    gff3 = args.gff3
    fasta = args.fasta

    df = extract_cds(gff3)
    #stats_calc(fasta, df)        
    compute_stats(fasta, df)

if __name__ == '__main__':
    main()