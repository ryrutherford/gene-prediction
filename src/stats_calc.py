import os.path as osp
from Bio import SeqIO
import pandas as pd
import argparse
import json

directory = osp.normpath(osp.dirname(osp.dirname(__file__)))

gff3_file = osp.join(directory, "Vibrio_cholerae.GFC_11.37.gff3")
fasta_file = osp.join(directory, "Vibrio_cholerae.GFC_11.dna.nonchromosomal.fa")

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

def stats_calc(fasta_file, df):
    #store the nucleotide frequency in intergenic regions
    nuc_freq_ig = {
        "A": 0,
        "C": 0,
        "T": 0,
        "G": 0,
    }
    #total number of nucleotides in intergenic regions
    total_nuc_ig = 0

    #store the codon frequency in genic regions
    codon_freq = {}
    #total number of codons
    total_codons = 0

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

    g_test_len = 0
    g_test_av = 0

    #iterating over every sequence in the fasta file
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        position = 1
        #getting the subset of the dataframe that corresponds to this section of the fasta sequence
        df_seq = df.loc[df["seqid"] == seq_record.id]
        row_iterator = df_seq.iterrows()
        all_nucs_are_ig = False
        
        try:
            _, row = next(row_iterator)
            """
            g_test_av += int(row["end"]) + 1 - int(row["start"])
            g_test_len += 1
            """
        except StopIteration:
            all_nucs_are_ig = True

        
        codon = ""
        done = False
        phase = int(row["phase"])
        for nucleotide in seq_record.seq:
            #moving to the next cds region since we've passed the previous one
            #we only enter this condition if the row is not None
            if(done == False and all_nucs_are_ig == False):
                if(position > int(row["end"])):
                    try:
                        #print(f"{codon} + {phase}")
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
                        done = True

            #if we haven't entered the next cds yet we update the stats for the ig regions
            if(all_nucs_are_ig == True or done == True or position < int(row["start"]) or row["strand"] == "-"):
                #checking if we've entered a new intergenic region
                if(prev_region == None or prev_region == "Genic"):
                    total_ig += 1
                    prev_region = "Intergenic"
                nuc_freq_ig[nucleotide] += 1
                total_nuc_ig += 1
            #otherwise we are in a cds region
            elif(position >= int(row["start"])):
                #checking if we've entered a new genic region
                if(prev_region == None or prev_region == "Intergenic"):
                    total_g += 1
                    #we would need to account for overlapping regions here
                    prev_region = "Genic"
                #if the phase isn't 0 we need to skip phase nucleotides until it's 0
                if(phase != 0):
                    phase -= 1
                else:
                    codon += nucleotide
                    if(len(codon) == 3):
                        #update the codon frequency dictionary
                        if(codon in codon_freq):
                            codon_freq[codon] += 1
                        else:
                            codon_freq[codon] = 1
                        total_codons += 1
                        #reset the codon var to the empty string
                        codon = ""
                #for each nucleotide in the genic region we update the average length variable
                av_length_g += 1
            position += 1

    #Normalize averages
    for key in nuc_freq_ig.keys():
        nuc_freq_ig[key] = round(nuc_freq_ig[key]/total_nuc_ig, 4)

    for key in codon_freq.keys():
        codon_freq[key] = round(codon_freq[key]/total_codons, 4)

    av_length_ig = total_nuc_ig/total_ig
    l = av_length_g + total_nuc_ig
    num_g = av_length_g
    av_length_g = av_length_g/total_g

    final_json = {
        "ig_nucleotides_freq": nuc_freq_ig,
        "g_codons_freq": codon_freq,
        "average_length_ig": av_length_ig,
        "average_length_g": av_length_g
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
    #rint(total_nuc_ig)
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