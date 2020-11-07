import pandas as pd
from stats_calc import extract_cds, average_dict
import argparse
import json
from Bio import SeqIO

def compare_results(df_predicted, df_official, fasta_file):
    predicted_stats ={
        "match": 0,
        "start_match_only": 0,
        "end_match_only": 0,
        "no_match": 0,
    }
    official_stats = {
        "match": 0,
        "start_match_only": 0,
        "end_match_only": 0,
        "no_match": 0,
    }
    #variables to track stats for mismatches (both predicted and official)
    missed_stats_start_pred = {}
    missed_stats_start_pred_len = [0, 0]
    missed_stats_end_pred = {}
    missed_stats_end_pred_len = [0, 0]
    missed_stats_both_pred = {}
    missed_stats_both_pred_len = [0, 0]
    missed_stats_start_off = {}
    missed_stats_start_off_len = [0, 0]
    missed_stats_end_off = {}
    missed_stats_end_off_len = [0, 0]
    missed_stats_both_off = {}
    missed_stats_both_off_len = [0, 0]

    #getting all the fasta sequences and their ids into a dictionary for easier access
    seq_record_dict = {}
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        seq_record_dict[seq_record.id] = seq_record.seq
    #calculating stats for predicted gff3 in a very inefficient manner
    for _, row in df_predicted.iterrows():
        df_for_seqid = df_official.loc[df_official["seqid"] == row["seqid"]]
        num_iters = 0
        for _, row_official in df_for_seqid.iterrows():
            num_iters += 1
            if(row_official["strand"] == "-"):
                continue
            if(row_official["start"] == row["start"] and row_official["end"] == row["end"]):
                predicted_stats["match"] += 1
                num_iters = -1
                break
            if(row_official["start"] == row["start"] and row_official["end"] != row["end"]):
                predicted_stats["start_match_only"] += 1
                num_iters = -1
                row_end = int(row["end"])
                codon = str(seq_record_dict[row["seqid"]][row_end - 3: row_end])
                #we add 1 to the tally of the end codon that was wrongly predicted
                missed_stats_end_pred[codon] = 1 if codon not in missed_stats_end_pred else missed_stats_end_pred[codon] + 1
                missed_stats_end_pred_len[0] += int(row["end"]) - int(row["start"])
                missed_stats_end_pred_len[1] += 1
                break
            if(row_official["start"] != row["start"] and row_official["end"] == row["end"]):
                predicted_stats["end_match_only"] += 1
                num_iters = -1
                row_start = int(row["start"])
                codon = str(seq_record_dict[row["seqid"]][row_start - 1: row_start + 2])
                #we add 1 to the tally of the start codon that was wrongly predicted
                missed_stats_start_pred[codon] = 1 if codon not in missed_stats_start_pred else missed_stats_start_pred[codon] + 1
                missed_stats_start_pred_len[0] += int(row["end"]) - int(row["start"])
                missed_stats_start_pred_len[1] += 1
                break
        if(num_iters != -1):
            predicted_stats["no_match"] += 1
            row_end = int(row["end"])
            codon_end = str(seq_record_dict[row["seqid"]][row_end - 3: row_end])
            row_start = int(row["start"])
            codon_start = str(seq_record_dict[row["seqid"]][row_start - 1: row_start + 2])
            if(len(codon_end) == 3 and len(codon_start) == 3):
                codon = f"{codon_start}, {codon_end}"
                missed_stats_both_pred[codon] = 1 if codon not in missed_stats_both_pred else missed_stats_both_pred[codon] + 1
            missed_stats_both_pred_len[0] += int(row["end"]) - int(row["start"])
            missed_stats_both_pred_len[1] += 1

    #calculating stats for official gff3 in a very inefficient manner
    for _, row in df_official.iterrows():
        df_for_seqid = df_predicted.loc[df_predicted["seqid"] == row["seqid"]]
        num_iters = 0
        for _, row_predicted in df_for_seqid.iterrows():
            num_iters += 1
            if(row["strand"] == "-"):
                continue
            #complete match
            if(row_predicted["start"] == row["start"] and row_predicted["end"] == row["end"]):
                official_stats["match"] += 1
                num_iters = -1
                break
            #start match
            if(row_predicted["start"] == row["start"] and row_predicted["end"] != row["end"]):
                official_stats["start_match_only"] += 1
                num_iters = -1
                row_end = int(row["end"])
                codon = str(seq_record_dict[row["seqid"]][row_end - 3: row_end])
                if(codon == "GTT"):
                    row_seqid = row["seqid"]
                    print(f"Non stop codon {codon} found in official gff3 at {row_seqid} end position {row_end}")
                #we add 1 to the tally of the end codon that was wrongly predicted
                missed_stats_end_off[codon] = 1 if codon not in missed_stats_end_off else missed_stats_end_off[codon] + 1
                missed_stats_end_off_len[0] += int(row["end"]) - int(row["start"])
                missed_stats_end_off_len[1] += 1
                break
            #end match
            if(row_predicted["start"] != row["start"] and row_predicted["end"] == row["end"]):
                official_stats["end_match_only"] += 1
                num_iters = -1
                row_start = int(row["start"])
                codon = str(seq_record_dict[row["seqid"]][row_start - 1: row_start + 2])
                if(codon not in ["ATG", "GTG", "TTG"]):
                    row_seqid = row["seqid"]
                    print(f"Non start codon {codon} found in official gff3 at {row_seqid} start position {row_start}")
                #we add 1 to the tally of the start codon that was wrongly predicted
                missed_stats_start_off[codon] = 1 if codon not in missed_stats_start_off else missed_stats_start_off[codon] + 1
                missed_stats_start_off_len[0] += int(row["end"]) - int(row["start"])
                missed_stats_start_off_len[1] += 1
                break
        if(num_iters != -1):
            official_stats["no_match"] += 1
            row_end = int(row["end"])
            codon_end = str(seq_record_dict[row["seqid"]][row_end - 3: row_end])
            if(len(codon_end) == 3):
                row_start = int(row["start"])
                codon_start = str(seq_record_dict[row["seqid"]][row_start - 1: row_start + 2])
                codon = f"{codon_start}, {codon_end}"
                missed_stats_both_off[codon] = 1 if codon not in missed_stats_both_off else missed_stats_both_off[codon] + 1
            missed_stats_both_off_len[0] += int(row["end"]) - int(row["start"])
            missed_stats_both_off_len[1] += 1

    predicted_stats = average_dict(predicted_stats)
    official_stats = average_dict(official_stats)
    missed_stats_start_pred = average_dict(missed_stats_start_pred)
    missed_stats_end_pred = average_dict(missed_stats_end_pred)
    missed_stats_both_pred = average_dict(missed_stats_both_pred)
    missed_stats_start_off = average_dict(missed_stats_start_off)
    missed_stats_end_off = average_dict(missed_stats_end_off)
    missed_stats_both_off = average_dict(missed_stats_both_off)
    output_dict = {
        "predictions/number of genes predicted": predicted_stats, 
        "predictions/number of actual genes": official_stats,
        "Gene length at predicted genes that didn't match": {"start only": missed_stats_start_pred_len[0]/missed_stats_start_pred_len[1], "end only": missed_stats_end_pred_len[0]/missed_stats_end_pred_len[1], "both": missed_stats_both_pred_len[0]/missed_stats_both_pred_len[1]},
        "Gene length at actual genes that didn't match": {"start only": missed_stats_start_off_len[0]/missed_stats_start_off_len[1], "end only": missed_stats_end_off_len[0]/missed_stats_end_off_len[1], "both": missed_stats_both_off_len[0]/missed_stats_both_off_len[1]},
        "Codon makeup at predicted genes that didn't match": {"start only": missed_stats_start_pred, "end only": missed_stats_end_pred ,"both": missed_stats_both_pred},
        "Codon makeup at actual genes that didn't match": {"start only": missed_stats_start_off, "end only": missed_stats_end_off, "both": missed_stats_both_off}
        }
    with open("compare_stats.json", "w") as fp:
        json.dump(output_dict, fp, indent=2)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("gff3_predicted", help="The GFF3 file that was outputted by a program")
    parser.add_argument("gff3_official", help="The GFF3 file that comes from ftp")
    parser.add_argument("fasta_file", help="The fasta file that accompanies the official gff3 file")

    args = parser.parse_args()

    #path to the gff3 and fasta files
    gff3_predicted = args.gff3_predicted
    gff3_official = args.gff3_official
    fasta_file = args.fasta_file

    df_predicted = extract_cds(gff3_predicted)
    df_official = extract_cds(gff3_official)
    compare_results(df_predicted, df_official, fasta_file)

if __name__ == '__main__':
    main()