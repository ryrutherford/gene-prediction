import pandas as pd
from stats_calc import extract_cds, average_dict
import argparse
import json

def compare_results(df_predicted, df_official):
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
                break
            if(row_official["start"] != row["start"] and row_official["end"] == row["end"]):
                predicted_stats["end_match_only"] += 1
                num_iters = -1
                break
        if(num_iters <= 0):
            predicted_stats["no_match"] += 1

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
                break
            #end match
            if(row_predicted["start"] != row["start"] and row_predicted["end"] == row["end"]):
                rp = int(row_predicted["start"])
                r = int(row["start"])
                if(rp - r >= -1 and rp - r <= 1):
                    print(f"Row predicted start: {rp}. Row start: {r}")
                official_stats["end_match_only"] += 1
                num_iters = -1
                break
        if(num_iters <= 0):
            official_stats["no_match"] += 1
    predicted_stats = average_dict(predicted_stats)
    official_stats = average_dict(official_stats)
    output_dict = {"predicted_stats": predicted_stats, "official_stats": official_stats}
    print(output_dict)
    with open("compare_stats.json", "w") as fp:
        json.dump(output_dict, fp, indent=2)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("gff3_predicted", help="The GFF3 file that was outputted by a program")
    parser.add_argument("gff3_official", help="The GFF3 file that comes from ftp")

    args = parser.parse_args()

    #path to the gff3 and fasta files
    gff3_predicted = args.gff3_predicted
    gff3_official = args.gff3_official

    df_predicted = extract_cds(gff3_predicted)
    df_official = extract_cds(gff3_official)
    compare_results(df_predicted, df_official)

if __name__ == '__main__':
    main()