from Bio import SeqIO # type: ignore
import argparse
import os
import sys
from pathlib import Path
from collections import defaultdict
import re
import subprocess
import pandas as pd

DELTA_BITS_THRESH = 20.0  

def get_top_hits(target_genes, isolates_path, output_folder): 
    # make the output directory if it does not exist
    output_folder = os.path.join(output_folder, "GCA", Path(target_genes).name)
    os.makedirs(output_folder, exist_ok=True)
    databases_folder = f"{output_folder}/databases_GCA/"
    os.makedirs(databases_folder, exist_ok=True)
    top_five_folder = f"{output_folder}/top_five_GCA/"
    os.makedirs(top_five_folder, exist_ok=True)

    # iterate through the isolates one at a time
    for iso_file in os.listdir(isolates_path):
        if iso_file.endswith(".gz"):
            print(f"Processing GCA file: {iso_file}")
        else: 
            print(f"skipping file {iso_file}")
            print("\n" + "#@*"* 35 + '\n')
            continue
        isolate_faa = os.path.join(isolates_path, iso_file)
        
        database_path = f"{databases_folder}{iso_file}_database.dmnd"
        command = ["diamond", "makedb", "--in", isolate_faa, "-d", database_path]
        print(">>", " ".join(command), "\n")
        try:
            result = subprocess.run(command, check=True)
            print(f"{iso_file} database created successfully!\n")
        except subprocess.CalledProcessError as e:
            print("There was an error creating database.\n")
            print(e)

        # Iterate through the target genes
        # and output top 5 results from diamond in tsv file  
        for target_folder in os.listdir(target_genes): 
            folder_path = os.path.join(target_genes, target_folder)
            if not os.path.isdir(folder_path):
                print(f"{folder_path} is not a directory")
                continue
            # iterate through categories of target genes
            for file in os.listdir(folder_path): 
                if file.endswith(".faa"): 
                    command = ["diamond", "blastp", "-d", database_path, "-q", f"{folder_path}/{file}", "-o", f"{top_five_folder}{target_folder}_{iso_file}.tsv", "--query-cover", "70", "--subject-cover", "70", "--max-target-seqs", "5", "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"]
                    print("\n>>", " ".join(command), "\n")
                    try: 
                        subprocess.run(command, check=True)
                        print(f"Category: {target_folder}, for isolate {iso_file} processed!")
                    except subprocess.CalledProcessError as e:
                        print("There was an error with diamond blastp.\n")
                        print(e)
        
        print("\n" + "#@*"* 35 + '\n')
    # check the difference in bitscore between the top and second best hit 
    get_top_hit(top_five_folder, output_folder)


def get_top_hit(diamond_hits, output_folder): 
    # make the output directory if it does not exist
    top_hits_dir = os.path.join(output_folder, "top_hits")
    top_hits_dir_consistent = os.path.join(output_folder, "top_hits_consistent")
    os.makedirs(top_hits_dir, exist_ok=True)
    os.makedirs(top_hits_dir_consistent, exist_ok=True)

    
    print("\n")
    for file in os.listdir(diamond_hits): 
        if file.endswith(".tsv"): 
            print(f"Processing file: {file}")
        else: 
            print(f"[skip]: {file}")
        print("\n")

        file_path = os.path.join(diamond_hits, file)

        base = Path(file).stem
        out_winners   = Path(top_hits_dir) / f"{base}.winners.tsv"
        out_ambiguous = Path(top_hits_dir) / f"{base}.ambiguous.tsv"

        out_winners_consistent = Path(top_hits_dir_consistent) / f"{base}.winners.tsv"
        out_ambiguous_consistent = Path(top_hits_dir_consistent) / f"{base}.ambiguous.tsv"

        
        if os.path.getsize(file_path) == 0:
            (Path(out_ambiguous_consistent)).write_text("")
            (Path(out_winners_consistent)).write_text("")
            print(f"[empty] {file_path}")
            print("\n" + "#@*"* 35 + '\n')
            continue

        ncols = len(pd.read_csv(file_path, sep="\t", nrows=1, header=None).columns)
        cols = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore",
        "qlen", "slen"
        ][:ncols]

        df = pd.read_csv(file_path, sep="\t", names=cols,
                         engine="c", comment="#")
        # sort so top row per qseqid has the highest bitscore
        df = df.sort_values(["qseqid", "bitscore"], ascending=[True, False])

        # top1 per query
        top1 = df.groupby("qseqid", as_index=False).nth(0).reset_index(drop=True)   
        
        # top2 per query (may be missing)
        top2 = df.groupby("qseqid", as_index=False).nth(1).reset_index(drop=True)
        # rename cols to avoid collision
        top2 = top2.rename(columns={c: f"{c}_2" for c in top2.columns if c != "qseqid"})

        # print(top1.head())
        # print(top2.head())

        # merge side-by-side to compute Δbits
        top = top1.merge(top2, on="qseqid", how="left")
        top["bitscore_2"] = top.get("bitscore_2", pd.Series([0.0]*len(top)))
        top["delta_bits"] = top["bitscore"] - top["bitscore_2"].fillna(0.0)

        # split by Δbits
        winners   = top[top["delta_bits"] >= DELTA_BITS_THRESH].copy()
        ambiguous = top[top["delta_bits"] <  DELTA_BITS_THRESH].copy()


        cols_win = [c for c in [
            "qseqid","sseqid","pident","length","qstart","qend","sstart","send",
            "evalue","bitscore","qlen","slen","query_cov","subject_cov","delta_bits"
        ] if c in winners.columns]
        if not winners.empty: 
            winners.to_csv(out_winners, sep="\t", index=False, columns=cols_win)
        winners.to_csv(out_winners_consistent, sep="\t", index=False, columns=cols_win)

        ### Ambiguous: stack top1 and top2 into rows
    
        # define a canonical set of "base" columns (for one hit)
        base_cols = [c for c in [
            "qseqid","sseqid","pident","length","qstart","qend",
            "sstart","send","evalue","bitscore","qlen","slen",
            "query_cov","subject_cov","delta_bits"
        ] if c in ambiguous.columns]
            # top1 rows
        amb1 = ambiguous.loc[:, base_cols].copy()
        amb1["rank"] = 1

        # build a mapping for top2 -> base col names
        map2 = {}
        for c in base_cols:
            c2 = (c + "_2") if c != "qseqid" else "qseqid"
            if c2 in ambiguous.columns:
                map2[c2] = c  # e.g., sseqid_2 -> sseqid

        # select and rename top2 to base col names
        amb2 = ambiguous.loc[:, list(map2.keys())].rename(columns=map2).copy()
        amb2["rank"] = 2

        # if some base columns don't exist for top2 (e.g., qlen/slen missing), add them
        for c in base_cols:
            if c not in amb2.columns:
                amb2[c] = pd.NA

        # ensure same column order and unique names
        amb1 = amb1.loc[:, base_cols + ["rank"]]
        amb2 = amb2.loc[:, base_cols + ["rank"]]

        # stack and sort
        amb_rows = pd.concat([amb1, amb2], ignore_index=True)
        amb_rows = amb_rows.sort_values(["qseqid","rank"]).reset_index(drop=True)

        if not ambiguous.empty: amb_rows.to_csv(out_ambiguous, sep="\t", index=False)
        amb1.to_csv(out_ambiguous_consistent, sep="\t", index=False)


        print("\n" + "#@*"* 35 + '\n')



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="use DIAMOND to find top hit proteins for each target gene from reference genomes")

    parser.add_argument("--target_genes", required=True, help="Path to the folder with subdirectories containing .faa files with reference target genes.")
    parser.add_argument("--isolates_path", required=True, help="Path to folder containing .faa file of all proteins")
    parser.add_argument("--output_folder", required=True, help="Path to the output folder where results will be saved.")
    # if all 4 arguments are not provided, print help message
    if len(sys.argv) < 2:
        parser.print_help(sys.stderr)
        print("\n")
        sys.exit(1)
    args = parser.parse_args()

    get_top_hits(args.target_genes, args.isolates_path, args.output_folder)