from Bio import SeqIO # type: ignore
import argparse
import os
import sys
from pathlib import Path
from collections import defaultdict
import re
import subprocess
import pandas as pd


def compare_hits(target_genes1, target_genes2, path_tophits1, path_tophits2, output_folder): 
    os.makedirs(output_folder, exist_ok=True)

    accession_lib = defaultdict(str)

    # iterate through all folders in the first provided target gene directory
    # and check that the same target folders/files exist in the second target gene dir 
    print("\n" + "#@" * 35 + "\n")

    # Iterate through folders in the provided target directory for target genes 
    # should be subdirs like "pili", "porin", "efflux"
    for target_folder1 in os.listdir(target_genes1): 
        target_folder_path1 = os.path.join(target_genes1, target_folder1)
        if not os.path.isdir(target_folder_path1): 
            print(f"{target_folder_path1} is not a dir")
            continue 
        if not target_folder1 in os.listdir(target_genes2): 
            print(f"{target_folder1} not found in {target_genes2}\n")
            print("#@" * 35 + "\n")
            continue 
        target_folder_path2 = os.path.join(target_genes2, target_folder1)
        if not os.path.isdir(target_folder_path2): 
            print(f"{target_folder_path2} is not a dir")
            continue 

        print(f"\ntarget folder 1: {target_folder_path1}")
        print(f"target folder 2: {target_folder_path2}\n")

        # iterate through all files in current sub dir and check for .faa 
        # e.g. "pili.faa", "efflux.faa", "porin.faa"
        for file1 in os.listdir(target_folder_path1): 
            if not file1.endswith(".faa"): 
                print(f"\n[skip] {file1}\n")
                continue
            if file1 not in os.listdir(target_folder_path2):
                print(f"[skip] {file1} not found in {target_folder_path2}\n")
                print("#@" * 35 + "\n")
                break

            file1_path = os.path.join(target_folder_path1, file1)
            file2_path = os.path.join(target_folder_path2, file1)
            
            accession_lib = populate_accessions(file1_path, file2_path, accession_lib)

        print("#@" * 35 + "\n")

    print("Full Accession Lib:")
    for id1, id2 in accession_lib.items(): 
            print(f"ID1 -> {id1}, ID2 -> {id2}")
    print("\n" + "#@" * 35 + "\n")

    # Now check that for winners, the top hits are the same, and output them in a single file
    print(f"Now crosschecking top hit files {path_tophits1} and {path_tophits2}\n")
    cross_check_top_hits(accession_lib, path_tophits1, path_tophits2, output_folder)


def cross_check_top_hits(accession_lib, path_tophits1, path_tophits2, output_folder):
    if not os.path.exists(path_tophits1) or not os.path.exists(path_tophits2): 
        print(f"One of the top hits dirs does not exist\n")
        return
    path_hits1 = os.path.join(path_tophits1, "top_hits_consistent")
    path_hits2 = os.path.join(path_tophits2, "top_hits_consistent")
    if not os.path.exists(path_hits1) or not os.path.exists(path_hits2): 
        print(f"One of the top hits subdirs doesn't exist: {path_hits1}, {path_hits2}\n")
    else: 
        print(f"Both top hits subirs found: {path_hits1}, {path_hits2}\n")
    # create a non-redundant dictionary of all isolate file prefixes
    target_iso_list1, target_iso_list2 = isolate_prefix_list_initial(path_hits1), isolate_prefix_list_initial(path_hits2)
    prefix_list, lone_prefixes = create_prefix_list(target_iso_list1, target_iso_list2)

    print("\nLone_prefixes only in one dir: ")
    for item in lone_prefixes: 
        print(item)

    # if the hit is a winner for both 
    # or if the hit is a winner for one but ambiguous for another, 
    # check if the top hits match and if they do then accept as canonical 
    for prefix in prefix_list:
        # iterate through one isolate at a time 
        winner_path1, winner_path2 = f"{path_hits1}/{prefix}winners.tsv", f"{path_hits2}/{prefix}winners.tsv"
        ambig_path1, ambig_path2 = f"{path_hits1}/{prefix}ambiguous.tsv", f"{path_hits2}/{prefix}ambiguous.tsv"
        if not os.path.exists(winner_path1) and not os.path.exists(winner_path2):
            print(f"ISSUE WITH {prefix} AMBIG")
            if not os.path.exists(ambig_path1) and not os.path.exists(ambig_path2): 
                print(f"ISSUE WITH {prefix} WINNERS")
        # print(f"\n{winner_path1}\n{winner_path2}\n{ambig_path1}\n{ambig_path2}")

        winPAO1, winPA14, ambPAO1, ambPA14 = load_top_hits(winner_path1), load_top_hits(winner_path2), load_top_hits(ambig_path1), load_top_hits(ambig_path2)
        # check for matches within current isolate hits 
        dfsure, dfunsure = check_for_matches(winPAO1, winPA14, ambPAO1, ambPA14, accession_lib, prefix)
        dfsure.to_csv(f"{output_folder}/{prefix}canonical_crosscheck.tsv", sep="\t", index=False)
        dfunsure.to_csv(f"{output_folder}/{prefix}unsure_crosscheck.tsv", sep="\t", index=False)
        print("\n" +"*" * 50 + "\n")



def check_for_matches(winPAO1, winPA14, ambPAO1, ambPA14, accession_lib, prefix):
    """
    accession lib: 
        PAO1 qseqid -> PA14 qseqid 
    winPAO1: 
        PAO1 qseqid -> sseqid 
    winPA14: 
        PA14 qseqid -> sseqid
    ambPAO1: 
        PAO1 qseqid -> sseqid
    ambPA14: 
        PA14 qseqid -> sseqid 

    """
    canonical, unsure, checked = [], [], []
    # list of all the qseqs with hits in PAO1 (win and ambig)
    combkeys_PAO1 = list(winPAO1.keys()) + list(ambPAO1.keys())
    # list of all the qseqs with hits in PA14 (win and ambig)
    combkeys_PA14 = list(winPA14.keys()) + list(ambPA14.keys())
    print(f"Now processing {prefix}: \n")
    # iterate through PAO1 keys first 
    for qseq1 in combkeys_PAO1: 
        qseq2 = accession_lib[qseq1]
        # if the gene found in both PAO1 and PA14
        if qseq2 in combkeys_PA14:
            # Both hits are winners, accept if they match 
            if qseq1 in winPAO1 and qseq2 in winPA14: 
                canonical, unsure = check_sseqid_match(qseq1, qseq2, winPAO1, winPA14, 
                                                       canonical, unsure, prefix, 
                                                       "WIN PAO1, WIN PA14")
                checked.append(qseq1)
                checked.append(qseq2)
            # in PAO1 the gene is a winner but in PA14 is ambiguous 
            elif qseq1 in winPAO1 and qseq2 in ambPA14:
                canonical, unsure = check_sseqid_match(qseq1, qseq2, winPAO1, ambPA14, 
                                                       canonical, unsure, prefix, 
                                                       "WIN PAO1, AMBIG PA14")
                checked.append(qseq1)
                checked.append(qseq2)
            # in PAO1 gene is ambiguous but in PA14 its a winner 
            elif qseq1 in ambPAO1 and qseq2 in winPA14: 
                canonical, unsure = check_sseqid_match(qseq1, qseq2, ambPAO1, winPA14, 
                                                       canonical, unsure, prefix, 
                                                       "AMBIG PAO1, WIN PA14")
                checked.append(qseq1)
                checked.append(qseq2)
            else: # if gene was ambiguous in both 
                print(f"Unsure: \n\t{qseq1} -> {ambPAO1[qseq1]}\n\t{qseq2} -> {ambPA14[qseq2]}")
                unsure.append({"PAO1 qseqid": qseq1, 
                                "PA14 qseqid": qseq2, 
                                "PAO1 sseqid": ambPAO1[qseq1],
                                "PA14 sseqid": ambPA14[qseq2],
                                "hit_type": "AMBIG PAO1, AMBIG PA14"})
                checked.append(qseq1)
                checked.append(qseq2)
        # if gene was only found in PAO1 as a hit 
        elif qseq2 not in combkeys_PA14:
            print(f"{qseq2} FOUND IN PAO1 but NOT IN PA14") 
            # if is a winner in PAO1 accept 
            win_or_ambig(qseq1, winPAO1, ambPAO1, qseq1, qseq2,
                          canonical, unsure, "PAO1", "PA14", prefix)
            checked.append(qseq1)
            checked.append(qseq2)
    # now iterate through PA14 genes to find hits 
    # that were found in PA14 but not PAO1
    inverted_accession_lib = {value:key for key, value in accession_lib.items()}
    for qseq2 in combkeys_PA14:
        qseq1 = inverted_accession_lib[qseq2]
        # if the qseqid has already been looked at in PAO1
        if qseq2 in checked: 
            continue 
        # if not
        else: 
            print(f"{qseq2} FOUND IN PA14 but NOT IN PAO1") 
            # if it is a PA14 winner accept it 
            win_or_ambig(qseq2, winPA14, ambPA14, qseq1, qseq2,
                          canonical, unsure, "PA14", "PAO1", prefix)
            checked.append(qseq1)
            checked.append(qseq2)

    dfcanon = pd.DataFrame(canonical)
    dfunsure = pd.DataFrame(unsure)
    return dfcanon, dfunsure


def win_or_ambig(checkqseq, winlib, ambiglib, qseq1, qseq2, canonical, unsure, in_genome, not_in_genome, prefix): 
    if checkqseq in winlib: 
        print(f"\t{qseq1} is a winner in {in_genome}!")
        canonical.append({"PAO1 qseqid": qseq1, 
                        "PA14 qseqid": qseq2, 
                        "sseqid": winlib[checkqseq],
                        "hit_type": f"WIN {in_genome} and NOT IN {not_in_genome}"})
    if checkqseq in ambiglib: 
        print(f"Unsure: \n\t{qseq1} -> {ambiglib[checkqseq]}\n")
        unsure.append({"PAO1 qseqid": qseq1, 
                        "PA14 qseqid": qseq2, 
                        "PAO1 sseqid": ambiglib[checkqseq] if ambiglib[qseq1] else "-",
                        "PA14 sseqid": ambiglib[checkqseq] if ambiglib[qseq2] else "-",
                        "hit_type": f"AMBIG {in_genome} and NOT IN {not_in_genome}"})


def win_or_ambig_PA14(checkqseq, winlib, ambiglib, qseq1, qseq2, canonical, unsure, in_genome, not_in_genome, prefix): 
    if checkqseq in winlib: 
        print(f"\t{qseq1} is a winner in {in_genome}!")
        canonical.append({"PAO1 qseqid": qseq1, 
                        "PA14 qseqid": qseq2, 
                        "sseqid": winlib[checkqseq],
                        "hit_type": f"WIN {in_genome} and NOT IN {not_in_genome}"})
    if checkqseq in ambiglib: 
        print(f"Unsure: \n\t{qseq1} -> {ambiglib[checkqseq]}\n")
        unsure.append({"PAO1 qseqid": qseq1, 
                        "PA14 qseqid": qseq2, 
                        "PAO1 sseqid": "-",
                        "PA14 sseqid": ambiglib[checkqseq],
                        "hit_type": f"AMBIG {in_genome} and NOT IN {not_in_genome}"})


def check_sseqid_match(qseq1, qseq2, lib1, lib2, canonical, unsure, prefix, message): 
    if lib1[qseq1] == lib2[qseq2]: 
        print(f"Found winners for {qseq1} and {qseq2}!")
        print (f"\t{lib1[qseq1]} matched {lib2[qseq2]}")
        canonical.append({"PAO1 qseqid": qseq1, 
                            "PA14 qseqid": qseq2, 
                            "sseqid": lib1[qseq1],
                            "hit_type": message})
    #sseqid does not match 
    else: 
        print(f"mismatch found for {qseq1} and {qseq2}:")
        print (f"\t{lib1[qseq1]} did not match {lib2[qseq2]}")
        unsure.append({"PAO1 qseqid": qseq1, 
                    "PA14 qseqid": qseq2, 
                    "PAO1 sseqid": lib1[qseq1],
                    "PA14 sseqid": lib2[qseq2],
                    "hit_type": message})
    return canonical, unsure


def isolate_prefix_list_initial(path_hits): 
    """
    get an inital list from a dir of all .tsv files
    no repeats in outputted list  

    """
    isolate_targ_list = []
    for file in os.listdir(path_hits): 
        if file.endswith(".tsv"):
            isolate_targ_list.append(file.replace("winners.tsv","").replace("ambiguous.tsv", ""))
    isolate_targ_list = list(set(isolate_targ_list))
    isolate_targ_list.sort()
    return isolate_targ_list


def create_prefix_list(target_iso_list1, target_iso_list2): 
    """
    compare lists and output a list with prefixes from both 
    and a list ith lone prefixes only in one list 
    
    """
    prefix_list, lone_prefixes = [], []
    combined = target_iso_list1 + target_iso_list2
    for prefix in combined: 
        if prefix in target_iso_list1 and prefix in target_iso_list2: 
            prefix_list.append(prefix)
        else: lone_prefixes.append(prefix)
    prefix_list, lone_prefixes = list(set(prefix_list)), list(set(lone_prefixes))
    prefix_list.sort()
    lone_prefixes.sort()
    return prefix_list, lone_prefixes


def load_top_hits(filepath): 
    """
    Reads .tsv file and returns a list of top hit accessions
    (assumes qseqid is column 1 and sseqid is column 2).
    no repeats allowed of qseqids (if ambig file shows top 2 for qsequid, only top kept)

    """
    hits = defaultdict(str)
    with open(filepath) as f:
        for line in f:
            if line.startswith("qseqid"):  # skip header
                continue
            parts = line.strip().split("\t")
            if len(parts) > 1:
                if not hits[parts[0]]: 
                    hits[parts[0]] = parts[1]  # qseqid -> sseqid
    return hits


def populate_accessions(file1_path, file2_path, accession_lib): 
    # read faa files from both directories 
    # and link record IDs by order (first one in file 1 linked with first in file 2)
    list1 = []
    list2 = []
    print(f"Now populating {file1_path} to {file2_path}")
    for record in SeqIO.parse(file1_path, "fasta"):
        list1.append(record.id)

    for record in SeqIO.parse(file2_path, "fasta"):
        list2.append(record.id)

    # link record IDs by order
    for id1, id2 in zip(list1, list2):
        accession_lib[id1] = (id2)

    return accession_lib


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="use DIAMOND to find top hit proteins for each target gene from reference genomes")

    parser.add_argument("--target_genes1", required=True, help="Path to the folder with subdirectories containing .faa files with first set of reference target genes.")
    parser.add_argument("--target_genes2", required=True, help="Path to the folder with subdirectories containing .faa files with second set of reference target genes. must be in SAME ORDER as first set.")
    parser.add_argument("--path_tophits1", required=True, help="Path to top hits folder for reference set 1 (generated by find_top_hits[_GCA].py)")
    parser.add_argument("--path_tophits2", required=True, help="Path to top hits folder for reference set 2 (generated by find_top_hits[_GCA].py)")
    parser.add_argument("--output_folder", required=True, help="Path to the output folder where results will be saved.")
    # if all 4 arguments are not provided, print help message
    if len(sys.argv) < 2:
        parser.print_help(sys.stderr)
        print("\n")
        sys.exit(1)
    args = parser.parse_args()

    compare_hits(args.target_genes1, args.target_genes2, args.path_tophits1, args.path_tophits2, args.output_folder)


