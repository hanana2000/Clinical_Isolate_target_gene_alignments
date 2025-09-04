from Bio import SeqIO # type: ignore
import argparse
import os
import sys
from pathlib import Path
from collections import defaultdict
import subprocess
import pandas as pd # type: ignore
import gzip
from Bio.SeqRecord import SeqRecord # type: ignore
from Bio import AlignIO # type: ignore
import csv
from pathlib import Path


"""
using MAFFT align the target genes to the top hits found 
for each isolate

"""

# CLUSTAL similarity groups (classic definition)
STRONG_GROUPS = [
    set("STA"), set("NEQK"), set("NHQK"), set("NDEQ"),
    set("QHRK"), set("MILV"), set("MILF"), set("HY"), set("FYW")
]
WEAK_GROUPS = [
    set("CSA"), set("ATV"), set("SAG"), set("STNK"), set("STPA"),
    set("SGND"), set("SNDEQK"), set("NDEQHK"), set("NEQHRK"),
    set("FVLIM"), set("HFY")
]
UNKNOWN = set("XBZJ")


def align_tophits(target_genes1, target_genes2, path_tophits, isolate_proteins, output_folder):
    os.makedirs(output_folder, exist_ok = True)
    print("")

    print("\n" + "*" * 70)

    # iterate through first target gene dirs
    for dir1 in os.listdir(target_genes1):
        current_dir1, current_dir2 = os.path.join(target_genes1, dir1), ""
        output_gene_folder = os.path.join(output_folder, f"{dir1}_alignments")
        tmp_folder_path1, tmp_folder_path2 = os.path.join(output_gene_folder, "temp1"), os.path.join(output_gene_folder, "temp2")
        current_file1, current_file2 = "", ""
        os.makedirs(output_gene_folder, exist_ok=True)
        os.makedirs(tmp_folder_path1, exist_ok=True)
        os.makedirs(tmp_folder_path2, exist_ok=True)
        if not os.path.isdir(current_dir1): 
            print(f"\n[skip] {current_dir1} is not a dir\n")
            print("*" * 60)
            continue
        print("\nCurrently processing: " + dir1)
        # iterate through second target gene dir for this category 
        for dir2 in os.listdir(target_genes2): 
            # find the matching dirs
            current_dir1, current_dir2, current_file1, current_file2 = check_dir_match(dir1, dir2, current_dir1, current_dir2, current_file1, current_file2, target_genes2)
            if current_dir1 and  current_dir2 and current_file1 and current_file2: break
        # check if could not find matching files in both target gene dirs
        if not current_file1 or not current_file2: 
            print(f"Could not find matching target gene files for {current_dir1} and {current_dir2}\n")
            print("*" * 60)
            continue
        print(f"Found matching files {current_file1} and {current_file2} successfully!\n")
        # create records lists for all target genes in this category 
        targ_sequences1, targ_sequences2 = get_target_seqs(current_file1, current_file2)
        if not targ_sequences1 or not targ_sequences2: 
            print (f"Error: The file '{current_file1}' or '{current_file2}' was not found.")
            print("\n" + "*" * 70)
            continue
        
        print(f"\n{dir1} PAO1 library:")
        for item in targ_sequences1: print("\t" + item.description)
        print("")

        print(f"\n{dir1} PA14 library:")
        for item in targ_sequences2: print("\t" + item.description)
        print("")

        if len(targ_sequences1) != len(targ_sequences2): 
            print(f"[skip] {dir1} has different number of sequences between target sets\n")
            print("\n" + "*" * 70)
            continue
        
        # iterate through all isolate files in the current target gene category 
        category_align(dir1, path_tophits, targ_sequences1, targ_sequences2, tmp_folder_path1, tmp_folder_path2, isolate_proteins, output_gene_folder)
        print("\n" + "*" * 70)


def get_target_seqs(current_file1, current_file2): 
    targ_sequences1 = []
    targ_sequences2 = []
    try:
        with open(current_file1, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                targ_sequences1.append(record)
        with open(current_file2, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                targ_sequences2.append(record)
    except FileNotFoundError:
        print(f"Error: The file '{current_file1}' or '{current_file2}' was not found.")
        return [], []
    return targ_sequences1, targ_sequences2


def get_iso_prots_paths(root):
    """
    Find proteome files. Try flat search first; if none, fall back to recursive.
    Supports:
      - *.aa.gz
      - *.faa
      - *.fasta
      - NCBI Datasets: **/*protein.faa
    Returns: list[Path]

    """

    pats = [
        "*aa.gz", 
        "*.faa", 
        "*.fasta",
    ]
    
    for pat in pats: 
        flat = list(root.glob(f"*aa.gz"))  # get list with current isolate protein path
        if flat:
            print(f"\n[info] Found {len(flat)} {pat} in {root} (flat search).")
            print("\n" + "*-" * 40)
            return flat
    # if nothing found, try recursive search for NCBI datasets formatted proteomes
    if not flat: 
        print(f"Resorting to recursive search")
        flat = list(root.rglob(f"*.faa"))
        print(f"\n[info] Found {len(flat)} .faa in {root} (flat search).")
    return flat 


def get_iso_prots(iso_prots_paths, isolate_id): 
    iso_prots = []
    for file_path in iso_prots_paths: 
        file = os.path.basename(file_path)
        if isolate_id in file and file.endswith("aa.gz"):
            print(file)
            print(file_path)
            try: 
                with gzip.open(file_path, "rt", encoding="utf-8") as iso:
                    for record in SeqIO.parse(iso, "fasta"):
                        iso_prots.append(record)
            except FileNotFoundError:
                print(f"Error: The file '{file_path}' was not found.")
                return []
        elif isolate_id in str(file_path) and (file.endswith(".faa") or file.endswith(".fasta")): 
            print(f"ISOPROT_FILE: {file}")
            for record in SeqIO.parse(file_path, "fasta"):
                iso_prots.append(record)
    return iso_prots


def category_align(dir1, path_tophits, targ_sequences1, targ_sequences2, tmp_folder_path1, tmp_folder_path2, isolate_proteins, output_gene_folder):
    cat_summary = f"{output_gene_folder}/{dir1}_summary_all.txt"
    with open(cat_summary, "w") as summary: summary.write(f"{dir1} summary: \n\n")
    # tsv file per gene category
    root = Path(isolate_proteins)
    iso_prots_paths = get_iso_prots_paths(root) # get all isolate protein paths 
    gene_tsv1, gene_tsv2 =  f"{output_gene_folder}/{dir1}_PAO1_all.tsv", f"{output_gene_folder}/{dir1}_PA14_all.tsv"
    for isolate_file in os.listdir(path_tophits): 
        num_hits1, num_hits2, stats1, stats2 = 0, 0, [], []
        if isolate_file.endswith("canonical_crosscheck.tsv") and isolate_file.startswith(dir1): 
            isolate_path = os.path.join(path_tophits, isolate_file)
            isolate_id = isolate_file.split(".")[0].replace(f"{dir1}_","")
            current_output_folder1, current_output_folder2 = os.path.join(output_gene_folder, f"{isolate_id}_PAO1"), os.path.join(output_gene_folder, f"{isolate_id}_PA14")
            os.makedirs(current_output_folder1, exist_ok=True)
            os.makedirs(current_output_folder2, exist_ok=True)
            combined_file_path1, combined_file_path2 = os.path.join(current_output_folder1, f"{dir1}_combined.txt"), os.path.join(current_output_folder2, f"{dir1}_combined.txt")
            summary1, summary2 = f"{current_output_folder1}/summary_{isolate_id}_{dir1}.txt", f"{current_output_folder2}/summary_{isolate_id}_{dir1}.txt"
            print("\n" + isolate_id)
            print(isolate_file)
            # check if it is an empty file and skip 
            if os.path.getsize(isolate_path) < 2: 
                print(f"[skip] {isolate_file} is empty")
                write_to_cat_summary(num_hits1, num_hits2, cat_summary, isolate_id, targ_sequences1)
                continue 
            iso_df = pd.read_csv(isolate_path, sep='\t') # isolates top hits pandas dataframe 
            iso_prots = get_iso_prots(iso_prots_paths, isolate_id)
            if not iso_prots: 
                print(f"[skip] {isolate_id} has no proteins")
                continue 
            # align target genes from both sets to proteins of isolate
            num_hits1 = targ_iso_align(iso_prots, iso_df, isolate_id, summary1, dir1, combined_file_path1, targ_sequences1, tmp_folder_path1, current_output_folder1, gene_tsv1, gene_tsv2, 0)
            num_hits2 = targ_iso_align(iso_prots, iso_df, isolate_id, summary2, dir1, combined_file_path2, targ_sequences2, tmp_folder_path2, current_output_folder2, gene_tsv1, gene_tsv2, 1)
            # check if all genes were found and write to category summary file 
            write_to_cat_summary(num_hits1, num_hits2, cat_summary, isolate_id, targ_sequences1)
    print("\n")


def write_to_cat_summary(num_hits1, num_hits2, cat_summary, isolate_id, targ_sequences1): 
    if num_hits1 != num_hits2: 
        with open(cat_summary, "a") as summary: 
            summary.write(f"MAJOR ISSUE: {isolate_id} has mismatched number of hits between PAO1 and PA14")
    elif num_hits1 < len(targ_sequences1): 
        with open(cat_summary, "a") as summary: 
            summary.write(f"{isolate_id} is missing {len(targ_sequences1) - num_hits1} genes\n")


def targ_iso_align(iso_prots, iso_df, isolate_id, summary_path, dir1, combined_file_path, targ_sequences, tmp_folder_path, current_output_folder, gene_tsv1, gene_tsv2, col): 
    print(f"Number of isolate proteins: {len(iso_prots)}")
    # write/clear the summary file and combined alignment file
    with open(summary_path, "w") as summary: summary.write(f"{dir1} {isolate_id} summary: \n\n")
    with open(combined_file_path, "w") as combined: combined.write("")
    # create a folder for just the per gene alignments
    per_gene_aln = os.path.join(current_output_folder, "per_gene_alns")
    os.makedirs(per_gene_aln, exist_ok=True)
    # iterate through target genes and see if there is a match with isolate top hits 
    num_hits = 0 
    for i, targ in enumerate(targ_sequences):
        # select the isolate accession where the PAO1 id is the target gene
        match = iso_df.loc[iso_df.iloc[:, col] == targ.id, "sseqid"].tolist()
        if match: # if there is a top hit for the PAO1 protein 
            num_hits += 1 
            temp_file_path = f"{tmp_folder_path}/temp_{isolate_id}_{i}.fasta"
            rec1 = next((r for r in targ_sequences if r.id == targ.id), None)
            rec2 = next((r for r in iso_prots if r.id == match[0]), None)
            clustal_file = f"{per_gene_aln}/{os.path.basename(current_output_folder)}_gene_{i}.aln"
            fasta_file = f"{per_gene_aln}/{os.path.basename(current_output_folder)}_gene_{i}.fasta"
            with open(temp_file_path, "w") as f:
                SeqIO.write([rec1, rec2], f, "fasta")
            subprocess.run(["mafft", "--clustalout", "--auto", temp_file_path], stdout=open(clustal_file, "w"))
            subprocess.run(["mafft", "--auto", temp_file_path], stdout=open(fasta_file, "w"))
            diffs = get_diffs(clustal_file, fasta_file)
            with open(summary_path, "a") as summary: 
                for key in diffs: 
                    summary.write(f"{key}: {diffs[key]}\n")
                summary.write("\n")
            combine_alignments(clustal_file, combined_file_path)
            # write the current isolate+gene info to the common tsv file 
            write_tsv(gene_tsv1, gene_tsv2, diffs, isolate_id, col)
        else: # there is no PAO1 top hit 
            print(f"{targ.description} MISSING {current_output_folder.replace(f"{isolate_id}_","")}")
            continue 

    return num_hits


def write_tsv(gene_tsv1, gene_tsv2, stats, isolate_id, col):
    # reorder so description is at end of row and add isolate id at beginning 
    stats_ord = reorder_summary(stats, isolate_id)
    # append to the common tsv file based on whether PAO1 or PA14
    if not col: append_to_tsv(gene_tsv1, stats_ord)
    else: append_to_tsv(gene_tsv2, stats_ord)


def append_to_tsv(tsv_path, row_dict):
    write_header = not Path(tsv_path).exists()
    with open(tsv_path, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=row_dict.keys(), delimiter="\t")
        if write_header:
            writer.writeheader()
        writer.writerow(row_dict)


def reorder_summary(row_dict, isolate_id):
    # keys in the order you want, moving ref_desc and iso_desc to the end
    key_order = [k for k in row_dict.keys() if k not in ("ref_desc", "iso_desc")] \
                + ["ref_desc", "iso_desc"]
    # use dict unpacking to prepend isolate_id
    new_dict = {
        "isolate_id": isolate_id,
        **{k: row_dict[k] for k in key_order}
    }
    return new_dict


def combine_alignments(current_output_file, combined_file_path): 
    """
    write to a combined alignment file for this gene category 

    """
    with open(combined_file_path, "a") as combined:
        with open (current_output_file, "r") as alignment:  
            text = alignment.read()
            if not text.endswith("\n"):
                text += "\n\n"
            else: text += "\n"
            combined.write(text)


def get_diffs(aln_file, fasta_file):
    """
    get differences between the sequences with helper functions
    
    """
    # Read the alignment (CLUSTAL output from MAFFT)
    aln = AlignIO.read(aln_file, "clustal")
    assert len(aln) == 2, "Expected 2 sequences in alignment"
    s1, s2 = str(aln[0].seq), str(aln[1].seq)
    id1, id2 = aln[0].id, aln[1].id
    # print(f"s1: {s1}\ns2: {s2}")
    if len(s1) != len(s2): 
        raise ValueError("ALIGNMENTS DIFF LENGTHS\n")
    # iterate through each pair of aas in the sequences
    summary_lib =  check_pair_aa(s1, s2, id1, id2, fasta_file)
    return summary_lib
    

def check_pair_aa(s1, s2, id1, id2, fasta_file):
    """
    count number of gaps, conserved substitutions, weak substitutions, 
    identical aa's, complete mismatches, and unknown residues there are in
    a seq and return as a library

    """
    records = list(SeqIO.parse(fasta_file, "fasta"))
    s1 = s1.upper(); s2 = s2.upper()
    gaps1, gaps2 = s1.count('-'), s2.count('-') # find number of gaps in seq1 and seq2 
    unk_res1, unk_res2 = sum(c in UNKNOWN for c in s1), sum(c in UNKNOWN for c in s2) 
    cons_sub, weak_sub, identical, mismatch = 0, 0, 0, 0
    in_indel, indel_count = False, 0
    # walk the two sequences at the same time and find the type of sub per AA pair 
    for i, (a, b) in enumerate(zip(s1, s2), start=1):
        if a == '-' or b == '-': 
            in_indel = True 
            continue
        elif in_indel == True: 
            in_indel = False
            indel_count += 1
        if a == b: 
            identical += 1
            continue                    # B==B/J==J/X==X/N==N
        if a in UNKNOWN or b in UNKNOWN:
            mismatch += 1 
            continue                    # if is X or J or B and doesnt match (" ")                         
        if any(a in g and b in g for g in STRONG_GROUPS):
            cons_sub += 1  
            continue                    # conservative substitution (:)
        if any(a in g and b in g for g in WEAK_GROUPS):
            weak_sub += 1 
            continue                    # weak substitution (.)
        mismatch += 1 
    if in_indel:
        indel_count += 1

    return {
        "ref_desc": records[0].description, "iso_desc": records[1].description, 
        "id1": id1, "id2": id2, 
        "gaps_seq1": gaps1, "gaps_seq2": gaps2,
        "conservative_subs": cons_sub, 
        "weak_subs": weak_sub,
        "complete_mismatches": mismatch,
        "TOTAL_individual_aa_diffs": gaps1 + gaps2 + cons_sub + weak_sub + mismatch, 
        # "total_individual_aa_diffs_check": len(s1) - identical, 
        "indel_count": indel_count,
        "identical_matches": identical,
        "%_identity": identical / len(s1) * 100,
        "alignment_length": f"{len(s1)}", 
        "coverage_seq1": (len(s1) - gaps1) / len(s1),
        "coverage_seq2":  (len(s2) - gaps2) / len(s2)
        # "unknown_residues_seq1": unk_res1, "unknown_residues_seq2": unk_res2
    }


def check_dir_match(dir1, dir2, current_dir1, current_dir2, current_file1, current_file2, target_genes2): 
    """
    check if the directory names match and if there are matching 
    faa files in each dir
    
    """
    if dir1 == dir2: 
        current_dir2 = os.path.join(target_genes2, dir2)
        for file1 in os.listdir(current_dir1): 
            if file1.endswith(".faa"): 
                current_file1 = os.path.join(current_dir1, file1)
        for file2 in os.listdir(current_dir2): 
            if file2.endswith(".faa"):
                current_file2 = os.path.join(current_dir2, file2)
    return current_dir1, current_dir2, current_file1, current_file2



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Align isolate top hits from top_hits_isolates/compare_hits_between_refs.py to reference target genes in TWO SETS of target proteins.")

    parser.add_argument("--target_genes1", required=True, help="Path to the folder with subdirectories containing .faa files with first set of reference target genes.")
    parser.add_argument("--target_genes2", required=True, help="Path to the folder with subdirectories containing .faa files with second set of reference target genes.")
    parser.add_argument("--path_tophits", required=True, help="Path to top hits folder generated by compare_hits_between_refs.py, with canonical and unsure files.")
    parser.add_argument("--isolate_proteins", required=True, help="Path to the folder with all isolate proteins GCA.")
    parser.add_argument("--output_folder", required=True, help="Path to the output folder where results will be saved.")
    # if all 4 arguments are not provided, print help message
    if len(sys.argv) < 2:
        parser.print_help(sys.stderr)
        print("\n")
        sys.exit(1)
    args = parser.parse_args()

    align_tophits(args.target_genes1, args.target_genes2, args.path_tophits, args.isolate_proteins, args.output_folder)