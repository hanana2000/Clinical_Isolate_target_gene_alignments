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
TODO:

README usage

"""

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


def align_tophits(target_genes, path_tophits, isolate_proteins, output_folder):
    os.makedirs(output_folder, exist_ok = True)
    print("")

    print("\n" + "*" * 70)
    target_genome = os.path.basename(target_genes)
    print(f"Target genome: {target_genome}")

    # iterate through first target gene dirs
    for dir in os.listdir(target_genes):
        current_dir = os.path.join(target_genes, dir)
        output_gene_folder = os.path.join(output_folder, f"{dir}_alignments")
        tmp_folder_path = os.path.join(output_gene_folder, "temp")
        current_file = ""
        os.makedirs(output_gene_folder, exist_ok=True)
        os.makedirs(tmp_folder_path, exist_ok=True)
        if not os.path.isdir(current_dir): 
            print(f"\n[skip] {current_dir} is not a dir\n")
            print("*" * 60)
            continue
        print("\nCurrently processing: " + dir)
        # get the .faa file for this dir. there should only be one per category (efflux.faa) 
        for file in os.listdir(current_dir): 
            if file.endswith(".faa"): current_file = os.path.join(current_dir, file)
        if not current_file: 
            print(f"Could not find target gene files for {current_dir}\n")
            print("*" * 60)
            continue
        print(f"Found files {current_file} successfully!\n")
        # create records lists for all target genes in this category 
        targ_sequences = get_target_seqs(current_file)
        if not targ_sequences: 
            print (f"Error: The file '{current_file}' was not found.")
            print("\n" + "*" * 70)
            continue
        
        print(f"\n{dir} PAO1 library:")
        for item in targ_sequences: print("\t" + item.description)
        print("")


        category_align(target_genome, dir, path_tophits, targ_sequences, tmp_folder_path, isolate_proteins, output_gene_folder)
        print("\n" + "*" * 70)


def get_target_seqs(current_file): 
    targ_sequences = []
    try:
        with open(current_file, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                targ_sequences.append(record)
    except FileNotFoundError:
        print(f"Error: The file '{current_file}' was not found.")
        return []
    return targ_sequences


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


def category_align(target_genome, dir, path_tophits, targ_sequences, tmp_folder_path, isolate_proteins, output_gene_folder):
    cat_summary = f"{output_gene_folder}/{dir}_summary_all.txt"
    with open(cat_summary, "w") as summary: summary.write(f"{dir} summary: \n\n")
    # tsv file per gene category
    root = Path(isolate_proteins)
    iso_prots_paths = get_iso_prots_paths(root) # get all isolate protein paths 
    gene_tsv =  f"{output_gene_folder}/{dir}_{target_genome}_all.tsv"
    for dir_folder in os.listdir(path_tophits): 
        if dir_folder == "top_hits_consistent":
            for isolate_file in os.listdir(os.path.join(path_tophits, dir_folder)): 
                num_hits = 0
                if isolate_file.endswith("winners.tsv") and isolate_file.startswith(dir): 
                    isolate_path = os.path.join(path_tophits, dir_folder, isolate_file)
                    isolate_id = isolate_file.split(".")[0].replace(f"{dir}_","")
                    current_output_folder = os.path.join(output_gene_folder, f"{isolate_id}_{target_genome}")
                    os.makedirs(current_output_folder, exist_ok=True)
                    combined_file_path = os.path.join(current_output_folder, f"{dir}_combined.txt")
                    summary = f"{current_output_folder}/summary_{isolate_id}_{dir}.txt"
                    print("\n" + isolate_id)
                    print(isolate_file)
                    # check if it is an empty file and skip 
                    print(f"winners file size: {os.path.getsize(isolate_path)}\n")
                    df = pd.read_csv(isolate_path, sep="\t")
                    if df.empty:
                        print("\n[skip] File has only a header (or is empty).")
                        write_to_cat_summary(num_hits, cat_summary, isolate_id, targ_sequences)
                        continue 
                    iso_df = pd.read_csv(isolate_path, sep='\t') # isolates top hits pandas dataframe 
                    iso_prots = get_iso_prots(iso_prots_paths, isolate_id)
                    if not iso_prots: 
                        print(f"\n[skip] {isolate_id} has no proteins")
                        continue 
                    # align target genes from both sets to proteins of isolate
                    num_hits = targ_iso_align(iso_prots, iso_df, isolate_id, summary, dir, combined_file_path, targ_sequences, tmp_folder_path, current_output_folder, gene_tsv)
                    # check if all genes were found and write to category summary file 
                    write_to_cat_summary(num_hits, cat_summary, isolate_id, targ_sequences)
    print("\n")


def write_to_cat_summary(num_hits, cat_summary, isolate_id, targ_sequences): 
    # if num_hits1 != num_hits2: 
    #     with open(cat_summary, "a") as summary: 
    #         summary.write(f"MAJOR ISSUE: {isolate_id} has mismatched number of hits between PAO1 and PA14")
    if num_hits < len(targ_sequences): 
        with open(cat_summary, "a") as summary: 
            summary.write(f"{isolate_id} is missing {len(targ_sequences) - num_hits} genes\n")


def targ_iso_align(iso_prots, iso_df, isolate_id, summary_path, dir, combined_file_path, targ_sequences, tmp_folder_path, current_output_folder, gene_tsv): 
    print(f"Number of isolate proteins: {len(iso_prots)}")
    # write/clear the summary file and combined alignment file
    with open(summary_path, "w") as summary: summary.write(f"{dir} {isolate_id} summary: \n\n")
    with open(combined_file_path, "w") as combined: combined.write("")
    # create a folder for just the per gene alignments
    per_gene_aln = os.path.join(current_output_folder, "per_gene_alns")
    os.makedirs(per_gene_aln, exist_ok=True)
    # iterate through target genes and see if there is a match with isolate top hits 
    num_hits = 0 
    for i, targ in enumerate(targ_sequences):
        # select the isolate accession where the PAO1 id is the target gene
        match = iso_df.loc[iso_df.iloc[:, 0] == targ.id, "sseqid"].tolist()
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
            write_tsv(gene_tsv, diffs, isolate_id)
        else: # there is no PAO1 top hit 
            print(f"{targ.description} MISSING {current_output_folder.replace(f"{isolate_id}_","")}")
            continue 

    return num_hits


def write_tsv(gene_tsv, stats, isolate_id):
    # reorder so description is at end of row and add isolate id at beginning 
    stats_ord = reorder_summary(stats, isolate_id)
    # append to the common tsv file based on whether PAO1 or PA14
    append_to_tsv(gene_tsv, stats_ord)


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



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Align isolate top hits from top_hits_isolates/find_top_hits_recurs.py to reference target genes in ONE SET of target proteins.")

    parser.add_argument("--target_genes", required=True, help="Path to the folder with subdirectories containing .faa files with first set of reference target genes.")
    parser.add_argument("--path_tophits", required=True, help="Path to top hits folder generated by find_top_hits_recurs.py, with winners and ambiguous files.")
    parser.add_argument("--isolate_proteins", required=True, help="Path to the folder with all isolate proteins GCA.")
    parser.add_argument("--output_folder", required=True, help="Path to the output folder where results will be saved.")
    # if all 4 arguments are not provided, print help message
    if len(sys.argv) < 2:
        parser.print_help(sys.stderr)
        print("\n")
        sys.exit(1)
    args = parser.parse_args()

    align_tophits(args.target_genes, args.path_tophits, args.isolate_proteins, args.output_folder)