import argparse
import os
from Bio import SeqIO # type: ignore
import sys
import subprocess
from Bio import AlignIO # type: ignore


"""
alignment notes: 

* = identical across sequences
: = conserved substitution (similar amino acids)
. = weakly similar substitution

(blank) = mismatch

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


def mafft_align(target_genes1, target_genes2, output_folder):
    # make the output directory if it does not already exist
    os.makedirs(output_folder, exist_ok=True)
    print("\n" + "*" * 60)

    # iterate through first target gene dir 
    for dir1 in os.listdir(target_genes1):
        current_dir1, current_dir2 = os.path.join(target_genes1, dir1), ""
        output_gene_folder = os.path.join(output_folder, f"{dir1}")
        tmp_folder_path = os.path.join(output_gene_folder, "temp")
        current_file1, current_file2 = "", ""
        os.makedirs(output_gene_folder, exist_ok=True)
        os.makedirs(tmp_folder_path, exist_ok=True)
        if not os.path.isdir(current_dir1): 
            print(f"\n[skip] {current_dir1} is not a dir\n")
            print("*" * 60)
            continue
        print("\nCurrently processing: " + dir1)
        # iterate through second target gene dir 
        for dir2 in os.listdir(target_genes2): 
            # find the matching dirs
            current_dir1, current_dir2, current_file1, current_file2 = check_dir_match(dir1, dir2, current_dir1, current_dir2, current_file1, current_file2, target_genes2)
        # check if could not find matching files in both target gene dirs
        if not current_file1 or not current_file2: 
            print(f"Could not find matching target gene files for {current_dir1} and {current_dir2}\n")
            print("*" * 60)
            continue
        print(f"Found matching files {current_file1} and {current_file2} successfully!\n")
        
        # write/clear the summary file 
        with open(f"{output_gene_folder}/summary_{dir1}.txt", "w") as summary:
            summary.write(f"{dir1} summary: \n\n")
        # parse the current fasta target gene file for sequences 
        align_output(current_file1, current_file2, tmp_folder_path, output_gene_folder, dir1)

        print("*" * 60)

    print("")        


def align_output(current_file1, current_file2, tmp_folder_path, output_gene_folder, dir1): 
    """
    align the two current sequences and output the results to a clustal alignment fasta
    
    """
    combined_file_path = os.path.join(output_gene_folder, f"{dir1}_combined.txt")
    with open(combined_file_path, "w") as combined: combined.write("")
    records1 = list(SeqIO.parse(current_file1, "fasta"))
    records2 = list(SeqIO.parse(current_file2, "fasta"))
    if len(records1) != len(records2):
        raise ValueError("Files have different number of sequences!\n")
    for i, (rec1, rec2) in enumerate(zip(records1, records2), start=1):
        temp_file_path = f"{tmp_folder_path}/temp_{i}.fasta"
        with open(temp_file_path, "w") as f:
            SeqIO.write([rec1, rec2], f, "fasta")
        # mafft align the current sequence from the temp file 
        current_output_file = f"{output_gene_folder}/{dir1}_gene_{i}.aln"
        subprocess.run(["mafft", "--clustalout", "--auto", temp_file_path], stdout=open(current_output_file, "w"))
        # update to a summary file:
        diffs = get_diffs(current_output_file, rec1.description, rec2.description)
        with open(f"{output_gene_folder}/summary_{dir1}.txt", "a") as summary: 
            for key in diffs: 
                summary.write(f"{key}: {diffs[key]}\n")
            summary.write("\n")
        combine_alignments(current_output_file, combined_file_path)


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


def get_diffs(aln_file, desc1, desc2):
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
    summary_lib =  check_pair_aa(s1, s2, id1, id2, desc1, desc2)
    return summary_lib
    

def check_pair_aa(s1, s2, id1, id2, desc1, desc2):
    """
    count number of gaps, conserved substitutions, weak substitutions, 
    identical aa's, complete mismatches, and unknown residues there are in
    a seq and return as a library

    """
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
        "desc_PAO1": desc1, "desc_PA14": desc2, 
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
    parser = argparse.ArgumentParser(description="align target genes from two sets of reference data using mafft")

    parser.add_argument("--target_genes1", required=True, help="Path to the folder with subdirectories containing .faa files with first set of reference target genes.")
    parser.add_argument("--target_genes2", required=True, help="Path to the folder with subdirectories containing .faa files with second set of reference target genes. must be in SAME ORDER as first set.")
    parser.add_argument("--output_folder", required=True, help="Path to the output folder where results will be saved.")
    # if all 4 arguments are not provided, print help message
    if len(sys.argv) < 2:
        parser.print_help(sys.stderr)
        print("\n")
        sys.exit(1)
    args = parser.parse_args()

    mafft_align(args.target_genes1, args.target_genes2, args.output_folder)

