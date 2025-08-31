# 🧬 Aligning two sets of target genes 

Aligns matched sets of target genes from two reference genomes and summarizes per-pair differences. The number of target genes per reference genome must be the same. 

### Why this exists

Given two matched target lists (e.g., P. aeruginosa PAO1 vs PA14), this tool:

1. builds pairwise MAFFT alignments (CLUSTAL format) per gene,
2. concatenates category-level alignments, and
3. outputs a concise per-pair summary of differences.

### Requirements

- Python ≥ 3.9
- MAFFT ≥ 7.5 (CLI on PATH)
- Biopython ≥ 1.8 (for Bio.AlignIO/Bio.SeqIO)

Install (example):

```bash 
# macOS (brew) or Linux (conda)
brew install mafft     # or: conda install -c bioconda mafft
python -m pip install biopython

```

## 🧫 Test Data

To test the script, test data is provided. 
download "fake_target_genesPAO1" and "fake_target_genesPA14". 

Quick start: 

```bash 
python3 mafft_target_genes.py --target_genes1 fake_target_genesPAO1 --target_genes2 fake_target_genesPA14 --output_folder fake_results

```

## 🧪 Usage of mafft_target_genes.py

This script will align two sets of target genes from different reference genomes. The number of target genes per ref genome must be the same. 

Both sets of reference genes must be formated in subdirectories representing categories of genes. Each .faa file can contain multiple sequences. For example: 

```bash
└── 📁PAO1_target_genes
    └── 📁efflux
        ├── efflux.faa
        ├── efflux.fna
    └── 📁pili
        ├── pili.faa
        ├── pili.fna
    └── 📁porin
        ├── porin.faa
        └── porin.fna

```

If you run mafft_target_genes.py without required arguments you will get this message: 

```bash 
usage: mafft_target_genes.py [-h] --target_genes1 TARGET_GENES1 --target_genes2 TARGET_GENES2
                             --output_folder OUTPUT_FOLDER

align target genes from two sets of reference data using mafft

options:
  -h, --help            show this help message and exit
  --target_genes1 TARGET_GENES1
                        Path to the folder with subdirectories containing .faa files with first set of reference
                        target genes.
  --target_genes2 TARGET_GENES2
                        Path to the folder with subdirectories containing .faa files with second set of reference
                        target genes. must be in SAME ORDER as first set.
  --output_folder OUTPUT_FOLDER
                        Path to the output folder where results will be saved.

```

## 📊 Output of mafft_target_genes.py

The output will be stored in the output folder specified. Each category of genes will be stored in a subdir with that category title. There will be temp folders for mafft to process which can be ignored. 

There will be one alignment file for each gene with the two reference genes aligned. There will be a combined.fasta with all the alignments combined. There will also be a summary file with an overview of differences found for each gene alignment. 

Toy data output: 

```bash 
└── 📁fake_results
    └── 📁fake_attack_genes
        ├── 📁temp
        ├── fake_attack_genes_combined.fasta
        ├── fake_attack_genes_gene_1.fasta
        ├── fake_attack_genes_gene_2.fasta
        ├── ...
        ├── summary_fake_attack_genes.txt
    └── 📁fake_colorful_genes
        ├── 📁temp
        ├── fake_colorful_genes_combined.fasta
        ├── fake_colorful_genes_gene_1.fasta
        ├── fake_colorful_genes_gene_2.fasta
        ├── ...
        ├── summary_fake_colorful_genes.txt
    └── 📁fake_defense_genes
        ├── 📁temp
        ├── fake_defense_genes_combined.fasta
        ├── fake_defense_genes_gene_1.fasta
        ├── fake_defense_genes_gene_2.fasta
        ├── fake_defense_genes_gene_3.fasta
        ├── ...
        └── summary_fake_defense_genes.txt

```

The alignments are made by MAFFT CLUSTAL align. 
Note: 

- \* = identical residues
- : = conserved substitution (strong similarity group)
- . = semi-conserved substitution (weak similarity group)
- [space] = no similarity symbol at that amino acid (non-similar or a gap)


Example alignment file: 

```bash 
CLUSTAL format alignment by MAFFT L-INS-i (v7.505)


BLUE456         GTNCIAAGCAFKLLPWWMWHARCYKDVEQFCFHMYCGEKMDCDKWFTEQVCDSRQHKDYP
BLUE123         GTNCIAAGCAFKLLPWWMWHARCYKDVEQFCFHMYCGEKMDCDKWFTEQVCDSRQHKDYP
                ************************************************************

BLUE456         EMSCFIHCSTVQLQYRPLDIKDSIDWSRELCFSFEHEDYTMLVIYTNHQVNMPPHEEMHK
BLUE123         EMSCFIHCSTVQLQYRPLDIKDSIDWSRELCFSFEHEDYTMLVIYTNHQVNMPPHEEMHK
                ************************************************************

BLUE456         KMTGQCVTTFHIGPMDHRTHHLGLRKHNNYLLEMRPSGMSKFSHIEQYTSHKKLHTGKGS
BLUE123         KMTGQCVTTFHIGPMDHRTHHLGLRKHNNYLLEMRPSGMSKFSHIEQYTSHKKLHTGKGS
                ************************************************************

BLUE456         WHSQSCEHCNKFVPPMDAFLYDIMRMYHDMTMETTERKDGLWSAMCYCVCEFNWVTMIYA
BLUE123         WHSQSCEHCNKFVPPMDAFLYDIMRMYHDMTMETTERKDGLWSAMCYCVCEFNWVTMIYA
                ************************************************************

BLUE456         GAHHTCVPPMPNWTVHEAMVLLQAKNKSDAPVMYTKDEQDTVDGIASHSVCAESYRLDTH
BLUE123         GAHHTCVPPMPNWTVHEAMVLLQAKNKSDAPVMYTKDEQDTVDGIASHSVCAESYRLDTH
                ************************************************************

BLUE456         GAHHTCVPPMPNW-------TVHEAMVLLQAKNKSDAPVM-YTKDEQDTVDGIASHSVCA
BLUE123         SQAKDQFKQYLNLGHLYGGVVNYYQRSLVMNQDHFRMPNMDYQRKVYPDVD---MHPECG
                .  :  .    *        . :    *:  :::   * * * :.    **    *. *.

BLUE456         ESYRLDTHHQLYQNKFGADRGRGVMSLYTPALLIPPQWWYYCWTLNYLLPDAFLEHKEID
BLUE123         P-----YDHQLYQNKFGADRGRGVMSLYTPALLIPPQWWYYCWTLNYLLPDAFLEHKEID
                       .****************************************************

BLUE456         VVYLCKKCHGNQRSCSGIMCSSRKYVKWKTGKIMAFCHQYEYARGNMWVGMMDGWDVVEI
BLUE123         VVYLCKKCHGNQRSCSGIMCSSRKYVKWKTGKIMAFCHQYEYARGNMWVGMMDGWDVVEI
                ************************************************************

BLUE456         KDNVSNHNLQLCSSKMFCWNVMCIDSRRDRHKNIFHLIVHQRNEYSAKGTCSFCWHTKAD
BLUE123         KDNVSNHNLQLCSSKMFCWNVMCIDSRRDRHKNIFHLIVHQRNEYSAKGTCSFCWHTKAD
                ************************************************************

BLUE456         CSWLNPALWEKKLLHIAF
BLUE123         CSWLNPALWEKKLLHIAF
                ******************

```

Similarity rule source: CLUSTAL similarity groups (not BLOSUM).
Strong groups: STA, NEQK, NHQK, NDEQ, QHRK, MILV, MILF, HY, FYW
Weak groups: CSA, ATV, SAG, STNK, STPA, SGND, SNDEQK, NDEQHK, NEQHRK, FVLIM, HFY

Note: 

- gaps_seq1 / gaps_seq2: number of "-" characters in each aligned sequence (PAO1, PA14)
- conservative_subs (:) / weak_subs (.): substitutions with a similar amino acid, or substitution with a weakly similar amino acid 
- complete_mismatches: residue–residue pairs with no CLUSTAL symbol (space), meaning complete amino acid mismatch excluding gaps
- Total individual aa differences:  conservative_subs + weak_subs + complete_mismatches + gaps in either sequence. This is individual amino acid differences, so it excludes indels (running gap count) 
- indel_count: number of gap runs (e.g. "---" counts as one event instead of 3)

Ambiguous residues policy: X/B/Z/J equal pairs (e.g., X==X) count as identity; non-equal pairs (e.g., B vs D) are treated as mismatches; they are not counted toward : or ..

Example summary file: 

```bash 
fake_colorful_genes summary: 

id1: RED456
id2: RED123
gaps_seq1: 26
gaps_seq2: 26
conservative_subs: 39
weak_subs: 18
complete_mismatches: 72
TOTAL_individual_aa_diffs: 181
indel_count: 12
identical_matches: 745

id1: BLUE456
id2: BLUE123
gaps_seq1: 8
gaps_seq2: 8
conservative_subs: 7
weak_subs: 7
complete_mismatches: 29
TOTAL_individual_aa_diffs: 59
indel_count: 4
identical_matches: 499

```

## Cite / Acknowledge

Please cite:

- Katoh & Standley (2013), MAFFT multiple sequence alignment software v7
- Thompson et al. (1994), CLUSTAL W/X (for similarity groups)


# 🙋‍♀️ Author/ 📬 Contact

For questions or suggestions, contact: 

Hannah Kapoor
📧 hannahkapoor00@gmail.com 