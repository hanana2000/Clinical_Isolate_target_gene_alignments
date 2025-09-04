# 🧬 Finding top hits for each target using diamond Blastp

requires: 
- DIAMOND 
- biopython 
- pandas

Logic Diagram: 
![gene search pipline diagram](https://github.com/hanana2000/Clinical_Isolate_target_gene_alignments/blob/501dad9a5ad1674b30afb49399bffdd9fe88ffe7/PA_Isolate_target_proteins.jpg)

## 🧫 Test Data

To test the scripts, test data is provided. 
download "fake_target_genesPAO1", "fake_target_genesPA14", "fake_isolate_proteins_GCA" and "fake_isolate_proteins_GCF". 

To test GCA data, run: 

```bash 
python3 find_top_hits_recurs.py --target_genes fake_target_genesPAO1 --isolates_path fake_isolate_proteins_GCA --output_folder fake_isolate_tophits

python3 find_top_hits_recurs.py --target_genes fake_target_genesPA14 --isolates_path fake_isolate_proteins_GCA --output_folder fake_isolate_tophits

python3 compare_hits_between_refs.py --target_genes1 fake_target_genesPAO1 --target_genes2 fake_target_genesPA14 --path_tophits1 fake_isolate_tophits/GCA/fake_target_genesPAO1 --path_tophits2 fake_isolate_tophits/GCA/fake_target_genesPA14 --output_folder fake_tophits_crosscheck_CGA

```

To test the GCF data, run: 

```bash
python3 find_top_hits_recurs.py --target_genes fake_target_genesPAO1 --isolates_path fake_isolate_proteins_GCF/isolates/ncbi_dataset/data/ --output_folder fake_isolate_tophits

python3 find_top_hits_recurs.py --target_genes fake_target_genesPA14 --isolates_path fake_isolate_proteins_GCF/isolates/ncbi_dataset/data/ --output_folder fake_isolate_tophits

python3 compare_hits_between_refs.py --target_genes1 fake_target_genesPAO1 --target_genes2 fake_target_genesPA14 --path_tophits1 fake_isolate_tophits/GCF/fake_target_genesPAO1 --path_tophits2 fake_isolate_tophits/GCF/fake_target_genesPA14 --output_folder fake_tophits_crosscheck_GCF

```


## 🔍 find_top_hits_recurs.py

First, either find_top_hits_recurs.py must be run on the isolate + target protein data. 

### 🧪 Usage of find_top_hits_recurs.py 

The find_top_hits_recurs.py file can take data formatted as the NCBI datasets command-line retrieval tool formats it, or as a flat directory. The datasets tool retreives the GCF data (RefSeq curated) with the directory structure as outlined below: 

```bash 
└── 📁isolates
    └── 📁ncbi_dataset
        └── 📁data
            └── 📁GCF_003968045.1
            │   ├── genomic.gff
            │   ├── protein.faa
            └── 📁GCF_003968115.1
            │   ├── genomic.gff
            │   ├── protein.faa
            ├── assembly_data_report.jsonl
            ├── dataset_catalog.json
    ├── md5sum.txt
    └── README.md

```

You want to pass the find_top_hits_recurs.py script the /ncbi_dataset/data folder. 

The find_top_hits_recurs.py script can also take GCA (or GCF) .fsa_aa.gz, .fasta, or .faa files that might have been manually downloaded. it expects flat data in this stucture: 

```bash
└── 📁isolate_proteins_GCA 
    ├── RXTE01P.1.fsa_aa.gz
    ├── RXTF01P.1.fsa_aa.gz
    ├── RXTH01P.1.fsa_aa.gz
    └── RXTQ01P.1.fsa_aa.gz
or 

└── 📁isolate_proteins_GCA 
    ├── RXTE01P.1.faa
    ├── RXTF01P.1.faa
    ├── RXTH01P.1.faa
    └── RXTQ01P.1.faa
or 

└── 📁isolate_proteins_GCA 
    ├── RXTE01P.1.fasta
    ├── RXTF01P.1.fasta
    ├── RXTH01P.1.fasta
    └── RXTQ01P.1.fasta
    
```

You want to pass the find_top_hits_recurs.py script the /isolate_proteins_GCA folder. 

Both folders expect the target genes to be formated in subdirectories representing categories of genes. For example: 

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

you would pass find_top_hits_recurs.py the PAO1_target_genes folder path. 

If you run find_top_hits.py without required arguments you will get this message: 

```bash
usage: find_top_hits.py [-h] --target_genes TARGET_GENES
                        --isolates_path ISOLATES_PATH
                        --output_folder OUTPUT_FOLDER

use DIAMOND to find top hit proteins for each target gene
from reference genomes

options:
  -h, --help            show this help message and exit
  --target_genes TARGET_GENES
                        Path to the folder with
                        subdirectories containing .faa files
                        with reference target genes.
  --isolates_path ISOLATES_PATH
                        Path to NCBI datasets generated
                        folder, usually called
                        'ncbi_dataset/data', containing each
                        clinical isolates subdirectory with
                        .faa file of all proteins
  --output_folder OUTPUT_FOLDER
                        Path to the output folder where
                        results will be saved.

```

If output folders specified do not exist, they will be created. 

### 📊 Output of find_top_hits_recurs.py

The output will be stored in the output folder specified. It will be stored in a subdir with the same dir name at the genome's proteome folder name with "results" appended at the front.  

the subdir will look like this: 

```bash 
└── 📁results_isolate_proteins_GCA
    └── 📁PA14_target_genes
        └── 📁databases_GCA
            ├── RXTE01P.1.fsa_aa.gz_database.dmnd
            ├── RXTF01P.1.fsa_aa.gz_database.dmnd
            ├── RXTH01P.1.fsa_aa.gz_database.dmnd
            ├── RXTQ01P.1.fsa_aa.gz_database.dmnd
        └── 📁top_five
            ├── efflux_RXTE01P.1.fsa_aa.gz.tsv
            ├── efflux_RXTF01P.1.fsa_aa.gz.tsv
            ├── efflux_RXTH01P.1.fsa_aa.gz.tsv
            ├── efflux_RXTQ01P.1.fsa_aa.gz.tsv
            ├── pili_RXTE01P.1.fsa_aa.gz.tsv
            ├── pili_RXTF01P.1.fsa_aa.gz.tsv
            ├── pili_RXTH01P.1.fsa_aa.gz.tsv
            ├── pili_RXTQ01P.1.fsa_aa.gz.tsv
            ├── porin_RXTE01P.1.fsa_aa.gz.tsv
            ├── porin_RXTF01P.1.fsa_aa.gz.tsv
            ├── porin_RXTH01P.1.fsa_aa.gz.tsv
            ├── porin_RXTQ01P.1.fsa_aa.gz.tsv
        └── 📁top_hits
            ├── efflux_RXTE01P.1.fsa_aa.gz.winners.tsv
            ├── efflux_RXTF01P.1.fsa_aa.gz.winners.tsv
            ├── efflux_RXTH01P.1.fsa_aa.gz.winners.tsv
            ├── efflux_RXTQ01P.1.fsa_aa.gz.winners.tsv
            ├── pili_RXTE01P.1.fsa_aa.gz.winners.tsv
            ├── pili_RXTF01P.1.fsa_aa.gz.winners.tsv
            ├── pili_RXTH01P.1.fsa_aa.gz.winners.tsv
            ├── pili_RXTQ01P.1.fsa_aa.gz.ambiguous.tsv
            ├── pili_RXTQ01P.1.fsa_aa.gz.winners.tsv
            ├── porin_RXTE01P.1.fsa_aa.gz.ambiguous.tsv
            ├── porin_RXTF01P.1.fsa_aa.gz.ambiguous.tsv
            ├── porin_RXTH01P.1.fsa_aa.gz.ambiguous.tsv
            ├── porin_RXTQ01P.1.fsa_aa.gz.ambiguous.tsv
        └── 📁top_hits_consistent
            ├── efflux_RXTE01P.1.fsa_aa.gz.ambiguous.tsv
            ├── efflux_RXTE01P.1.fsa_aa.gz.winners.tsv
            ├── efflux_RXTF01P.1.fsa_aa.gz.ambiguous.tsv
            ├── efflux_RXTF01P.1.fsa_aa.gz.winners.tsv
            ├── efflux_RXTH01P.1.fsa_aa.gz.ambiguous.tsv
            ├── efflux_RXTH01P.1.fsa_aa.gz.winners.tsv
            ├── efflux_RXTQ01P.1.fsa_aa.gz.ambiguous.tsv
            ├── efflux_RXTQ01P.1.fsa_aa.gz.winners.tsv
            ├── pili_RXTE01P.1.fsa_aa.gz.ambiguous.tsv
            ├── pili_RXTE01P.1.fsa_aa.gz.winners.tsv
            ├── pili_RXTF01P.1.fsa_aa.gz.ambiguous.tsv
            ├── pili_RXTF01P.1.fsa_aa.gz.winners.tsv
            ├── pili_RXTH01P.1.fsa_aa.gz.ambiguous.tsv
            ├── pili_RXTH01P.1.fsa_aa.gz.winners.tsv
            ├── pili_RXTQ01P.1.fsa_aa.gz.ambiguous.tsv
            ├── pili_RXTQ01P.1.fsa_aa.gz.winners.tsv
            ├── porin_RXTE01P.1.fsa_aa.gz.ambiguous.tsv
            ├── porin_RXTE01P.1.fsa_aa.gz.winners.tsv
            ├── porin_RXTF01P.1.fsa_aa.gz.ambiguous.tsv
            ├── porin_RXTF01P.1.fsa_aa.gz.winners.tsv
            ├── porin_RXTH01P.1.fsa_aa.gz.ambiguous.tsv
            ├── porin_RXTH01P.1.fsa_aa.gz.winners.tsv
            ├── porin_RXTQ01P.1.fsa_aa.gz.ambiguous.tsv
            └── porin_RXTQ01P.1.fsa_aa.gz.winners.tsv
    └── 📁PAO1_target_genes
        └── ...

```

The databases folder contains databases made by diamond to blastp against for all isolates. 

The top_five folder contains the top five hits for each target protein.

The top_hits folder contains the top hit for each protein if the bitscore threshhold was surpassed. if not, then it contains an "ambiguous" folder with the top two hits for that protein. It does not contain any empty files. 

The top_hits_consistent folder contains a winners and ambiguous file for each isolate, regardless of if the file is empty or not. Each file will only include one hit for each protein (ambiguous or not). 

If you chose to blastp against another set of target proteins, another folder (e.g. PAO1_target_genes) will be generated with the same subdirs under the GCA folder. 

If you use find_top_hits.py the main folder will be named GCF instead of GCA. 


## ⚖️ compare_hits_between_refs.py

If you have decided to run the isolates against two sets of target genes, then you can use the compare_hits_between_refs.py script to compare the top hits between either, and determine the top hit using the logic outlined in the repo diagram.  

### 🧪 Usage of compare_hits_between_refs.py

This script expects to be passed the same target genes folders as the previous script, one for each set of target genes blasted. Keep the order consisent (target_genes1 vs. target_genes2). 

This script also expects to be passed the folders generated by the previous script. Within the results_* folder, there will be a subdir for each set of target genes run. For example: 

```bash 
└── 📁results_isolate_proteins_GCA
    └── 📁PA14_target_genes
        └── 📁databases_GCA
        └── 📁top_five
        └── 📁top_hits
        └── 📁top_hits_consistent
    └── 📁PAO1_target_genes
        └── 📁databases_GCA
        └── 📁top_five
        └── 📁top_hits
        └── 📁top_hits_consistent

```

Pass the script the PA14_target_genes and the PAO1_target_genes folders. 

If you run this script without the required arguments, it will display this message: 

```bash 
usage: compare_hits_between_refs.py [-h] --target_genes1 TARGET_GENES1 --target_genes2 TARGET_GENES2
                                    --path_tophits1 PATH_TOPHITS1 --path_tophits2 PATH_TOPHITS2
                                    --output_folder OUTPUT_FOLDER

compare data from top hit proteins for each target gene from two reference genomes

options:
  -h, --help            show this help message and exit
  --target_genes1 TARGET_GENES1
                        Path to the folder with subdirectories containing .faa files with first set of
                        reference target genes.
  --target_genes2 TARGET_GENES2
                        Path to the folder with subdirectories containing .faa files with second set
                        of reference target genes. must be in SAME ORDER as first set.
  --path_tophits1 PATH_TOPHITS1
                        Path to top hits folder for reference set 1 (generated by
                        find_top_hits[_GCA].py)
  --path_tophits2 PATH_TOPHITS2
                        Path to top hits folder for reference set 2 (generated by
                        find_top_hits[_GCA].py)
  --output_folder OUTPUT_FOLDER
                        Path to the output folder where results will be saved.

```

If the output folder specified does not exist, it will be created. 

### 📊 Output of compare_hits_between_refs.py

The output will have the file structure below: 

```bash
└── 📁final_top_hit_crosscheck
    ├── efflux_RXTE01P.1.fsa_aa.gz.canonical_crosscheck.tsv
    ├── efflux_RXTE01P.1.fsa_aa.gz.unsure_crosscheck.tsv
    ├── efflux_RXTF01P.1.fsa_aa.gz.canonical_crosscheck.tsv
    ├── efflux_RXTF01P.1.fsa_aa.gz.unsure_crosscheck.tsv
    ├── pili_RXTE01P.1.fsa_aa.gz.canonical_crosscheck.tsv
    ├── pili_RXTE01P.1.fsa_aa.gz.unsure_crosscheck.tsv
    ├── pili_RXTF01P.1.fsa_aa.gz.canonical_crosscheck.tsv
    ├── pili_RXTF01P.1.fsa_aa.gz.unsure_crosscheck.tsv
    ├── porin_RXTE01P.1.fsa_aa.gz.canonical_crosscheck.tsv
    ├── porin_RXTE01P.1.fsa_aa.gz.unsure_crosscheck.tsv
    ├── porin_RXTF01P.1.fsa_aa.gz.canonical_crosscheck.tsv
    ├── porin_RXTF01P.1.fsa_aa.gz.unsure_crosscheck.tsv
    └── summary.txt

```

The canonical .tsv files will have the confirmed top hits (using the logic in the repo diagram) with the header: 

```bash
PAO1 qseqid	PA14 qseqid	sseqid	hit_type

```

PAO1 qseqid will have the id of the protein in the first reference genome. 
The PA14 qseqid will have the id of 0he protein in the second reference genome. 
The ssequid will have the id of the protein in the isolate itself. 
The "hit_type" will have whether it was a WIN or AMBIG in each reference genome, and whether it was NOT IN one or the other. 

The unsure .tsv files will have the ambiguous top hits (using the logic in the repo diagram) with the header: 

```bash 
PAO1 qseqid	PA14 qseqid	PAO1 sseqid	PA14 sseqid	hit_type

```

PAO1 qseqid will have the id of the protein in the first reference genome. 
The PA14 qseqid will have the id of the protein in the second reference genome. 
The PAO1 sseqid	PA14 sseqid will have the id of the protein in the isolate itself based on the first reference genome hit and the second reference genome hit (in case they did not match).
The "hit_type" will have whether it was a WIN or AMBIG in each reference genome, and whether it was NOT IN one or the other. 

The summary.txt file will contain info about the number of each type of gene, number of isolates tested, number of unsure protein hits, and number of missing hits (not found in the isolate for either set of target genes). For example: 

```bash
efflux genes: 10
pili genes: 6
porin genes: 1
number of isolates: 39
number of unsure hits: 17

efflux_RXWN01P.1.fsa_aa.gz. is missing 1 genes
porin_RXTE01P.1.fsa_aa.gz. is missing 1 genes
porin_RXTF01P.1.fsa_aa.gz. is missing 1 genes
porin_RXTH01P.1.fsa_aa.gz. is missing 1 genes
porin_RXTQ01P.1.fsa_aa.gz. is missing 1 genes

Unsure: porin_RXTE01P.1.fsa_aa.gz.
CAA78448.1, UXO62020.1, RTR86907.1, AMBIG PAO1, AMBIG PA14
Unsure: porin_RXTF01P.1.fsa_aa.gz.
CAA78448.1, UXO62020.1, RTR82200.1, AMBIG PAO1, AMBIG PA14
Unsure: porin_RXTH01P.1.fsa_aa.gz.
CAA78448.1, UXO62020.1, RTS25620.1, AMBIG PAO1, AMBIG PA14
Unsure: porin_RXTQ01P.1.fsa_aa.gz.
CAA78448.1, UXO62020.1, RTS78609.1, AMBIG PAO1, AMBIG PA14

```

# 🙋‍♀️ Author/ 📬 Contact

For questions or suggestions, contact: 

Hannah Kapoor
📧 hannahkapoor00@gmail.com 
