# 🧬 Finding top hits for each target using diamond Blastp

requires: 
- DIAMOND 
- biopython 
- pandas

Logic Diagram: 
![gene search pipline diagram](https://github.com/hanana2000/Clinical_Isolate_target_gene_alignments/blob/f0dfb14a9100767b402b663f975e6a75a8342d3a/PA_Isolate_target_proteins.jpg)

## Test Data




## 🔍 find_top_hits.py or find_top_hits_GCA

First, either find_top_hits.py or find_top_hits_GCA.py must be run on the isolate/ target protein data. 

### 🧪 Usage of find_top_hits.py or find_top_hits_GCA.py

The find_top_hits.py file and find_top_hits_GCA.py file will both give the same output. they both find the top 5 hits, and from this the top hits (using logic outlined in diagram found in this repo based on bitscore), including top hit ambiguity. 

the only difference in these files is the required input format. 

The find_top_hits.py file expects data formatted as the NCBI datasets command-line retrieval tool formats it. The datasets tool retreives the GCF data (RefSeq curated), but it should be in the same format (.faa files) as the GCA data. The directory structure is outlined below: 

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

You want to pass the find_top_hits.py script the /ncbi_dataset/data folder. 

The find_top_hits_GCA.py script can take GCA (or GCF) .fsa_aa.gz files that might have been manually downloaded. it expects data in this stucture: 

```bash
└── 📁isolate_proteins_GCA
    ├── RXTE01P.1.fsa_aa.gz
    ├── RXTF01P.1.fsa_aa.gz
    ├── RXTH01P.1.fsa_aa.gz
    └── RXTQ01P.1.fsa_aa.gz
    
```

You want to pass the find_top_hits_GCA.py script the /isolate_proteins_GCA folder. 

both folders expect the target genes to be formated in subdirectories representing categories of genes. For example: 

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

you would pass either find_top_hits.py or find_top_hits_GCA.py the PAO1_target_genes folder path. 

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

If you run find_top_hits_GCA.py without required arguments you will get this message: 

```bash 
usage: find_top_hits_GCA.py [-h] --target_genes TARGET_GENES
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
                        Path to folder containing .faa file
                        of all proteins
  --output_folder OUTPUT_FOLDER
                        Path to the output folder where
                        results will be saved.

```

If output folders specified do not exist, they will be created. 

### 📊 Output of find_top_hits.py and find_top_hits_GCA.py

The output will be stored in the output folder specified. It will either be stored in a GCF subdir (find_top_hits.py) or a GCF subdir (find_top_hits_GCA.py). 

the subdir will look like this: 

```bash 
└── 📁GCA
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

The databases folder contains databases made by diamon to blastp against for all isolates. 

The top_five folder contains the top five hits for each target protein

The top_hits folder contains the top hit for each protein if the bitscore threshhold was surpassed. if not, then it contains an "ambiguous" folder with the top two hits for that protein. It does not contain any empty files. 

The top_hits_consistent folder contains a winners and ambiguous file for each isolate, regardless of if the file is empty or not. Each file will only include one hit for each protein (ambiguous or not). 

If you chose to blastp against another set of target proteins, another folder (e.g. PAO1_target_genes) will be generated with the same subdirs under the GCA folder. 

If you use find_top_hits.py the main folder will be named GCF instead of GCA. 


## ⚖️ compare_hits_between_refs.py

If you have decided to run the isolates against two sets of target genes, then you can use the compare_hits_between_refs.py script to compare the top hits between either, and determine the top hit using the logic outlined in the repo diagram.  

### 🧪 Usage of compare_hits_between_refs.py

This script expects to be passed the same target genes folders as the previous scripts, one for each set of target genes blasted. Keep the order consisent (target_genes1 vs. target_genes2). 

This script also expects to be passed the folders generated by the previous scripts. Within the GCA (or GCF) folder, there will be a subdir for each set of target genes run. For example: 

```bash 
└── 📁GCA
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

The summary.txt file will have the following info: 

```bash
number of isolates: 39
number of unsure hits: 32
efflux genes: 10
pili genes: 6
porin genes: 1

Unsure: pili_RXTH01P.1.fsa_aa.gz.
		WKE25143.1 ABJ13792.1 RTS18461.1 RTS16263.1 WIN PAO1, WIN PA14
Unsure: pili_RXTU01P.1.fsa_aa.gz.
		WKE25143.1 ABJ13792.1 RTS87398.1 RTS87646.1 WIN PAO1, WIN PA14
Unsure: pili_RXUB01.1.fsa_aa.gz.
		WKE25143.1 ABJ13792.1 RTT37168.1 RTT34533.1 WIN PAO1, WIN PA14

```

# 🙋‍♀️ Author/ 📬 Contact

For questions or suggestions, contact: 

Hannah Kapoor
📧 hannahkapoor00@gmail.com 
