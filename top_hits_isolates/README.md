# Finding top hits for each target using diamond Blastp

requires: 
- DIAMOND 
- biopython 
- pandas

## find_top_hits.py or find_top_hits_GCA.

First, either find_top_hits.py or find_top_hits_GCA.py must be run on the isolate/ target protein data. 

### Usage of find_top_hits.py or find_top_hits_GCA.py

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

### Output of find_top_hits.py and find_top_hits_GCA.py

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


