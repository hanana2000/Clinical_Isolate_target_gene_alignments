# 🧬 Pseudomonas Aeruginosa Target Gene Alignment 

This repository is for a project looking for and aligning known phage target genes (pili, efflux, and porins) from PAO1 and PA14 in 39 Pseudomonas clinical isolates. 


This diagram depicts the logic for determining top hits for each target protein in each clinical isolate: 

![gene search pipline diagram](https://github.com/hanana2000/Clinical_Isolate_target_gene_alignments/blob/501dad9a5ad1674b30afb49399bffdd9fe88ffe7/PA_Isolate_target_proteins.jpg)


## 🏆 top_hits_isolates folder 

- top_hits_isolates folder 
    - find_top_hits[_GCA].py: 
        - The find_top_hits.py file and find_top_hits_GCA.py file will both give the same output. They both find the top 5 hits, and from this the top hits (using logic outlined in diagram found in this repo based on bitscore), including top hit ambiguity. 
    - compare_hits_between_refs.py: 
        - If you have decided to run the isolates against two sets of target genes, then you can use the compare_hits_between_refs.py script to compare the top hits between either, and determine the top hit using the logic outlined in the repo diagram.  
- align_target_genes folder 
    - mafft_target_genes.py:
        - Aligns matched sets of target genes from two reference genomes and summarizes per-pair differences.

# 🙋‍♀️ Author/ 📬 Contact

For questions or suggestions, contact: 

Hannah Kapoor
📧 hannahkapoor00@gmail.com 
