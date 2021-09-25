# DeBekker_and_Das_2021

REVIEW PAPER

*Parasite Immunology*

### Datasets required

1. Time-course daily transcriptomes for your host-of-interest 
1.1 Host: *Camponotus floridanus* or Florida carpenter ants
1.2 Datasets: normalized gene expression for forager and nurse ant brains
1.3 Resoultion: collected every 2h, over a 24h period

2. Genes underlying caste differentiation (behavioral plasticity) in ants
2.1 Genes sig. differentially expressed between foragers and nurse ant brains

3. Gene underlying parasite-induced behavioral manipiulation in ants
3.1 Host: *Camponotus floridanus*
3.2 Parasite: *Ophiocordyceps camponoti-floridani*
3.3 Genes sig. differentially expressed in ant heads during parasite-induced manipulated biting behavior as compared to uninfected controls.

### Goals:

Step 1
- **Build the circadian gene co-expression network (GCN)** 
    - clean up the data
    - use only genes that show “expression” throughout the day

Step 2
- **Identify clusters (modules) of highly co-expressed genes**

Step 3
- **Annotate the network**
    - identify where your genes-of-interest (GOI) are located in the network
        - Step 3.1.1: 01_GOI
            - genes underlying caste differentiation in ants
            - i.e., genes differentially expressed between foragers and nurses
        - Step 3.1.2: 02_GOI
            - genes underlying parasitic behavioral manipulation
            - i.e., genes differentially expressed in forager heads during manipulation (parasite induced biting behavior)

Step 4
- **Identify the functions of these genes using GO/PFAM enrichments**

Step 5
- **Form a data-driven hypothesis that can be empirically tested**

