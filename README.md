# Pectobacterium pangenome analysis

### IMP NCBI link
https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/

summarize ANI assembly report: visualize


### Other notes

```bash
## generate heap dump
jmap -dump:live,format=b,file=heap_dump.bin 27422

## dump garbage collection logs
java -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails -Xloggc:<filename>

```

## Workflow

### Preprocessing

```mermaid
flowchart TB
    eutils[NCBI Eutils]
    fasta[download FASTA]
    metadata[extract metadata]
    prokka[GFF]
    ani[ANI]
    interproscan[Functional annotation]
    subgraph raw_data [Raw data]
        direction TB
        subgraph d1[-]
            direction LR
            eutils --> fasta
            fasta --Prokka--> prokka
        end
        d1 ---->|InterProScan|interproscan
        fasta --> metadata
    end
    raw_data ------> QC[BUSCO and QUAST]
    subgraph data_cleanup[QC]
        direction TB
        fasta --fastANI--> ani
        ani --> taxonomy_check[Verify taxonomy]
        QC & metadata & taxonomy_check --> filtered[Filtered data]
    end
```

### Pangenome construction

```mermaid
flowchart TB
    filtered[Filtered data]
    subgraph pangenome[Pangenome construction]
        direction LR
        build_pangenome --> add_annotations
        add_annotations --> add_phenotypes
        add_phenotypes --> add_functions
    end
    
    filtered --> pangenome

    subgraph grouping[Homology grouping]
        direction LR
        optimal_grouping --> change_grouping
        change_grouping --> grouping_overview
    end

    pangenome --> grouping




```

### Analysis

```mermaid
flowchart TD
    ANI_tree[c01_pangenome_ANI_tree.R]
    ANI_distribution[c01_intraspecies_ANI_distribution.R]
    plot_trees[c02_plot_evolutionary_trees.R]

    ANI_tree --> ANI_distribution
    ANI_tree --> plot_trees
```
