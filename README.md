# Pectobacterium pangenome analysis

### IMP NCBI link
https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/


### Other notes

```bash
## generate heap dump
jmap -dump:live,format=b,file=heap_dump.bin 27422

## dump garbage collection logs
java -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails -Xloggc:<filename>

```

## Project workflow

```mermaid
%%{init: {"flowchart": {"htmlLabels": false}} }%%
flowchart TB
subgraph data["`*Pectobacterium* spp. genomes`"]
    direction LR
    inhouse[("Inhouse: 56")]
    ncbi[("NCBI: 450")]
    inhouse ~~~ ncbi
end
data --> Preprocessing
subgraph Preprocessing
    direction LR
    prokka[Prokka annotation]
    interproscan["InterProScan & COG annotation"]
    prokka --- interproscan
end
Preprocessing --> QC
subgraph "QC"
    direction LR
    ncbi_qc["fa:fa-ban Remove genomes with NCBI QC tags:\nAnomalous/Replaced/Excluded"]
    busco[BUSCO >= 99%]
    ani["ANI based\ntaxonomy correction"]
    duplicate["Remove duplicate \nassemblies"]
    ncbi_qc ~~~ busco
    ani ~~~ duplicate
end
subgraph pangenome[Pangenome with 454 genomes]
    direction LR
    construct[Build pangenome]
    group[Group orthologs]
    analyze[Association analysis]
    viz[Visualize & summarize]
    construct ~~~ group
    analyze ~~~ viz
end
QC --> pangenome
```

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

#### ANI data analysis

```mermaid
flowchart TD
    ANI_tree[c01_pangenome_ANI_tree.R]
    ANI_distribution[c01_intraspecies_ANI_distribution.R]
    plot_trees[c02_plot_evolutionary_trees.R]

    ANI_tree --> ANI_distribution
    ANI_tree --> plot_trees
```

#### *In*-sliico assay validation

```mermaid
flowchart TD
```

#### Homology group summary

```mermaid
flowchart TD
```

#### Phylogeny processing

```bash
Rscript scripts/c_analysis/phylo.process_newicks.R
Rscript scripts/c_analysis/phylo.correlate_trees.R
```

#### Phenotype association

```bash
Rscript scripts/c_analysis/phylo.explore_clades.R
Rscript scripts/c_analysis/phylo.make_clade_cmp_conf.R 
```

```mermaid
flowchart TD
```

#### Pangenome characterization for species group

```mermaid
flowchart TD
```

#### Temp

```bash
nohup bash scripts/b_construction/b03_build_pangenome.v2.sh pectobacterium.v2 1 >>logs/pantools/pheno_association.log 2>&1 &
```
