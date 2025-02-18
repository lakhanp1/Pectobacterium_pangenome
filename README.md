<p align="center">

![GitHub Release](https://img.shields.io/github/v/release/lakhanp1/pectobacterium_pangenome)
![GitHub Release Date](https://img.shields.io/github/release-date/lakhanp1/Pectobacterium_pangenome)
![GitHub commits since latest release](https://img.shields.io/github/commits-since/lakhanp1/Pectobacterium_pangenome/latest)
![GitHub License](https://img.shields.io/github/license/lakhanp1/pectobacterium_pangenome)
<br>
[![Static Badge](https://img.shields.io/badge/doi-10.5281%2Fzenodo.14882427-blue)](https://doi.org/10.5281/zenodo.12772014)
<p>

---

# Pectobacterium pangenome analysis

Please refer to https://lakhanp1.github.io/Pectobacterium_pangenome for data
preprocessing, pangenome construction and prophage analysis workflow documents.

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

## [Preprocessing](https://lakhanp1.github.io/Pectobacterium_pangenome/scripts/notebooks/preprocessing.html)

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

## [Pangenome construction](https://lakhanp1.github.io/Pectobacterium_pangenome/scripts/notebooks/pangenome_construction.html)

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

## [Prophage analysis](https://lakhanp1.github.io/Pectobacterium_pangenome/scripts/notebooks/prophage_analysis.html)

Tracing prophages in *Pectobacterium* genus pangenome to study their dynamics and contribution to the pangenome growth.