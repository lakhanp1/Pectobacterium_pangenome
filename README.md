<div align="center">
    <p>
        <img alt="GitHub Release" src="https://img.shields.io/github/v/release/lakhanp1/pectobacterium_pangenome">
        <img alt="GitHub Release Date" src="https://img.shields.io/github/release-date/lakhanp1/pectobacterium_pangenome">
        <img alt="GitHub commits since latest release" src="https://img.shields.io/github/commits-since/lakhanp1/pectobacterium_pangenome/latest">
        <img alt="GitHub License" src="https://img.shields.io/github/license/lakhanp1/pectobacterium_pangenome">
    </p>
    <p>
        <img alt="GitHub pages" src="https://img.shields.io/badge/Project-GitHub%20pages-green?link=https%3A%2F%2Flakhanp1.github.io%2FPectobacterium_pangenome%2Fscripts%2Fnotebooks%2F">
        <img alt="bioRxiv" src="https://img.shields.io/badge/bioRxiv-10.1101%2F2024.09.02.610764-green?link=https%3A%2F%2Fdoi.org%2F10.1101%2F2024.09.02.610764">
        <img alt="Zenodo" src="https://img.shields.io/badge/Zenodo-10.5281%2Fzenodo.12772014-blue?link=https%3A%2F%2Fdoi.org%2F10.5281%2Fzenodo.12772014">
    </p>
</div>

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