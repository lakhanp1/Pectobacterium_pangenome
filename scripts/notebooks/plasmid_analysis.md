# Plasmid analysis

## JGI IMG/PR data analysis

Search result for *Pectobacterium* genus plasmid query:

https://img.jgi.doe.gov/cgi-bin/plasmid/main.cgi?section=PlasmidBrowse&page=plasmiddetails&type=byhost&genus=g__Pectobacterium&domain=d__Bacteria&phylum=p__Proteobacteria&ir_class=c__Gammaproteobacteria&ir_order=o__Enterobacterales&family=f__Enterobacteriaceae&genus=g__Pectobacterium&col=Plasmids

## Compare plasmid identification in identical genomes

**Genomes with plasmid predicted:**

- g_401: Contigs NAK238_contig_15_consensus, NAK238_contig_7_consensus,
NAK238_contig_24_consensus, NAK238_contig_35_consensus are predicted as plasmids.

- g_13: BRCT01000006.1, BRCT01000014.1, BRCT01000028.1, BRCT01000026.1

Closest genomes with chromosomal assembly: g_214, g_188, g_187, g_190

Check if and where the plasmid homology groups map onto the chromosome level
genome assembly:

```bash
Rscript scripts/utils/HG_pav_genomes.R \
--genomes "g_401,g_13,g_214,g_188,g_187,g_190" --tree kmer_nj --node n25 \
--out analysis/pangenome_v2/plasmids/plasmid_evaluation_hgs.tab
```

Extract pairwise ANI for genomes of interest:

```bash
s1="NAK238|GCA_024343655.1_ASM2434365v1|GCF_003990515.2_ASM399051v2|GCF_016944595.1_ASM1694459v1"
s2="GCF_016949255.1_ASM1694925v1|NAK638|NAK644"
grep "($s1).fna.+($s2).fna" analysis/fastANI/ANI_results
```
