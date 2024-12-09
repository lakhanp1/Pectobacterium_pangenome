
# Carotovoricin (CTV) cluster or phage_grp_1 data analysis

## Extract all homology groups representing CTV region across the pangenome

Genes and their respective homology groups flanking CTV:
- ybiB: hg_22427603
- tolC_2: hg_22427641

Extract genomic coordinates for CTV region together with the haplotypes

```bash
Rscript scripts/utils/HG_range_coordinates.R --hgs hg_22427603,hg_22427641 \
--inner_region --haplotypes \
--out analysis/pangenome_v2/carotovoricin/ctv_region/hg_regions.tab
```

Process the output file from the previous step to build CTV homology group
table and then manually annotate the homology groups with broad functionl
categories based on the keywords in GO, COG and PFAM annotations.

Intact CTV is present in 366 genomes.

CTV is missing or deleted in the following 25 genomes:
- *P. brasiliense*: g_149, g_177, g_182, g_185, g_236
- *P. carotovorum*: g_15
- *P_atrosepticum*: g_124, g_200, g_27, g_283, g_311, g_322, g_336, g_348, g_382,
  g_385, g_389, g_392, g_396, g_50, g_51, g_55
- *P. betavasculorum*: g_383, g_386
- *P. cacticida*: g_451
- *P. versatile*: g_313 (deletion of 5 genes including tail sheath, tail tube
  and tape measure genes). Check g_313.vir_2

Genome g_149 CTV cluster has integration of another prophage and hence it is 
disrupted. It was detected as g_149.vir_1	and was part of phage_grp_89 with 3 members.

CTV fragmented into multiple contigs in rest of the genomes.

## Clustermap visualization of CTV clusters

### CTV across pangenome

Include *P. atrosepticum* (g_385), *P. betavasculorum* (g_386) and *P. cacticida*
(g_451) genomes that lack CTV cluster.

```bash
Rscript scripts/utils/HG_range_coordinates.R --hgs hg_22427603,hg_22427641 \
--inner_region --genomes "g_385,g_386,g_451"
```

```r
cluster_title <- "ctv_typeStrains"
outDir <- paste(confs$analysis$ctv$dir, "/cluster_viz", sep = "")
hg_color_categories <- confs$analysis$ctv$files$hg_broad_functions

# a vector of prophage identifiers that will be included in clustermap plot
region_cluster <- NA
other_regions <- c(
  "g_345.vir_1", "g_446.vir_4", "g_66.vir_3", "g_222.vir_2", "g_365.vir_3",
  "g_442.vir_1", "g_8.vir_2", "g_38.vir_2", "g_273.vir_2", "g_259.vir_4",
  "g_305.vir_1", "g_378.vir_6", "g_428.vir_1", "g_248.vir_1", "g_449.vir_1",
  "g_54.vir_1", "g_116.vir_3", "g_423.vir_3", "g_375.vir_2", "g_381.vir_2"
)

subSample <- FALSE
cutHeight <- 1.5
addFlankingRegions <- TRUE
flankingRegion <- 5000

# ordering factor for prophages: "host" phylogeny, "hg_pav" for prophage HG PAV,
# "cluster_mash" for prophage MASH and "default" to use the provided
clusterOrder <- "host" # host, hg_pav, cluster_mash, default

# whether to keep custom regions at the bottom or consider during phylogeny
# based ordering
regions_phy_ordered <- TRUE

# regions to append as list of list with following structure
# list(r1 = list(chr, start, end, genomeId), r2 = list(chr, start, end, genomeId))
customRegions <- list(
  g_385_reg = list(chr = "NZ_JQHK01000003.1", start = 203963, end = 207120, genomeId = "g_385"),
  g_386_reg = list(chr = "NZ_JQHM01000001.1", start = 553213, end = 555615, genomeId = "g_386"),
  g_451_reg = list(chr = "Contig_2_668.636", start = 191452, end = 191490, genomeId = "g_451")
)

```

### carotovoricin cluster absent in some *P. brasiliense* isolates

*P. brasiliense* lacking CTV cluster: g_149, g_177, g_182, g_185, g_236. Genome
g_149 CTV cluster has integration of another prophage and hence it is disrupted.

```bash
Rscript scripts/utils/HG_range_coordinates.R --hgs hg_22427641,hg_22427603 \
--inner_region --genomes "g_149,g_177,g_182,g_185,g_236,g_15,g_313"
```

```r
cluster_title <- "ctv_pbr"
outDir <- paste(confs$analysis$ctv$dir, "/cluster_viz", sep = "")
hg_color_categories <- confs$analysis$ctv$files$hg_broad_functions

# a vector of prophage identifiers that will be included in clustermap plot
region_cluster <- "phage_grp_1"
other_regions <- c("g_368.vir_3")

subSample <- TRUE
cutHeight <- 1.5
addFlankingRegions <- TRUE
flankingRegion <- 5000

# ordering factor for prophages: "host" phylogeny, "hg_pav" for prophage HG PAV,
# "cluster_mash" for prophage MASH and "default" to use the provided
clusterOrder <- "host" # host, hg_pav, cluster_mash, default

# whether to keep custom regions at the bottom or consider during phylogeny
# based ordering
regions_phy_ordered <- FALSE

# regions to append as list of list with following structure
# list(r1 = list(chr, start, end, genomeId), r2 = list(chr, start, end, genomeId))
customRegions <- list()

regionClusters <- suppressMessages(
  readr::read_tsv(confs$analysis$prophages$files$clusters)
)

# optional filter
regionClusters %<>%
  dplyr::filter(
    SpeciesName == "P. brasiliense",
    !prophage_id %in% c("g_408.vir_3", "g_403.vir_3", "g_399.vir_3")
  )

```

### Visualize CTV deletion

Extract genomic coordinates for the species/genomes in which CTV is deleted.

```bash
Rscript scripts/utils/HG_range_coordinates.R --hgs hg_22427603,hg_22427641 \
--inner_region --genomes "g_149,g_177,g_182,g_185,g_236,g_15,g_313,g_385,g_386"
```

```r
cluster_title <- "ctv_deletion"
outDir <- paste(confs$analysis$ctv$dir, "/cluster_viz", sep = "")
hg_color_categories <- confs$analysis$ctv$files$hg_broad_functions

# a vector of prophage identifiers that will be included in clustermap plot
region_cluster <- NA
other_regions <- c("g_345.vir_1", "g_66.vir_3", "g_8.vir_2", "g_381.vir_2")

subSample <- FALSE
cutHeight <- 1.5
addFlankingRegions <- TRUE
flankingRegion <- 5000

# ordering factor for prophages: "host" phylogeny, "hg_pav" for prophage HG PAV,
# "cluster_mash" for prophage MASH and "default" to use the provided
clusterOrder <- "host" # host, hg_pav, cluster_mash, default

# whether to keep custom regions at the bottom or consider during phylogeny
# based ordering
regions_phy_ordered <- TRUE

# regions to append as list of list with following structure
# list(r1 = list(chr, start, end, genomeId), r2 = list(chr, start, end, genomeId))
customRegions <- list(
  g_177_reg = list(chr = "NZ_JACGEP010000002.1", start = 102352, end = 103187, genomeId = "g_177"),
  g_182_reg = list(chr = "NZ_JACGZZ010000050.1", start = 81626, end = 82461, genomeId = "g_182"),
  g_185_reg = list(chr = "NZ_JACGEN010000006.1", start = 81658, end = 82493, genomeId = "g_185"),
  g_236_reg = list(chr = "NZ_JACDSF010000027.1", start = 89639, end = 90475, genomeId = "g_236"),
  g_149_reg = list(chr = "NZ_JACGFD010000001.1", start = 2861612, end = 2913419, genomeId = "g_149"),
  g_385_reg = list(chr = "NZ_JQHK01000003.1", start = 203963, end = 207120, genomeId = "g_385"),
  g_386_reg = list(chr = "NZ_JQHM01000001.1", start = 553213, end = 555615, genomeId = "g_386"),
  g_451_reg = list(chr = "Contig_2_668.636", start = 191452, end = 191490, genomeId = "g_451"),
  g_313_reg = list(chr = "NZ_CP024842.1", start = 2917505, end = 2929851, genomeId = "g_313"),
  g_15_reg = list(chr = "BRCR01000012.1", start = 138656, end = 139632, genomeId = "g_15")
)

```

Prophage clusters found in the ctv-lacking Pbr: phage_grp_30, phage_grp_6, phage_grp_29

### A prophage integrated within the CTV in *P. brasiliense* genome `g_149`

<!--  -->
Genome g_149 CTV cluster has integration of another prophage and hence it is 
disrupted. It was detected as g_149.vir_1	and was part of phage_grp_89 which has
3 members, `r c("g_149.vir_1", "g_116.vir_1", "g_46.vir_2")`

```r
cluster_title <- "prophage_in_ctv"
outDir <- paste(confs$analysis$ctv$dir, "/cluster_viz", sep = "")
hg_color_categories <- confs$analysis$ctv$files$hg_broad_functions

# a vector of prophage identifiers that will be included in clustermap plot
region_cluster <- NA
other_regions <- c("g_368.vir_3", "g_149.vir_1", "g_116.vir_1", "g_46.vir_2")

subSample <- FALSE
cutHeight <- 1.5
addFlankingRegions <- TRUE
flankingRegion <- 5000

# ordering factor for prophages: "host" phylogeny, "hg_pav" for prophage HG PAV,
# "cluster_mash" for prophage MASH and "default" to use the provided
clusterOrder <- "default" # host, hg_pav, cluster_mash, default

# whether to keep custom regions at the bottom or consider during phylogeny
# based ordering
regions_phy_ordered <- FALSE

# regions to append as list of list with following structure
# list(r1 = list(chr, start, end, genomeId), r2 = list(chr, start, end, genomeId))
customRegions <- list()

```

### CTV in *P. versatile*

```r
grpToView <- "ctv_pvs"
subSample <- TRUE
cutHeight <- 1.5
addFlankingRegions <- TRUE
flankingRegion <- 5000


grp <- list(
  phage_grp = grpToView,
  members = dplyr::filter(
    regionClusters,
    SpeciesName == "P. versatile", nFragments == 1, phage_grp == "phage_grp_1"
  ) %>%
    dplyr::pull(prophage_id)
)
```

### CTV cluster in *P. versatile* collected from France

```r
grpToView <- "ctv_pvs_fr"
subSample <- TRUE
cutHeight <- 1.5
addFlankingRegions <- TRUE
flankingRegion <- 5000

grp <- list(
  phage_grp = grpToView,
  members = dplyr::filter(
    regionClusters,
    SpeciesName == "P. versatile", nFragments == 1, phage_grp == "phage_grp_1",
    geo_loc_country == "France"
  ) %>%
    dplyr::pull(prophage_id)
)
```

## Carotovoricin tail fiber variation region MSA

Extract the region between the two homology groups, `hg_22427604` and `hg_22427603`
for the carotovoricin cluster. Additionally, generate `GFF3` files with the
homology group, COG, PFAM and other metadata information to visualize.

```bash
# Rscript scripts/analysis/clustermap_get_regions.R
Rscript scripts/utils/HG_range_coordinates.R --hgs hg_22427604,hg_22427603 \
--genomes "g_302,g_364,g_337,g_439,g_403,g_399,g_408,g_175,g_138,g_368,g_345,\
g_366,g_308,g_191,g_173,g_155,g_166,g_299,g_438,g_263,g_43,g_391" \
--out analysis/pangenome_v2/prophages/cluster_viz/ctv_pbr/ctv_pbr.variable_regions.tab

# create CTV regions GFF3 files for one species representative each
Rscript scripts/utils/HGs_gff_from_regions.R \
--regions analysis/pangenome_v2/carotovoricin/tail_region/hg_regions.tab
```

Extract DNA sequence for these 3' tail fiber regions.

```bash
mkdir analysis/pangenome_v2/prophages/cluster_viz/ctv_pbr/region_fasta

tail -n +2 analysis/pangenome_v2/prophages/cluster_viz/ctv_pbr/ctv_pbr.variable_regions.tab |
  cut -f 1,3-5 |
  while read name sampleId pos strand; do
    rc=""
    if [ ${strand} == "-" ]; then
        rc="--reverse-complement"
    fi
    outFa="analysis/pangenome_v2/prophages/cluster_viz/ctv_pbr/region_fasta/${name}.fasta"

    seq=">${name} ${pos}(${strand})\n"
    seq+=`samtools faidx ${rc} data/prokka_annotation/${sampleId}/${sampleId}.fna $pos | tail -n +2`
    printf "${seq}\n" > ${outFa}
    printf "${seq}\n"
  done >analysis/pangenome_v2/prophages/cluster_viz/ctv_pbr/ctv_pbr.variable_regions.fasta

```

Mauve aligns sequences provided in `fasta` and `genbank` format. However, the
sequence annotation is shown only if the alignment input was in `genbank`
format. Therefore, we need to combine the homology group annotation in `gff3`
format with the `fasta` sequences to generate a `genbank` formatted files. To
do this, Emboss tool `seqret` was used.

```bash
conda activate omics_py37

cd analysis/pangenome_v2/prophages/cluster_viz/ctv_pbr/region_fasta

for vir in `cut -f 1 ../ctv_pbr.variable_regions.tab | tail -n +2`
do
  seqret -sequence ${vir}.fasta -feature -fformat gff3 \
  -fopenfile ${vir}.gff3 -osformat genbank -osextension gbk \
  -osname_outseq ${vir} -ofdirectory_outseq gbk_file -auto 
done

sed -i.bak -r -e 's/(^\s+CDS\s+(complement\()?)[^:]+:([[:digit:]]+\.\.[[:digit:]]+\)?).*/\1\3/' \
  -e 's/\/note="\*([^:]+): /\/\1="/' *.gbk 

```

Align and visualize sequences using Mauve.

```bash
mauve_out="ctv_pbr"

/Applications/Mauve.app/Contents/MacOS/progressiveMauve \
--output="${mauve_out}.mauve.xmfa" \
--backbone-output="${mauve_out}.mauve.backbone" \
--output-guide-tree="${mauve_out}.mauve.guide_tree.newick" \
region_fasta/g_302.vir_3.gbk region_fasta/g_364.vir_2.gbk \
region_fasta/g_337.vir_1.gbk region_fasta/g_439.vir_2.gbk \
region_fasta/g_403.vir_3.gbk region_fasta/g_399.vir_3.gbk \
region_fasta/g_175.vir_3.gbk region_fasta/g_408.vir_3.gbk \
region_fasta/g_138.vir_2.gbk \
region_fasta/g_368.vir_3.gbk region_fasta/g_345.vir_1.gbk \
region_fasta/g_366.vir_3.gbk region_fasta/g_308.vir_2.gbk \
region_fasta/g_191.vir_1.gbk region_fasta/g_173.vir_2.gbk \
region_fasta/g_155.vir_1.gbk region_fasta/g_166.vir_2.gbk \
region_fasta/g_299.vir_1.gbk region_fasta/g_438.vir_2.gbk \
region_fasta/g_263.vir_2.gbk region_fasta/g_43.vir_1.gbk \
region_fasta/g_391.vir_3.gbk  > "${mauve_out}".log 2>&1

mauve_out="ctv_pbr.ein"

/Applications/Mauve.app/Contents/MacOS/progressiveMauve \
--output="${mauve_out}.mauve.xmfa" \
--backbone-output="${mauve_out}.mauve.backbone" \
--output-guide-tree="${mauve_out}.mauve.guide_tree.newick" \
region_fasta/g_345.vir_1.gbk region_fasta/g_366.vir_3.gbk \
region_fasta/g_308.vir_2.gbk region_fasta/g_191.vir_1.gbk \
region_fasta/g_173.vir_2.gbk region_fasta/g_155.vir_1.gbk \
region_fasta/g_166.vir_2.gbk  > "${mauve_out}".log 2>&1


mauve_out="ctv_pbr.rep"

/Applications/Mauve.app/Contents/MacOS/progressiveMauve \
--output="${mauve_out}.mauve.xmfa" \
--backbone-output="${mauve_out}.mauve.backbone" \
--output-guide-tree="${mauve_out}.mauve.guide_tree.newick" \
region_fasta/g_302.vir_3.gbk region_fasta/g_337.vir_1.gbk \
region_fasta/g_175.vir_3.gbk region_fasta/g_173.vir_2.gbk \
region_fasta/g_299.vir_1.gbk region_fasta/g_438.vir_2.gbk \
region_fasta/g_391.vir_3.gbk  > "${mauve_out}".log 2>&1

```

Perform MSA using `MAFFT`

```bash
conda activate pantools_v4_3
cd analysis/pangenome_v2/prophages/cluster_viz/ctv_pbr/

mafft --globalpair --quiet --maxiterate 1000 --treeout ctv_pbr.variable_regions.fasta
```

Smash++ pairwise sequence comparison

```bash
conda activate omics_py37

cd analysis/pangenome_v2/prophages/cluster_viz/ctv_pbr/region_fasta

function smashpp_compare
{
  smashpp -r $1 -t $2
  smashpp viz -o $(basename $1 .fasta)$(basename $2 .fasta)".svg" "${1}.${2}.pos"
}

export -f smashpp_compare

smashpp_compare g_345.vir_1.fasta g_345.vir_1.fasta 
smashpp_compare g_345.vir_1.fasta g_366.vir_3.fasta 
smashpp_compare g_345.vir_1.fasta g_302.vir_3.fasta
smashpp_compare g_345.vir_1.fasta g_403.vir_3.fasta
smashpp_compare g_345.vir_1.fasta g_173.vir_2.fasta
smashpp_compare g_345.vir_1.fasta g_138.vir_2.fasta
smashpp_compare g_345.vir_1.fasta g_263.vir_2.fasta
smashpp_compare g_345.vir_1.fasta g_155.vir_1.fasta

```

Left inverted repeat for Ein: CTCCCGCAAACCTCGGTTTTGGGGAC (CTCCCGCAAACCTCGGTTT)  
Left inverted repeat for Ein(rev-com): GTCCCCAAAACCGAGGTTTGCGGGAG (AAACCGAGGTTTGCGGGAG, AAACCGAGGTTTGCG)  

Right inverted repeat for Ein: TTCTCGCAAACCTCGGTTTTGGAGAA  
Right inverted repeat for Ein(rev): AAGAGGTTTTGGCTCCAAACGCTCTT  
Right inverted repeat for Ein(rev-com): TTCTCCAAAACCGAGGTTTGCGAGAA (AAACCGAGGTTTGCGAGAA, AAACCGAGGTTTGCG)  
Right inverted repeat for Ein(comp): AAGAGCGTTTGGAGCCAAAACCTCTT  

## Inter-species horizontal gene transfer of CTV

Clustermap visualization for some manually selected candidates. Change the code
blocks in the script `scripts/analysis/clustersmap_data_prophages.R` and run.

```r
grpToView <- "ctv_hgt"
subSample <- FALSE 
cutHeight <- 1.5
addFlankingRegions <- TRUE
flankingRegion <- 5000

# ordering factor for prophages: host phylogeny, prophage HG PAV, prophage MASH,
# completeness score
clusterOrder <- "hg_pav" # completeness, host, hg_pav, cluster_mash


# optionally, a custom region list can be provided to generate the plot
grp <- list(
  phage_grp = grpToView,
  members = c(
    "g_145.vir_1", "g_194.vir_1", "g_429.vir_1", "g_442.vir_1", "g_421.vir_1", 
    "g_150.vir_4", "g_447.vir_1", "g_434.vir_4",
    "g_221.vir_3", "g_53.vir_3", "g_106.vir_2", "g_57.vir_1", "g_125.vir_1"
  )
)
```

Run script `scripts/utils/HG_range_coordinates.R` to extract genomic coordinates
for CTV region in selected genomes.

Complete CTV region along with 3 flanking genes:

```bash
Rscript scripts/utils/HG_range_coordinates.R --hgs hg_22427643,hg_22427599 \
--inner_region --out analysis/pangenome_v2/carotovoricin/ctv_region/hg_regions.tab 

# create GFF3 files for CTV regions of interest
Rscript scripts/utils/HGs_gff_from_regions.R \
--regions analysis/pangenome_v2/carotovoricin/ctv_region/hg_regions.tab \
--genomes g_145,g_194,g_429,g_442,g_421,g_150,g_447,g_434,g_53,g_106,g_57,g_221,g_125

# create CTV regions GFF3 files for one species representative each
Rscript scripts/utils/HGs_gff_from_regions.R \
--regions analysis/pangenome_v2/carotovoricin/ctv_region/hg_regions.tab \
--genomes g_183,g_423,g_277,g_375,g_116,g_449,g_446,g_378,g_296,g_428,g_259,g_66,g_337,g_442,g_305,g_381,g_8,g_273,g_38
```

CTV conserved loci:

```bash
Rscript scripts/utils/HG_range_coordinates.R --hgs hg_22427640,hg_22427604 \
--out analysis/pangenome_v2/carotovoricin/ctv_conserved/hg_regions.tab
```

CTV variable loci:

```bash
Rscript scripts/utils/HG_range_coordinates.R --hgs hg_22427604,hg_22427603 \
--inner_region --haplotypes \
--out analysis/pangenome_v2/carotovoricin/ctv_tail/hg_regions.tab 
```

CTV tape measure protein gene:

```bash
Rscript scripts/utils/HG_range_coordinates.R --hgs hg_22427622,hg_22427616 \
--haplotypes --inner_region --overlapping \
--out analysis/pangenome_v2/carotovoricin/tape_measure/hg_regions.tab
```

Run `scripts/analysis/ctv_hgt.qmd` script to perform HGT analysis.

### DNA sequence comparison

#### Extract DNA sequence for the regions

```bash
# file_regions="analysis/pangenome_v2/carotovoricin/ctv_tail/hg_regions.tab"
# file_regions="analysis/pangenome_v2/carotovoricin/ctv_conserved/hg_regions.tab"
# file_regions="analysis/pangenome_v2/carotovoricin/ctv_tail/selected_haplotypes.tab"
# file_regions="analysis/pangenome_v2/carotovoricin/ctv_tail/selected_haplotypes_flanking.tab"
# file_regions="analysis/pangenome_v2/carotovoricin/ctv_conserved/selected_conserved.tab"
# file_regions="analysis/pangenome_v2/carotovoricin/upstream_core/selected_upstream.tab"

dir_path=$(dirname "${file_regions}")
out_fasta="${file_regions%.*}".fasta

tail -n +2 ${file_regions} |
  cut -f 1,3-5 |
  while read name sampleId pos strand; do
    rc=""
    if [ ${strand} == "-" ]; then
        rc="--reverse-complement"
    fi
    outFa="${dir_path}/fasta/${name}.fasta"

    seq=">${name}\n"
    seq+=`samtools faidx ${rc} data/prokka_annotation/${sampleId}/${sampleId}.fna $pos | tail -n +2`

    # printf "${seq}\n" > ${outFa}
    printf "${seq}\n"
  done > ${out_fasta}

```

#### Calculate the distance using `mash`

```bash
cd analysis/pangenome_v2/carotovoricin/ctv_tail
mash dist -p 8 -k 12 -s 2000 -i -S 124 hg_regions.fasta hg_regions.fasta > mash_dist.tab

cd analysis/pangenome_v2/carotovoricin/ctv_conserved
mash dist -p 8 -k 12 -s 2000 -i -S 124 hg_regions.fasta hg_regions.fasta > mash_dist.tab

cd analysis/pangenome_v2/carotovoricin/tape_measure
mash dist -p 8 -k 12 -s 2000 -i -S 124 hg_regions.fasta hg_regions.fasta > mash_dist.tab

grep -e 'g_(53|106|57|221|125).*g_(53|106|57|221|125)' ctv_*/ctv_dist.tab
grep -e 'g_(150|447|434).*g_(150|447|434)' ctv_*/ctv_dist.tab
grep -e 'g_(53|106|57).*g_(53|106|57)' ctv_*/ctv_dist.tab

```

#### Calculate Jaccard distance matrix using `sourmash`

MASH distance is low for haplotype pairs where one is subset of the other, eg,
A-B-C-D and A-B. Therefore, another tool was tried to calculate the Jaccard
index between sequences.

```bash
# file_fa="ctv_conserved/hg_regions.fasta"
# file_fa="ctv_tail/hg_regions.fasta"

file_sig="${file_fa%.*}".sig.zip
file_dist="${file_fa%.*}".sourmash.csv

# # sketch: scaled=5
# sourmash sketch dna -f -p k=7,k=9,k=11,scaled=5,abund,seed=124 --singleton \
# -o ${file_sig} ${file_fa}

# sketch: num=2000
sourmash sketch dna -f -p k=7,k=9,k=11,k=13,num=2000,abund,seed=124 --singleton \
-o ${file_sig} ${file_fa}

sourmash compare --ignore-abundance --distance-matrix -p 12 -k 13 --dna \
--csv ${file_dist} ${file_sig} 

```

#### Calculate Jaccard distance using `dashing`

```bash
# dir_fa="ctv_conserved/fasta"
# dir_fa="ctv_tail/fasta"

file_dist="${dir_fa%/*}/"distance.dashing.txt

dashing cmp -k 11 --seed 124 --nthreads 12 --full-mash-dist  \
-O ${file_dist} ${dir_fa}/*.fasta

```

### Tail fiber haplotype examples

Extract sequences for some example tail fiber loci haplotypes that are present
in multiple *Pectobacterium* species.

```bash
# tail fiber locus
Rscript scripts/utils/HG_range_coordinates.R --hgs hg_22427604,hg_22427603 \
--genomes "g_279,g_425,g_149,g_377,g_100,g_106,g_331,g_249,g_125,\
g_221,g_53,g_395,g_108,g_444,g_160" \
--out analysis/pangenome_v2/carotovoricin/ctv_tail/selected_haplotypes_flanking.tab

Rscript scripts/utils/HG_range_coordinates.R --hgs hg_22427604,hg_22427603 \
--genomes "g_279,g_425,g_149,g_377,g_100,g_106,g_331,g_249,g_125,\
g_221,g_53,g_395,g_108,g_444,g_160" \
--inner_region --out analysis/pangenome_v2/carotovoricin/ctv_tail/selected_haplotypes.tab

# conserved CTV locus
Rscript scripts/utils/HG_range_coordinates.R --hgs hg_22427640,hg_22427604 \
--genomes "g_279,g_425,g_149,g_377,g_100,g_106,g_331,g_249,g_125,\
g_221,g_53,g_395,g_108,g_444,g_160" \
--out analysis/pangenome_v2/carotovoricin/ctv_conserved/selected_conserved.tab

# upstream core
Rscript scripts/utils/HG_range_coordinates.R --hgs hg_22427642,hg_22427641 \
--genomes "g_279,g_425,g_149,g_377,g_100,g_106,g_331,g_249,g_125,\
g_221,g_53,g_395,g_108,g_444,g_160" \
--out analysis/pangenome_v2/carotovoricin/upstream_core/selected_upstream.tab
```

[DNA sequence comparison of these regions](#dna-sequence-comparison)

#### Perform MSA using `MAFFT`

```bash
nohup mafft --reorder --allowshift --unalignlevel 0.2 --leavegappyregion \
--maxiterate 0 --globalpair ctv_tail/selected_haplotypes.fasta \
> ctv_tail/selected_haplotypes.msa.fasta 2>nohup.out &

nohup mafft --reorder --allowshift --unalignlevel 0.2 --leavegappyregion \
--maxiterate 0 --globalpair ctv_conserved/selected_conserved.fasta \
> ctv_conserved/selected_conserved.msa.fasta 2>nohup.out &

nohup mafft --reorder --allowshift --unalignlevel 0.2 --leavegappyregion \
--maxiterate 0 --globalpair upstream_core/selected_upstream.fasta \
> upstream_core/selected_upstream.msa.fasta 2>nohup.out &
```

Generate a maximum-likelihood phylogenetic tree for MSAs

```bash
cd cd analysis/pangenome_v2/carotovoricin/

nohup nice iqtree2 -T 40 -s ctv_tail/selected_haplotypes.msa.fasta -B 1000 \
--prefix ctv_tail/iqtree/selected_haplotypes.msa.fasta >> iqtree_tree.log 2>&1 &

nohup nice iqtree2 -T 40 -s ctv_conserved/selected_conserved.msa.fasta -B 1000 \
--prefix ctv_conserved/iqtree/selected_conserved.msa.fasta >> iqtree_tree.log 2>&1 &

nohup nice iqtree2 -T 40 -s upstream_core/selected_upstream.msa.fasta -B 1000 \
--prefix upstream_core/iqtree/selected_upstream.msa.fasta >> iqtree_tree.log 2>&1 &

```



### Use Mauve to detect the structural variation and consevation

Mauve aligns sequences provided in `fasta` and `genbank` format. However, the
sequence annotation is shown only if the alignment input was in `genbank`
format. Therefore, we need to combine the homology group annotation in `gff3`
format with the `fasta` sequences to generate a `genbank` formatted files. To
do this, Emboss tool `seqret` was used.

```bash
conda activate omics_py37

for vir in `cut -f 1 hg_regions.tab | tail -n +2`
do
  seqret -sequence ${vir}.fasta -feature -fformat gff3 \
  -fopenfile ${vir}.gff3 -osformat genbank -osextension gbk \
  -osname_outseq ${vir} -ofdirectory_outseq gbk_file -auto 
done
```

Correct GenBank file as per Mauve requirement

```bash
sed -i.bak -r -e 's/(^\s+CDS\s+(complement\()?)[^:]+:([[:digit:]]+\.\.[[:digit:]]+\)?).*/\1\3/' \
  -e 's/\/note="\*([^:]+): /\/\1="/' *.gbk
```

Align and visualize sequences using Mauve.

```bash
mauve_out="ctv_hgt.defaults"

progressiveMauve --output="${mauve_out}.mauve.xmfa" \
--backbone-output="${mauve_out}.mauve.backbone" \
--output-guide-tree="${mauve_out}.mauve.guide_tree.newick" \
*.gbk  > "${mauve_out}".mauve.log 2>&1

# mauve with custom settings
mauve_out="ctv_hgt.opt"

progressiveMauve --output="${mauve_out}.mauve.xmfa" \
--seed-weight 10 \
--backbone-output="${mauve_out}.mauve.backbone" \
--output-guide-tree="${mauve_out}.mauve.guide_tree.newick" \
*.gbk  > "${mauve_out}".mauve.log 2>&1

# --seed-weight 5 \

```

### Visualize tail fiber locus HG sets on the pangenome

R script `scripts/analysis/HG_tandem_match.viz.qmd` to plot tail fiber locus HGs
for a set of genomes from different species.

```r
hgSets <- list(
  s1 = c("hg_22427604", "hg_22426828"),
  s2 = c("hg_22427604", "hg_22426818"),
  s3 = c("hg_22427604", "hg_22426820", "hg_22426802", "hg_22426815")
)

outDir <- paste(confs$analysis$prophages$dir, "/cluster_viz/ctv_hgt", sep = "")
```

## O-antigen coevolution

### O-antigen homology groups

#### ABC transporter dependent O-antigen synthesis

- *Wzm* (ABC transporter permease):
- *Wzt* (ABC transporter ATP binding protein):

#### G4C or O-antigen capsule cluster

- *Wza*: hg_22428205
- *Wzb*: hg_22428204
- *Wzc*: hg_22428201
- *WaaL* (ligase):
- *rfbB*: found & conserved, GCF_009873295.1
- *rfbA*: found & conserved
- *rfbC*: found & conserved

- ?: hg_22428182
- *GfcB*: hg_22428181
- *GfcC*: hg_22428167
- *YbjH*: hg_22428166

- *wecA* (hg_22429666)

#### *Wzx/Wzy*-dependent pathway

- *Wzx* (flipase):
- *Wxy* (O-antigen polymerase):
- *Wzz* (chain length regulator):

*Dickeya solani* genome : https://www.ncbi.nlm.nih.gov/nuccore/CP015137.1

GCF_002068115.1: genes in LPS biosynthesis

### O-antigen HGs

hg_22428205,hg_22428166

```bash
Rscript scripts/utils/HG_range_coordinates.R --hgs hg_22428205,hg_22428166 \
--out analysis/pangenome_v2/o_antigen/abc_transporter_hg_regions.tab
```
