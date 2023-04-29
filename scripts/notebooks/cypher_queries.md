# Cypher queries for pangenome exploration

## Explore basic feature links

### Gene, mRNA and CDS relationships

```cypher
MATCH (m:mRNA) RETURN m.id, m.genome, m.protein_ID LIMIT 25
```

mRNA to gene relationship

```cypher
## mRNA<-gene
MATCH (g:gene)-[:codes_for]->(m:mRNA) 
RETURN g.genome AS genome, g.sequence AS chr_num, g.begin AS start,
g.end AS end, g.strand AS strand, g.id AS gene_name, m.id as mRNA_id
LIMIT 20
```

```cypher
## CDS->mRNA<-gene
MATCH (c:CDS)-[:contributes_to]->(m:mRNA)<-[:codes_for]-(g:gene) RETURN c, m, g LIMIT 20
```

> **Note**
> This running for very long time 
> Need to debug 

```cypher
## 
MATCH (g:gene)-[:codes_for]->(m:mRNA) 
WITH * 
MATCH (chr:sequence) 
WHERE chr.genome = g.genome AND chr.sequence = g.sequence
WITH chr, g, m
RETURN g.genome, g.sequence, g.begin, g.end, g.strand, g.id, m.id, chr.title
LIMIT 20

```

### Chromosome information

```cypher
MATCH (chr:sequence)
RETURN chr.genome AS genome, chr.number AS chrNum, chr.identifier AS chrId,
chr.title AS chrName
LIMIT 20
```

### Homology groups

```cypher
MATCH (m:mRNA)<-[:has_homolog]-(h:homology_group)
WHERE h.group_version = 1
RETURN m.id AS mRNA_id, m.genome AS genome, m.sequence AS chr,
id(h) AS hg_id, h.group_version AS hg_ver
LIMIT 25
```

### mRNA - homology_group - GO link data

```cypher
## explore structure
MATCH (m:mRNA)<-[:has_homolog]-(hg:homology_group) 
WHERE hg.group_version = 1
OPTIONAL MATCH (m)-[:has_go]->(go:GO)
OPTIONAL MATCH (m)<-[:is_parent_of]-(f:feature)
RETURN m, go, hg
LIMIT 25


MATCH (m:mRNA)<-[:has_homolog]-(hg:homology_group) 
WHERE hg.group_version = 1
OPTIONAL MATCH (m)-[:has_go]->(g:GO) 
RETURN m.id AS mRNA_id, m.genome AS genome, m.sequence AS chr, 
g.id AS go_id, id(hg) AS hg_id, hg.group_version AS hg_ver,
m.COG_id AS COG_id, m.COG_description AS COG_description,
m.COG_category AS COG_category
LIMIT 25
```

### mRNA - COG

```cypher
MATCH (m:mRNA)
RETURN m.id AS mRNA_id, m.genome AS genome, m.sequence AS chr, 
m.COG_id AS COG_id, m.COG_description AS COG_description,
m.COG_category AS COG_category
LIMIT 25 
```


## Test 

```cypher

MATCH (g:GO) RETURN g.id LIMIT 25

MATCH (m:mRNA)-->(g:GO) WHERE (m)-[:has_go]->(g) RETURN m.id, m.genome, m.protein_ID, g.id LIMIT 25

MATCH (m:mRNA)-[r:has_go]->(g:GO) WITH * RETURN m.id, m.genome, m.protein_ID, g.id LIMIT 25

MATCH p = (m:mRNA)-[r:has_go]->(g:GO) RETURN p LIMIT 5




MATCH (m:mRNA {id: 'KJEENCHO_04179'})-->(g:GO) RETURN g

MATCH (n:GO) WHERE any(f IN n.frequency WHERE f = 1) RETURN n.id, size(n.frequency), n.frequency LIMIT 20

MATCH (m:mRNA)<-[:has_homolog]-(h:homology_group) WHERE h.group_version = 1
OPTIONAL MATCH (m)-[:has_go]->(g:GO)
RETURN m.id AS mRNA_id, m.genome AS genome, m.sequence AS chr, g.id AS go_id, id(h) AS hg_id, h.group_version AS hg_ver

```