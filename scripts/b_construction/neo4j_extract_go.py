#!/usr/bin/env python3

from typing import List
from neo4j import GraphDatabase
import pandas as pd

# https://neo4j.com/docs/api/python-driver/current/#example-application
# https://neo4j.com/docs/api/python-driver/current/api.html


class QueryPangenome:
    def __init__(self, uri: str, user: str, password):
        self.driver = GraphDatabase.driver(uri, auth=(user, password))

    def close(self):
        self.driver.close()

    # get GO and homology group informaiton
    def get_GO(self, out: str, version: int) -> None:
        query = (
            "MATCH (m:mRNA)<-[:has_homolog]-(hg:homology_group) "
            f"WHERE hg.group_version = {version} "
            "OPTIONAL MATCH (m)-[:has_go]->(g:GO) "
            "RETURN m.id AS mRNA_id, m.genome AS genome, m.sequence AS chr, "
            "g.id AS go_id, id(hg) AS hg_id, hg.group_version AS hg_ver "
            # "LIMIT 25 "
        )

        col_names = ["mRNA_id", "genome", "chr", "go_id", "hg_id", "hg_ver"]

        # run transaction
        with self.driver.session() as session:
            session.execute_read(
                self._query_pangenome, query=query, cols=col_names, out=out
            )

    # get homology group information
    def get_homology_groups(self, version: int, out: str) -> None:
        query = (
            "MATCH (m:mRNA)<-[:has_homolog]-(h:homology_group) "
            f"WHERE h.group_version = {version} "
            "RETURN m.id AS mRNA_id, m.genome AS genome, m.sequence AS chr, "
            "id(h) AS hg_id, h.group_version AS hg_ver "
            # "LIMIT 25"
        )

        col_names = ["mRNA_id", "genome", "chr", "hg_id", "hg_ver"]

        # run transaction
        with self.driver.session() as session:
            session.execute_read(self._query_pangenome, query, cols=col_names, out=out)

    def get_chr_info(self, out: str) -> None:
        query = (
            "MATCH (chr:sequence) "
            "RETURN chr.genome AS genome, chr.number AS chr_num, "
            "chr.identifier AS chr_id, chr.title AS chr_name "
            # "LIMIT 20"
        )

        col_names = ["genome", "chr_num", "chr_id", "chr_name"]

        # run transaction
        with self.driver.session() as session:
            session.execute_read(self._query_pangenome, query, cols=col_names, out=out)

    def get_gene_info(self, out: str) -> None:
        query = (
            "MATCH (g:gene)-[:codes_for]->(m:mRNA) "
            "OPTIONAL MATCH (m)-[:has_pfam]->(p:pfam) "
            "RETURN g.genome AS genome, g.sequence AS chr_num, g.begin AS start, "
            "g.end AS end, g.strand AS strand, g.id AS gene_name, m.id as mRNA_id, "
            "m.COG_id AS COG_id, m.COG_description AS COG_description, "
            "m.COG_category AS COG_category "
            "p.id as pfam_id, p.description as pfam_description "
            # "LIMIT 20"
        )

        col_names = [
            "genome",
            "chr_num",
            "start",
            "end",
            "strand",
            "gene_name",
            "mRNA_id",
            "COG_id",
            "COG_description",
            "COG_category",
            "pfam_id",
            "pfam_description",
        ]

        # run transaction
        with self.driver.session() as session:
            session.execute_read(self._query_pangenome, query, cols=col_names, out=out)

    # Static methods, much like class methods, are methods that are bound to a
    # class rather than its object.
    @staticmethod
    def _query_pangenome(
        tx,
        query: str,
        cols: List[chr],
        out: str,
    ) -> None:
        # file to store the output
        fo = open(out, mode="w")

        # print("\t".join(cols))
        fo.write("\t".join(cols) + "\n")

        # df = pd.DataFrame(columns=cols)

        for record in tx.run(query):
            # vals = record.data(*cols)     # dict
            vals = record.values(*cols)  # list
            # print(type(record), type(vals), vals)
            fo.write("\t".join(map(str, vals)) + "\n")
            # print("\t".join(map(str, vals)))
            # df.loc[len(df)] = vals

        # print(df)
        if out is not None:
            fo.close()


if __name__ == "__main__":
    uri = "bolt://eddy.bioinformatics.nl:1241"
    user = "neo4j"
    password = "neo4j"

    panQ = QueryPangenome(uri, user, password)

    panQ.driver.verify_connectivity()

    ## get mRNA - GO - homology group data
    try:
        panQ.get_GO(out="pangenome_GO.tab", version=1)
    finally:
        panQ.close()

    ## get gene - mRNA - COG
    try:
        panQ.get_gene_info(out="pangenome_gene_info.tab")
    finally:
        panQ.close()

    ## get genome - chr# -chr_id data
    try:
        panQ.get_chr_info(out="pangenome_chr.tab")
    finally:
        panQ.close()
