#!/usr/bin/env python3

from neo4j import GraphDatabase

# https://neo4j.com/docs/api/python-driver/current/#example-application


class QueryPangenome:
    def __init__(self, uri: str, user: str, password):
        self.driver = GraphDatabase.driver(uri, auth=(user, password))

    def close(self):
        self.driver.close()

    # get GO and homology group informaiton
    def get_GO(self, out: str, version: int) -> None:

        query = (
            "MATCH (m:mRNA)<-[:has_homolog]-(h:homology_group) "
            f"WHERE h.group_version = {version} "
            "OPTIONAL MATCH (m)-[:has_go]->(g:GO) "
            "RETURN m.id AS mRNA_id, m.genome AS genome, m.sequence AS chr, "
            "g.id AS go_id, id(h) AS hg_id, h.group_version AS hg_ver "
            # "LIMIT 25"
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
            "LIMIT 25"
        )

        col_names = ["mRNA_id", "genome", "chr", "hg_id", "hg_ver"]

        # run transaction
        with self.driver.session() as session:
            session._query_pangenome(query, cols=col_names, out=out)

    # Static methods, much like class methods, are methods that are bound to a
    # class rather than its object.
    @staticmethod
    def _query_pangenome(
        tx,
        query: str,
        cols: list[chr],
        out: str,
    ) -> None:
        # file to store the output

        print("\t".join(cols))

        for record in tx.run(query):
            # print(record)
            vals = record.values(*cols)
            print("\t".join(map(str, vals)))


if __name__ == "__main__":
    uri = "bolt://eddy.bioinformatics.nl:1241"
    user = "neo4j"
    password = "neo4j"

    panQ = QueryPangenome(uri, user, password)

    panQ.driver.verify_connectivity()

    try:
        panQ.get_GO(out="pangenome_GO.tab", version=1)
    finally:
        panQ.close()
