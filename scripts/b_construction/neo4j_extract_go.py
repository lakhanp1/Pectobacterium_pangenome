#!/usr/bin/env python3

from neo4j import GraphDatabase
# from neo4j.exceptions import ServiceUnavailable


def get_GO(tx):

    query = (
        "MATCH (m:mRNA)<-[:has_homolog]-(h:homology_group) "
        "WHERE h.group_version = 1 "
        "OPTIONAL MATCH (m)-[:has_go]->(g:GO) "
        "RETURN m.id AS mRNA_id, m.genome AS genome, m.sequence AS chr, "
        "g.id AS go_id, id(h) AS hg_id, h.group_version AS hg_ver"
    )

    colNames = ['mRNA_id', 'genome', 'chr', 'go_id', 'hg_id', 'hg_ver']
    print("\t".join(colNames))

    for record in tx.run(query):
        # print(record)
        vals = record.values(*colNames)
        print("\t".join(map(str, vals)))


if __name__ == "__main__":
    URI = "bolt://eddy.bioinformatics.nl:1241"
    AUTH = ("neo4j", "neo4j")
    driver = GraphDatabase.driver(URI, auth=AUTH)

    # verify connection
    driver.verify_connectivity()

    with driver.session() as session:
        session.execute_read(get_GO)

    driver.close()