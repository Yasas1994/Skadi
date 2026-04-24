#!/usr/bin/env python
"""
author: Yasas Wijesekara (yasas.wijesekara@uni-greifswald.de)

This script extracts taxonomic information from the cluster members (used to build a profile)
add assigns each cluster to the lca of all members.

This information is written to mmseqs_pprofiles_lca.tsv -- later used to postprocess profile-protein
search results

mmseqs_pprofiles_lca.tsv incudes - cluster_representative (accession id), root_p1 (1st taxid above the root, usually the realm),
taxid (taxid of the lca of the cluster)

"""

import pandas as pd
import taxopy
from sys import argv
from vcat.color_logger import logger


# You can also use your own set of taxonomy files:
DATABASE_DIR = argv[1]

taxdb = taxopy.TaxDb(
    nodes_dmp=f"{DATABASE_DIR}/ictv-taxdump/nodes.dmp",
    names_dmp=f"{DATABASE_DIR}/ictv-taxdump/names.dmp",
    merged_dmp=f"{DATABASE_DIR}/ictv-taxdump/merged.dmp",
)

acc2tx = f"{DATABASE_DIR}/VMR_latest/virus_protein.accession2taxid"
clust = f"{DATABASE_DIR}/VMR_latest/mmseqs_pclusters/mmseqs_pclusters.tsv"

acc2tx = pd.read_table(acc2tx)
acc2tx = acc2tx.drop("gi", axis=1)

clusters = pd.read_table(clust, names=["target", "cluster_mem"])


def lca(taxa, taxdb=taxdb, fraction=0.6):
    t = taxa.to_list()
    if len(t) > 1:
        return taxopy.find_majority_vote(t, taxdb=taxdb, fraction=fraction).taxid
    return taxa.iloc[0].taxid


clusters = pd.merge(
    left=clusters, right=acc2tx, left_on="cluster_mem", right_on="accession.version"
)
clusters["lineage"] = clusters.apply(
    lambda x: taxopy.Taxon(x["taxid"], taxdb=taxdb), axis=1
)
clusters["root_p1"] = clusters.apply(lambda x: x["lineage"].taxid_lineage[-2], axis=1)

clustlca = (
    clusters.groupby(["target", "root_p1"])["lineage"]
    .agg(lca)
    .reset_index(name="taxid")
)
logger.info("writing cluster to LCA mapping file")
clustlca.to_csv(
    f"{DATABASE_DIR}/VMR_latest/mmseqs_pprofiles/mmseqs_pprofiles_lca.tsv",
    index=False,
    sep="\t",
)
