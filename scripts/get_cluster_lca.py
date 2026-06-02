#!/usr/bin/env python3
"""
author: Yasas Wijesekara (yasas.wijesekara@uni-greifswald.de)

This script extracts taxonomic information from the cluster members (used to build a profile)
and assigns each cluster to the LCA of all members.

Output: mmseqs_pprofiles_lca.tsv — used to postprocess profile-protein search results.

Columns:
    - cluster_representative (accession id)
    - root_p1 (1st taxid above the root, usually the realm)
    - taxid (taxid of the LCA of the cluster)
"""

import sys
import click
import pandas as pd
import taxopy

from skadi.color_logger import logger


def lca(taxa, taxdb, fraction=0.6):
    t = taxa.to_list()
    if len(t) > 1:
        return taxopy.find_majority_vote(t, taxdb=taxdb, fraction=fraction).taxid
    return taxa.iloc[0].taxid


@click.command(context_settings=dict(show_default=True))
@click.argument(
    "database_dir",
    type=click.Path(exists=True, file_okay=False, readable=True),
)
@click.option(
    "--fraction",
    default=0.6,
    type=float,
    help="Fraction of cluster members required for majority vote LCA.",
)
def main(database_dir: str, fraction: float) -> None:
    """Build cluster-to-LCA mapping from MMseqs2 protein clusters."""
    taxdb = taxopy.TaxDb(
        nodes_dmp=f"{database_dir}/ictv-taxdump/nodes.dmp",
        names_dmp=f"{database_dir}/ictv-taxdump/names.dmp",
        merged_dmp=f"{database_dir}/ictv-taxdump/merged.dmp",
    )

    acc2tx_path = f"{database_dir}/VMR_latest/virus_protein.accession2taxid"
    clust_path = f"{database_dir}/VMR_latest/mmseqs_pclusters/mmseqs_pclusters.tsv"
    out_path = f"{database_dir}/VMR_latest/mmseqs_pprofiles/mmseqs_pprofiles_lca.tsv"

    logger.info("Loading accession-to-taxid mapping: %s", acc2tx_path)
    acc2tx = pd.read_table(acc2tx_path)
    acc2tx = acc2tx.drop("gi", axis=1)

    logger.info("Loading cluster assignments: %s", clust_path)
    clusters = pd.read_table(clust_path, names=["target", "cluster_mem"])

    clusters = pd.merge(
        left=clusters,
        right=acc2tx,
        left_on="cluster_mem",
        right_on="accession.version",
    )

    logger.info("Building taxon objects…")
    clusters["lineage"] = clusters.apply(
        lambda x: taxopy.Taxon(x["taxid"], taxdb=taxdb), axis=1
    )
    clusters["root_p1"] = clusters.apply(
        lambda x: x["lineage"].taxid_lineage[-2], axis=1
    )

    logger.info("Computing LCAs…")
    clustlca = (
        clusters.groupby(["target", "root_p1"])["lineage"]
        .agg(lambda t: lca(t, taxdb=taxdb, fraction=fraction))
        .reset_index(name="taxid")
    )

    logger.info("Writing cluster-to-LCA mapping: %s", out_path)
    clustlca.to_csv(out_path, index=False, sep="\t")
    logger.info("Done.")


if __name__ == "__main__":
    main()
