import os
# Set the input file and output directory as global variables
configfile: "config.yaml"
# this values are incuded in the snakefile
DBDIR = config["database_dir"]
CWD = os.path.dirname(__file__)
rule all:
    input:
        f"{DBDIR}/VMR_latest/virus_genome.accession2taxid",
        f"{DBDIR}/VMR_latest/virus_protein.accession2taxid",
        f"{DBDIR}/VMR_latest/proteins.faa",
        f"{DBDIR}/VMR_latest/genomes.fna",
        f"{DBDIR}/VMR_latest/genome2protein",
        f"{DBDIR}/VMR_latest/gff",
        f"{DBDIR}/VMR_latest/bbmap_index/index_done",
        f"{DBDIR}/VMR_latest/mmseqs_proteins/mmseqs_proteins",
        f"{DBDIR}/VMR_latest/diamond_proteins/diamond_proteins.dmnd",
        f"{DBDIR}/VMR_latest/mmseqs_genomes/mmseqs_genomes",
        f"{DBDIR}/VMR_latest/blast_genomes/blast_genomes",
        f"{DBDIR}/VMR_latest/mmseqs_pprofiles/mmseqs_pprofiles",
        f"{DBDIR}/VMR_latest/mmseqs_pprofiles/mmseqs_pprofiles_lca.tsv"

# Download the ICTV genomes
rule download_genbank:
    input:
        f"{DBDIR}/ictv.xlsx"
    output:
        outdir = directory(f"{DBDIR}/VMR_latest/tmp"),
        stamp = f"{DBDIR}/VMR_latest/tmp/download_complete"
    log:
        f"{DBDIR}/logs/download_genbank.log"
    shell:
        """
        download_gb.py -i {input} -o {output.outdir} &> {log}
        """
# generate gff files
rule generate_gff:
    input:
        gb_dir=f"{DBDIR}/VMR_latest/tmp"
    output:
        gff_dir=directory(f"{DBDIR}/VMR_latest/gff")
    log:
        f"{DBDIR}/logs/generate_gff.log"
    shell:
        r"""
        mkdir -p {output.gff_dir}
        for i in {input.gb_dir}/*.gb; do
            [ -e "$i" ] || continue
            gb2gff.py "$i" {output.gff_dir}
        done > {log} 2>&1
        """

# generate genome to protein map
rule generate_genome2protein:
    input:
        gb_dir=directory(f"{DBDIR}/VMR_latest/tmp")
    output:
        map=f"{DBDIR}/VMR_latest/genome2protein"
    log:
        f"{DBDIR}/logs/generate_genome2protein.log"
    shell:
        r"""
        genome_to_protein_map.py {input.gb_dir}/ -o {output.map} > {log} 2>&1
        """
# preprocess and prepare .gb files to build mmseqs datasets
rule prepare:
    input:
        DATADIR=f"{DBDIR}",
        XLTABLE=f"{DBDIR}/ictv.xlsx",
        MAN=f"{DBDIR}/VMR_latest/tmp/download_complete"

    output:
        f"{DBDIR}/VMR_latest/virus_genome.accession2taxid",
        f"{DBDIR}/VMR_latest/virus_protein.accession2taxid",
        f"{DBDIR}/VMR_latest/proteins.faa",
        f"{DBDIR}/VMR_latest/genomes.fna",

    log:
        f"{DBDIR}/logs/prepare.log"
    params:
        GBDIR=f"{DBDIR}/VMR_latest/tmp",
        SEQDIR=f"{DBDIR}/VMR_latest",
    shell:
        """
        prepare.py -d {input.DATADIR} -x {input.XLTABLE} -g {params.GBDIR} -s {params.SEQDIR} &> {log}
        """
     
rule make_mmseqs_proteindb:
    input:
        seqs=f"{DBDIR}/VMR_latest/proteins.faa",
        taxmap=f"{DBDIR}/VMR_latest/virus_protein.accession2taxid"
    output:
        f"{DBDIR}/VMR_latest/mmseqs_proteins/mmseqs_proteins"
    log:
        f"{DBDIR}/logs/make_mmseqs_proteindb.log"
    params:
        tmp=f"{DBDIR}/tmp/pdb",
        taxdump=f"{DBDIR}/ictv-taxdump",
        taxmap=f"{DBDIR}/VMR_latest/mmseqs_proteins/mmseqs_proteins_mapping"
    shell:
        """
        mkdir -p {params.tmp} &> {log}
        mmseqs createdb --dbtype 1 {input.seqs} {output}  &> {log}
        mmseqs createtaxdb {output} {params.tmp} --ncbi-tax-dump {params.taxdump} --tax-mapping-file {input.taxmap}  &>> {log}
        mmseqs nrtotaxmapping {input.taxmap} {output} {params.taxmap}  &>> {log}
        """

rule make_diamond_proteindb:
    input:
        seqs=f"{DBDIR}/VMR_latest/proteins.faa",
        taxmap=f"{DBDIR}/VMR_latest/virus_protein.accession2taxid"
    output:
        f"{DBDIR}/VMR_latest/diamond_proteins/diamond_proteins.dmnd"
    log:
        f"{DBDIR}/logs/make_diamond_proteindb.log"
    params:
        tmp=f"{DBDIR}/tmp/pdb",
        taxnodes=f"{DBDIR}/ictv-taxdump/nodes.dmp",
        taxnames=f"{DBDIR}/ictv-taxdump/names.dmp"
    shell:
        """
        diamond makedb --in {input.seqs} --db {output} --taxonmap {input.taxmap} --taxonnames {params.taxnames} --taxonnodes {params.taxnodes} > {log}
        """

rule make_mmseqs_genomedb:
    input:
        seqs=f"{DBDIR}/VMR_latest/genomes.fna",
        taxmap=f"{DBDIR}/VMR_latest/virus_genome.accession2taxid"
    output:
        f"{DBDIR}/VMR_latest/mmseqs_genomes/mmseqs_genomes"
    log:
        f"{DBDIR}/logs/make_mmseqs_genomedb.log"
    params:
        tmp=f"{DBDIR}/tmp/gdb",
        taxdump=f"{DBDIR}/ictv-taxdump",
        taxmap=f"{DBDIR}/VMR_latest/mmseqs_genomes/mmseqs_genomes_mapping"
    shell:
        """
        mkdir -p {params.tmp} &> {log}
        mmseqs createdb --dbtype 0 {input.seqs} {output}  &> {log}
        mmseqs createtaxdb {output} {params.tmp} --ncbi-tax-dump {params.taxdump} --tax-mapping-file {input.taxmap}  &>> {log}
        mmseqs nrtotaxmapping {input.taxmap} {output} {params.taxmap}  &>> {log}
        """

rule make_blast_genomedb:
    input:
        seqs=f"{DBDIR}/VMR_latest/genomes.fna",
        accmap=f"{DBDIR}/VMR_latest/virus_genome.accession2taxid"
    output:
        f"{DBDIR}/VMR_latest/blast_genomes/blast_genomes"
    log:
        f"{DBDIR}/logs/make_blast_genomedb.log"
    params:
        taxmap=f"{DBDIR}/VMR_latest/virus_genomes_mapping"
    shell:
        """
        cut -f 2,3 {input.accmap} | sed 1d > {params.taxmap}
        makeblastdb -in {input.seqs} -out {output} -parse_seqids -dbtype nucl -taxid_map {params.taxmap} &> {log} && \
        rm {params.taxmap} && \
        touch {output}
        """

rule make_bbmap_index:
    input:
        seqs=f"{DBDIR}/VMR_latest/genomes.fna",
    output:
        f"{DBDIR}/VMR_latest/bbmap_index/index_done"
    log:
        f"{DBDIR}/logs/make_bbmap_genomedb.log"
    params:
        build = 1,
        outdir=f"{DBDIR}/VMR_latest/bbmap_index"
    shell:
        """
        bbmap.sh ref={input} path={params.outdir} build={params.build} &> {log} && touch {output}

        """

rule make_mmseqs_profiledb:
    input:
        proteindb=f"{DBDIR}/VMR_latest/mmseqs_proteins/mmseqs_proteins",
        
    output:
        profiledb=f"{DBDIR}/VMR_latest/mmseqs_pprofiles/mmseqs_pprofiles",
        pclustdir=directory(f"{DBDIR}/VMR_latest/mmseqs_pclusters")

    log:
        f"{DBDIR}/logs/make_mmseqs_profiledb.log"
    params:
        tmp=f"{DBDIR}/tmp/pfdb",
        pclustdb=f"{DBDIR}/VMR_latest/mmseqs_pclusters/mmseqs_pclusters",
        clusters_tsv=f"{DBDIR}/VMR_latest/mmseqs_pclusters/mmseqs_pclusters.tsv",
        profiledir=f"{DBDIR}/VMR_latest/mmseqs_pprofiles",
        pclustdir=f"{DBDIR}/VMR_latest/mmseqs_pclusters",
        seqrep=f"{DBDIR}/VMR_latest/mmseqs_pclusters/sequenceRepDb",
        seqrep_h=f"{DBDIR}/VMR_latest/mmseqs_pclusters/sequenceRepDb_h",
        protein_h=f"{DBDIR}/VMR_latest/mmseqs_proteins/mmseqs_proteins_h",
        pclust_lookup='./sequenceRepDb.lookup',
        pclust_mapping='./sequenceRepDb_mapping',
        pclust_source='./sequenceRepDb.source',
        pclust_taxonomy= './sequenceRepDb_taxonomy',
        pprof_lookup='./sequenceRepDb.lookup',
        pprof_mapping='./sequenceRepDb_mapping',
        pprof_source='./sequenceRepDb.source',
        pprof_taxonomy='./sequenceRepDb_taxonomy',
        pprof_lookup_='./mmseqs_pprofiles.lookup',
        pprof_mapping_='./mmseqs_pprofiles_mapping',
        pprof_source_='./mmseqs_pprofiles.source',
        pprof_taxonomy_='./mmseqs_pprofiles_taxonomy',
        pprof_pprofiles_h = './mmseqs_pprofiles_h',
        pprof_pprofiles_h_dbtype = './mmseqs_pprofiles_h.dbtype',
        pprof_pprofiles_h_index = './mmseqs_pprofiles_h.index'
        
    shell:
        """
        mkdir -p {params.tmp} &> {log}
        mkdir -p {output.pclustdir} >> {log}
        mkdir -p {params.profiledir} >> {log}
        mmseqs cluster {input.proteindb} {params.pclustdb} {params.tmp} --cluster-reassign  -s 7 &> {log}
        mmseqs createtsv {input.proteindb} {input.proteindb} {params.pclustdb} {params.clusters_tsv} &>> {log}
        mmseqs createsubdb {params.pclustdb}  {input.proteindb} {params.seqrep}  &>> {log}
        mmseqs createsubdb {params.pclustdb}  {params.protein_h} {params.seqrep_h}  &>> {log}
        mmseqs result2profile {params.seqrep} {input.proteindb} {params.pclustdb} {output.profiledb}  &>> {log}
        cd {params.pclustdir} &>> {log}
        ln -sf ../mmseqs_proteins/mmseqs_proteins.lookup {params.pclust_lookup}  &>> {log}
        ln -sf ../mmseqs_proteins/mmseqs_proteins_mapping {params.pclust_mapping}  &>> {log}
        ln -sf ../mmseqs_proteins/mmseqs_proteins.source {params.pclust_source} &>> {log}
        ln -sf ../mmseqs_proteins/mmseqs_proteins_taxonomy {params.pclust_taxonomy} &>> {log}
        cd {params.profiledir} &>> {log}
        ln -sf ../mmseqs_proteins/mmseqs_proteins.lookup {params.pprof_lookup} &>> {log}
        ln -sf ../mmseqs_proteins/mmseqs_proteins_mapping {params.pprof_mapping} &>> {log}
        ln -sf ../mmseqs_proteins/mmseqs_proteins.source {params.pprof_source} &>> {log}
        ln -sf ../mmseqs_proteins/mmseqs_proteins_taxonomy {params.pprof_taxonomy} &>> {log}
        ln -sf ../mmseqs_proteins/mmseqs_proteins.lookup {params.pprof_lookup_} &>> {log}
        ln -sf ../mmseqs_proteins/mmseqs_proteins_mapping {params.pprof_mapping_} &>> {log}
        ln -sf ../mmseqs_proteins/mmseqs_proteins.source {params.pprof_source_} &>> {log}
        ln -sf ../mmseqs_proteins/mmseqs_proteins_taxonomy {params.pprof_taxonomy_} &>> {log}
        ln -sf ../mmseqs_pclusters/sequenceRepDb_h {params.pprof_pprofiles_h} &>> {log}
        ln -sf  ../mmseqs_pclusters/sequenceRepDb_h.dbtype {params.pprof_pprofiles_h_dbtype} &>> {log}
        ln -sf ../mmseqs_pclusters/sequenceRepDb_h.index {params.pprof_pprofiles_h_index} &>> {log}
        """


rule cluster_lca:
    input:
        DATADIR=DBDIR,
        profiledb=f"{DBDIR}/VMR_latest/mmseqs_pprofiles/mmseqs_pprofiles"

    output:
        f"{DBDIR}/VMR_latest/mmseqs_pprofiles/mmseqs_pprofiles_lca.tsv",

    log:
        f"{DBDIR}/logs/get_cluster_lca.log"

    shell:
        """
        get_cluster_lca.py {input} &> {log}
        """