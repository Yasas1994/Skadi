from pathlib import Path
import os

# Set the input file and output directory as global variables
configfile: "config.yaml"
# this values are incuded in the snakefile
DBDIR = config["database_dir"]

rule all:
    input:
        f"{DBDIR}/ictv-taxdump.tar.gz",
        f"{DBDIR}/ictv-taxdump"

# Download the ICTV file
rule download_file:
    output:
        f"{DBDIR}/ictv.xlsx"
    params:
        url = config['ictv_url'],
	outdir = DBDIR,
	outfile = "ictv.xlsx"
    shell:
        """
        aria2c -x 4 {params.url}  --dir {params.outdir}  -o {params.outfile}
        """

# Convert XLSX to TSV
rule convert_xlsx_to_tsv:
    input:
        f"{DBDIR}/ictv.xlsx"
    output:
        f"{DBDIR}/ictv.tsv"
    params:
        sheet=f"{config['sheet']}"  # Sheet name from config
    shell:
        """
        csvtk xlsx2csv {input} -i {params.sheet} | csvtk csv2tab > {output}
        """

# Remove M-BM- characters
rule clean_special_characters:
    input:
        f"{DBDIR}/ictv.tsv"
    output:
        f"{DBDIR}/ictv.cleaned.tsv"
    shell:
        """
        sed -i 's/\xc2\xa0/ /g' {input}
        cp {input} {output}
        """

# Clean leading/trailing spaces and accidental newline characters
rule clean_whitespace:
    input:
        f"{DBDIR}/ictv.cleaned.tsv"
    output:
        f"{DBDIR}/ictv.clean.tsv"
    shell:
        """
        csvtk replace -t -F -f "*" -p "^\s+|\s+$" {input} \
        | csvtk replace -t -F -f "*" -p "\n " -r "" > {output}
        """

# Choose columns, rename, and remove duplicates
rule process_taxonomy:
    input:
        f"{DBDIR}/ictv.clean.tsv"
    output:
        f"{DBDIR}/ictv.taxonomy.tsv"
    shell:
        """
        csvtk cut -t {input} -f "Realm,Subrealm,Kingdom,Subkingdom,Phylum,Subphylum,Class,Subclass,Order,Suborder,Family,Subfamily,Genus,Subgenus,Species" \
        | csvtk rename -t -f 1- -n "realm,subrealm,kingdom,subkingdom,phylum,subphylum,class,subclass,order,suborder,family,subfamily,genus,subgenus,species" \
        | csvtk uniq -t -f 1- > {output}
        """

# Create taxdump directory
rule create_taxdump:
    input:
        f"{DBDIR}/ictv.taxonomy.tsv"
    output:
        directory(f"{DBDIR}/ictv-taxdump")
    shell:
        """
        taxonkit create-taxdump {input} --out-dir {output}
        """

# Create taxdump.tar.gz
rule create_taxdump_tar:
    input:
        f"{DBDIR}/ictv-taxdump"
    output:
        f"{DBDIR}/ictv-taxdump.tar.gz"
    shell:
        """
        tar -czvf {output} -C {input} names.dmp merged.dmp names.dmp nodes.dmp
        """


