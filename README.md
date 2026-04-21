<picture>
  <source media="(prefers-color-scheme: dark)" srcset="./docs/_source/_static/figures/vcat_logo_dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="./docs/_source/_static/figures/vcat_logo_light.svg">
  <img alt="Fallback image description" src="./docs/_source/_static/figures/vcat_logo_dark.svg">
</picture>

![GitHub](https://img.shields.io/github/license/Yasas1994/vcat) 
![GitHub last commit (branch)](https://img.shields.io/github/last-commit/Yasas1994/vcat/main?color=8a35da)


The Virus contig annotation tool (Vcat) is a straightforward, homology-based application designed to provide taxonomic annotations for virus contigs, mapping reads to virus contigs and much more.

> [!NOTE] 
> Vcat [documentation](https://vcat.readthedocs.io/en/latest/
) is now online! 

### Changelog
---
<b>0.0.4</b>
* VMR_41v1 is now available to download
* added protein-level read annotation workflow 

<b>0.0.3</b>
* VMR_40v2 is now available to download
* added read annotation workflow 

<b>0.0.2</b>
* added downloaddb commad to download pre-built databases (VMR_40v1 and VMR_39v4)
* added a definition file to build Apptainer containers
* added subroutines to clean up tmp files generated duting the database build

---

### Quick start

```
git clone https://github.com/Yasas1994/vcat.git
cd vcat

# create a conda env
 mamba create -f environment.yml

# install vcat pipeline
pip install .

# test the installation
vcat --help

Usage: vcat [OPTIONS] COMMAND [ARGS]...

  Vcat: a command-line tool-kit for adding ICTV taxonomy annotations to virus
  contigs, mapping reads to virus genomes and much more.
  (https://github.com/Yasas1994/vcat)

Options:
  --version   Show the version and exit.
  -h, --help  Show this message and exit.

Commands:
  contigs    run contig annotation workflow
  downloaddb download pre-built reference databases
  preparedb  download and build reference databases
  reads      run read annotation workflow
  utils      tool chain for calculating ani, aai and visualizations

```
### Singularity (Now Apptainer)

```
git clone https://github.com/Yasas1994/vcat.git
cd vcat

# build container
apptainer build vcat.sif Apptainer

# test the container build
apptainer run vcat.sif vcat --help

Usage: vcat [OPTIONS] COMMAND [ARGS]...

  Vcat: a command-line tool-kit for adding ICTV taxonomy annotations to virus
  contigs, mapping reads to virus genomes and much more.
  (https://github.com/Yasas1994/vcat)

Options:
  --version   Show the version and exit.
  -h, --help  Show this message and exit.

Commands:
  contigs    run contig annotation workflow
  downloaddb download pre-built reference databases
  preparedb  download sequences and build reference databases
  reads      run read annotation workflow
  utils      tool chain for calculating ani, aai and visualizations

```
---
### downloading pre-built databases
---

```
# pulling pre-built databases from remote server [msl39v and masl40v1]
vcat downloaddb --dbversion masl39v4 -d <path to save the database> --cores 1
```

---

### downloading and preparing the databases

---

```
vcat preparedb -d <path to save the database>
```

---

### running contig annotation pipeline

---

```
vcat contigs -i <input>.fasta -o outdir
```

results can be found in the results directory within the ouput directory
```
.
├── logs
├── nuc
│   ├── input_genome_ani.tsv
│   └── input_genome.m8
├── prof
│   ├── input_fasta_prof_api.tsv
│   └── input_fasta_prof.m8
├── prot
│   ├── input_fasta.faa
│   ├── input_fasta.gff
│   ├── input_fasta_prot_aai.tsv
│   └── input_fasta_prot.m8
├── results
│   ├── *input_fasta_ictv.csv
│   └── input_fasta.tsv
└── tmp
```

#### Expected runtime ?

It takes ~4hrs to run vcat on the ICTV Taxonomy challenge dataset on a laptop computer.

### running read annotation pipeline

---

```
vcat reads -in <reads1>.fastq [-in2 <reads2.fastq>] -o outdir
```

---
### running other workflows
---
```
# calculate aai of query contigs to ICTV genomes
vcat utils aai  [OPTIONS] -i contigs.m8 -g configs.gff -d [DBDIR]

# calculate ani of query contigs to ICTV genomes
vcat utils ani [OPTIONS] -i contigs.m8

# creating genome comparision plots. i.e query sequence to highly similar ICTV genomes (comming soon)
vcat utils visualize --ani --taxa [taxname] -i contigs.m8 -o outdir

# create phage contig annotation plots (coming soon)
vcat utils visualize --phrogs -i contigs.fasta -o outdir 

# identify provirus (coming soon)
vcat utils provirus -i contigs.fasta -o outdir
```
---
#### Some additional stuff
---
```
# to view the contig length distribution of your contigs
seqkit fx2tab -lg {input.fasta} | awk -F "\t" '{print $4}' | tail -n +2 | hist -b 100 -s 10

# to view the length distribution of contigs in the ani calculation output (after apply a filter)
cat testing/gut_jaeger/nuc/gut_jaeger_virus_seqs_fasta_genome_ani.tsv | awk -F "\t" '$8 > 0.1' | awk -F "\t" '{print $3}' | tail -n +2 | hist -b 100 -x

```


---
If you use vcat please cite,

```

[MMSEQS2] MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets.
          Steinegger M & Söding J. 2017. Nat. Biotech. 35, 1026–1028. https://doi.org/10.1038/nbt.3988

[PRODIGAL] Prodigal: prokaryotic gene recognition and translation initiation site identification.
           Hyatt et al. 2010. BMC Bioinformatics 11, 119. https://doi.org/10.1186/1471-2105-11-119.

[TAXONKIT] TaxonKit: A practical and efficient NCBI taxonomy toolkit.
           Shen, W. & Ren, H. J. 2021. Genet. Genomics https://doi:10.1016/j.jgg.2021.03.006.

```
