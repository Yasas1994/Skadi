# Extras and Citation

## Helpful command-line snippets

### Contig length distribution of your input contigs

This snippet uses [seqkit](https://bioinf.shenwei.me/seqkit/) to extract sequence lengths, then pipes them to a histogram. The `hist` command is available from the [tqdm](https://github.com/tqdm/tqdm) Python package (`pip install tqdm`) or similar histogram utilities.

```bash
seqkit fx2tab -lg input.fasta \
  | awk -F "\t" '{print $4}' \
  | tail -n +2 \
  | hist -b 100 -s 10
```

### Length distribution of contigs in ANI output (after filtering)

```bash
cat testing/gut_jaeger/nuc/gut_jaeger_virus_seqs_fasta_genome_ani.tsv \
  | awk -F "\t" '$8 > 0.1' \
  | awk -F "\t" '{print $3}' \
  | tail -n +2 \
  | hist -b 100 -x
```

---

## Citation

If you use **skadi**, please cite the underlying tools:

```text
[MMSEQS2] MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets.
          Steinegger M & Söding J. 2017. Nat. Biotech. 35, 1026–1028.
          https://doi.org/10.1038/nbt.3988

[PRODIGAL] Prodigal: prokaryotic gene recognition and translation initiation site identification.
           Hyatt et al. 2010. BMC Bioinformatics 11, 119.
           https://doi.org/10.1186/1471-2105-11-119.

[TAXONKIT] TaxonKit: A practical and efficient NCBI taxonomy toolkit.
           Shen, W. & Ren, H. J. 2021. Genet. Genomics.
           https://doi.org/10.1016/j.jgg.2021.03.006.
```
