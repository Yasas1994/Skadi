## Annotating reads

This workflow can be used to map reads to ICTV geomes. This might come in handy in instances where the virus abundance is
extreamly low (clinical samples).


```bash
skadi reads -in <reads1>.fastq [-in2 <reads2.fastq>] -o outdir
```

Additional bbmap and pileup options can be used to customize this workflow.

```bash 
Usage: skadi reads [OPTIONS] [SNAKEMAKE_ARGS]...

  SKADI: Sequence-based Knowledgebase for Annotation, Detection, and Identification
  application designed to  provide taxonomy annotations to virus contigs and
  mapping reads directly to virus genomes.

  usage (paired-end): skadi reads [OPTIONS] -i1 pair1.fastq -i2 pair2.fastq -o
  mapping_results.tsv

  usage (single-end): skadi reads [OPTIONS] -i1 pair1.fastq -o
  mapping_results.tsv

Options:
  -in, --input PATH    input read file1 to run skadi on  [required]
  -in2, --input2 PATH  input read file2 (paired-end) to run skadi on
  -o, --output PATH    dir to store skadi results  [required]
  -j, --jobs INTEGER   use at most this many jobs in parallel (see cluster
                       submission for more details).  [default: 28]
  --bbmap_args TEXT    Extra arguments passed directly to BBMap (e.g.
                       'minid=0.95 maxindel=3')  [default: nodisk=t minid=0.95
                       maxindel=2 -Xmx20g]
  --pileup-args TEXT   Extra arguments passed to pileup.sh / CoveragePileup
                       (e.g. 'minmapq=20 minbaseq=20 mincov=2 secondary=f
                       delcov=f')  [default:  qtrim=t trimq=10 border=5
                       secondary=f delcoverage=f minmapq=20]
  --summary-args TEXT  coverage percentage = cpAverage fold = afMin total
                       reads = mtr-cp 50.0 af 1 mtr 100  [default:  -cp 0.5
                       -af 1 -mtr 100]
  --profile TEXT       snakemake profile e.g. for cluster execution.
  -n, --dryrun         Test execution.
  -h, --help           Show this message and exit.
  ```

**Results**

Results can be found in the `results` directory within the output directory. A typical layout is:

```text

├── input.coverage.tsv
├── input.pileup.tsv
├── input.sorted.bam
└── input.sorted.bam.bai
```
- `input.skadi.tsv` contains summary of all reads mapped to reference virus genomes.
- `input.pileup.tsv` contains pileup output.
- `input.sorted.bam` sorted bam file.
- `input.sorted.bam.bai` index of the sorted bam file.

