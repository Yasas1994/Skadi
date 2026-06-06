# Annotating contigs

This workflow can be used to transfer taxonomic annotations to virus contigs. The workflow is written in 
Snakemake. So, if you interrupt it in the middle you can start again from where you left off. The following command
runs skadi with default parameters.

```bash
skadi contigs -i <input>.fasta -o outdir
```
## Advanced options

However, you can customize skadi with the following command line options. For example, consider setting `--nuc-search mmseqs_tblastx`
for more sensitive tANI calculations.

```bash
Options:
  -i, --input PATH                input file/s to run skadi on  [required]
  -o, --output PATH               dir to store skadi results  [required]
  -d, --database PATH             dir to skadi database
  --nuc_search [blastn|tbalstx|mmseqs_blastn|mmseqs_tblastx]
                                  nucleotide search algorithm (blastn,
                                  tblastx, mmseqs_blastn, mmseqs_tblastx).
                                  [default: mmseqs_blastn]
  --prot_search [blast|diamond|mmseqs]
                                  protein search algorithm (blast, diamond,
                                  mmseqs).  [default: mmseqs]
  --prof_search [mmseqs|hmmer]    profile search algorithm (mmseqs, hmmer).
                                  [default: mmseqs]
  -j, --jobs INTEGER              use at most this many jobs in parallel (see
                                  cluster submission for more details).
                                  [default: 8]
  --tapif FLOAT                   assign sequences above this tapi threshold
                                  to families  [default: 0.3]
  --tapio FLOAT                   assign sequences above this tapi threshold
                                  to orders  [default: 0.15]
  --tapic FLOAT                   assign sequences above this tapi threshold
                                  to classes  [default: 0.15]
  --tapip FLOAT                   assign sequences above this tapi threshold
                                  to phyla  [default: 0.15]
  --tapik FLOAT                   assign sequences above this tapi threshold
                                  to kingdoms  [default: 0.15]
  --taaig FLOAT                   assign sequences above this taai threshold
                                  to genera  [default: 0.49]
  --taaif FLOAT                   assign sequences above this taai threshold
                                  to families  [default: 0.3]
  --taaio FLOAT                   assign sequences above this taai threshold
                                  to orders  [default: 0.3]
  --taaic FLOAT                   assign sequences above this taai threshold
                                  to classes  [default: 0.3]
  --taaip FLOAT                   assign sequences above this taai threshold
                                  to phyla  [default: 0.3]
  --taaik FLOAT                   assign sequences above this taai threshold
                                  to kingdoms  [default: 0.3]
  --tanis FLOAT                   assign sequences above this tani threshold
                                  to species  [default: 0.81]
  --tanig FLOAT                   assign sequences above this tani threshold
                                  to genera  [default: 0.49]
  --batch INTEGER                 number of records to process at a time
                                  [default: 5000]
  -n, --dryrun                    Test execution.
  -h, --help                      Show this message and exit.
```

## Results

Results can be found in the `results` directory within the output directory. A typical layout is:

```text
outdir/
├── logs
├── nuc
│   ├── input_genome_ani.tsv
│   └── input_genome.m8
├── prof
│   ├── input_fasta_prof_api.tsv
│   └── input_fasta_prof.m8
├── prot
│   ├── input_fasta.faa
│   ├── input_fasta.gff
│   ├── input_fasta_prot_aai.tsv
│   └── input_fasta_prot.m8
├── results
│   ├── *input_fasta_ictv.csv
│   └── input_fasta.tsv
└── tmp
```
- `results/input_fasta.tsv` contains the putative taxonomic annotations for the input contigs.
- `nuc/input_genome_ani.tsv` contains all the tANI calculations before applying taxon-level tANI cutoffs.
- `prot/input_fasta_prot_aai.tsv` contains all txAAI calculations before applying taxon-level txAAI cutoffs.
- `prof/input_fasta_prof_api.tsv` contains all txAPI calculations before applying taxon-level txAPI cutoffs.
- To know more about tANI, txAAI and txAPI, refer to {doc}`annotation-pipeline`


```{note}
As a rough guideline, it takes around **4 hours** to run skadi on the ICTV Taxonomy challenge dataset on a typical laptop.
```
## SKADI output

skadi reports the predicted taxonomic lineage for each query sequence in a structured results table. Each row corresponds to a single sequence and lists the assigned taxa across all ICTV hierarchical ranks, from realm down to species, whenever a match is detected.

The Score column contains the computed similarity metric—tANI, txAAI, or txAPI, depending on which database level produced the annotation. These values range from 0 to 1, with 1.0 indicating a perfect match to a reference at that rank. The Method column specifies which comparison level generated the score:

ani → nucleotide-level comparison (tANI)

aai → protein-level comparison (txAAI)

api → profile HMM comparison (txAPI)

Higher levels of the table represent broader taxonomic classifications (e.g., Realm, Kingdom, Phylum), while lower levels provide increasingly specific annotations (e.g., Family, Genus, Species). Empty fields indicate that no confident match was found at that particular rank.

Together, these columns provide a detailed view of how confidently and at what resolution skadi was able to assign taxonomy to each query sequence.

```{csv-table} Taxonomic assignment summary
:header: SequenceID,Seqlen,Score,Method,Realm (-viria),Realm_score,Subrealm (-vira),Subrealm_score,Kingdom (-virae),Kingdom_score,Subkingdom (-virites),Subkingdom_score,Phylum (-viricota),Phylum_score,Subphylum (-viricotina),Subphylum_score,Class (-viricetes),Class_score,Subclass (-viricetidae),Subclass_score,Order (-virales),Order_score,Suborder (-virineae),Suborder_score,Family (-viridae),Family_score,Subfamily (-virinae),Subfamily_score,Genus (-virus),Genus_score,Subgenus (-virus),Subgenus_score,Species (binomial),Species_score
:widths: auto

EU188799,7498,1.0,ani,Monodnaviria,,,,Shotokuvirae,,,,Cossaviricota,,,,Papovaviricetes,,,,Zurhausenvirales,,,,Papillomaviridae,,Firstpapillomavirinae,,Dyoepsilonpapillomavirus,,,,Dyoepsilonpapillomavirus 1,
JQ245696,8601,1.0,ani,Riboviria,,,,Orthornavirae,,,,Kitrinoviricota,,,,Alsuviricetes,,,,Tymovirales,,,,Betaflexiviridae,,Quinvirinae,,Carlavirus,,,,Carlavirus americanense,
MK170446,43094,1.0,ani,Duplodnaviria,,,,Heunggongvirae,,,,Uroviricota,,,,Caudoviricetes,,,,,,,,Drexlerviridae,,,,Webervirus,,,,Webervirus wv13,
U01060,4615,1.0,ani,Riboviria,,,,Orthornavirae,,,,Duplornaviricota,,,,Chrymotiviricetes,,,,Ghabrivirales,,Alphatotivirineae,,Orthototiviridae,,,,Totivirus,,,,Totivirus ni,
KU687349,39444,1.0,ani,Duplodnaviria,,,,Heunggongvirae,,,,Uroviricota,,,,Caudoviricetes,,,,Autographivirales,,,,Autotranscriptaviridae,,Studiervirinae,,Kayfunavirus,,,,Kayfunavirus SH3,
MF352410,7407,1.0,ani,Riboviria,,,,Orthornavirae,,,,Pisuviricota,,,,Pisoniviricetes,,,,Picornavirales,,,,Picornaviridae,,Caphthovirinae,,Mischivirus,,,,Mischivirus ehoushre,
OL436139,42600,1.0,ani,Duplodnaviria,,,,Heunggongvirae,,,,Uroviricota,,,,Caudoviricetes,,,,,,,,Sarkviridae,,Guernseyvirinae,,Jerseyvirus,,,,Jerseyvirus LP31,
MK896629,6949,1.0,ani,Riboviria,,,,Orthornavirae,,,,Negarnaviricota,,Polyploviricotina,,Bunyaviricetes,,,,Elliovirales,,,,Peribunyaviridae,,,,Orthobunyavirus,,,,Orthobunyavirus bellavistaense,
MT074450,110183,0.999,ani,Duplodnaviria,,,,Heunggongvirae,,,,Uroviricota,,,,Caudoviricetes,,,,,,,,Demerecviridae,,Markadamsvirinae,,Epseptimavirus,,,,Epseptimavirus faergetype,
CASDWV010003421,3859,1.0,ani,Riboviria,,,,Orthornavirae,,,,Lenarviricota,,,,Leviviricetes,,,,Norzivirales,,,,Atkinsviridae,,,,Lobdovirus,,,,Lobdovirus lutivicinum,
AY780919,3095,1.0,ani,Riboviria,,,,Orthornavirae,,,,,,,,,,,,,,,,Birnaviridae,,,,Aquabirnavirus,,,,Aquabirnavirus salmonidae,
DQ163914,43365,1.0,ani,Duplodnaviria,,,,Heunggongvirae,,,,Uroviricota,,,,Caudoviricetes,,,,,,,,Fredfastierviridae,,,,Jamesmcgillvirus,,,,Jamesmcgillvirus jv119X,
CASDWV010000512,3611,1.0,ani,Riboviria,,,,Orthornavirae,,,,Lenarviricota,,,,Leviviricetes,,,,Norzivirales,,,,Solspiviridae,,,,Wulvevirus,,,,Wulvevirus geoadaptatum,
```


## Utility workflows

### Calculate txAAI of query contigs to ICTV genomes

```bash
skadi utils aai [OPTIONS] -i contigs.m8 -g contigs.gff -d [DBDIR]
```

### Calculate tANI of query contigs to ICTV genomes

```bash
skadi utils ani [OPTIONS] -i contigs.m8
```

### Visualization of genome comparisons (coming soon)

Create genome comparison plots (query sequence vs highly similar ICTV genomes):

```bash
skadi utils visualize --ani --taxa <taxname> -i contigs.m8 -o outdir
```

### Phage contig annotation plots (coming soon)

```bash
skadi utils visualize --phrogs -i contigs.fasta -o outdir
```

### Provirus identification (coming soon)

```bash
skadi utils provirus -i contigs.fasta -o outdir
```
