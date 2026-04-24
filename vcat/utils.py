
from typing import List, Union, Dict, Any, OrderedDict, DefaultDict, Optional, Mapping
from tqdm import tqdm
import mmap
from collections import defaultdict
from pathlib import Path
# from collections import OrderedDict as ordereddict
from io import StringIO
import polars as pl
import pandas as pd
from taxopy.core import TaxDb
from taxopy import Taxon
import polars.selectors as cs
import numpy as np
from numpy.typing import NDArray

# pandas implementation: readable but slow
# def compute_ani(alns):
#     # added pandas support for checkv ani script
#     return round(sum(alns["alnlen"] * alns["fident"]) / sum(alns["alnlen"]), 2)
#     # return round(
#     #     sum(a["len"] * a["pid"] for a in alns) / sum(a["len"] for a in alns), 2
#     # )

# def merge_intervals(intervals):
#     starts = intervals[:,0]
#     ends = np.maximum.accumulate(intervals[:,1])
#     valid = np.zeros(len(intervals) + 1, dtype=bool)
#     valid[0] = True
#     valid[-1] = True
#     valid[1:-1] = starts[1:] >= ends[:-1]
#     return np.vstack((starts[:][valid[:-1]], ends[:][valid[1:]])).T

# def compute_cov(alns):
#     t = alns[["qstart", "qend"]].to_numpy()
#     # sort both axis of the interval tree and merge overlaps
#     merged =  merge_intervals(np.sort(t, axis=1)[np.argsort(t, axis=0)[:,0]])
#     # calculate query coverage
#     return np.round(np.diff(merged).sum()/alns["qlen"].iloc[0], 2)


# def ani_summary(infile, outfile, header):
#     mmseqs_nuc = pd.read_table(infile, names=header)
#     output = []
#     groups = mmseqs_nuc.groupby(["query", "target"])
#     for g in tqdm(groups.groups):
#         g_ = groups.get_group(g)
#         output.append(
#             {"query":g[0],
#             "target": g[1],
#             "coverage":compute_cov(g_),
#             "ani":compute_ani(g_),
#             "lineage": g_["taxlineage"].iloc[0],
#             "qlen": g_["qlen"].iloc[0],
#             "tlen": g_["tlen"].iloc[0]
#             }
#             )
#     ani = pd.DataFrame(output)
#     ani["score"] = ani["ani"] * ani["coverage"]
#     ani.iloc[ani.groupby("query")["score"].idxmax().values].to_csv(outfile, index=None, sep="\t")


# polars implementation: fast!

_PARENT_MAP: Mapping[str, Optional[str]] = {
    'species':   'subgenus',   # change to 'genus' if that's your convention
    'subgenus':  'genus',
    'genus':     'subfamily',
    'subfamily': 'family',
    'family':    'suborder',
    'suborder':  'order',
    'order':     'class',
    'class':     'subphylum',
    'subphylum': 'phylum',
    'phylum':    'kingdom',
    'kingdom':   'realm',      # top-most rank
    'realm':     None,
}

def get_rank_taxid(*, taxon, rank: str) -> int:
    """
    Return the taxid associated with a given rank, walking up the taxonomy
    if the exact rank is missing for this taxon.

    Parameters
    ----------
    taxon : object
        An object with attribute `rank_taxid_dictionary` mapping rank (str) -> taxid (int).
        Missing/unknown entries may be absent or set to -1.
    rank : str
        The taxonomic rank you are querying (e.g. 'species', 'genus').

    Returns
    -------
    int
        The taxid for the nearest available rank at or above `rank`, or -1 if none found.

    Raises
    ------
    ValueError
        If `rank` is not found in the internal rank mapping.
    """
    # normalize
    key = rank.strip().lower()

    # validate rank exists in our hierarchy map
    if key not in _PARENT_MAP:
        raise ValueError(f"Unknown taxonomic rank: {rank!r}")

    # fetch from the taxon's dictionary, climbing to parents if missing
    # Accept both missing and explicit -1 as "not available".
    while key is not None:
        taxid = getattr(taxon, "rank_taxid_dictionary", {}).get(key, -1)
        if taxid is not None and taxid != -1:
            return taxid
        key = _PARENT_MAP[key]

    # Reached the top with nothing found
    return -1
        

def merge_intervals_with_cutoff(intervals, cutoff=0):
    """
    Merge intervals if the distance between consecutive (sorted) ranges is
    less than or equal to `cutoff`.

    Parameters
    ----------
    intervals : array-like of shape (n, 2)
        Each row is [start, end]. Start may be > end; we'll sort each pair.
        Works with ints or floats.
    cutoff : float or int, default=0
        Maximum allowed gap between two ranges to still merge them.
        Gap is defined as next_start - current_merged_end.
        Overlapping ranges have negative or zero gap and always merge.

    Returns
    -------
        np.ndarray of shape (m, 2)
        Merged intervals, sorted by start.

    Notes
    -----
    * Distance/gap uses continuous semantics: gap = next_start - current_end.
      Example: [0, 2] and [3, 5] have gap 1. With cutoff >= 1 they merge.
    * If you want discrete (integer) semantics where touching intervals
      should merge (e.g., [0,2] and [3,5] merge with cutoff=0), pass cutoff=1.
    """
    arr = np.asarray(intervals, dtype=int).copy()
    if arr.size == 0:
        return arr.reshape(0, 2)

    # Ensure each interval is [min, max]
    arr.sort(axis=1)

    # Sort by start, then end for stability
    order = np.lexsort((arr[:, 1], arr[:, 0]))
    arr = arr[order]

    merged = []
    cur_start, cur_end = arr[0]

    for s, e in arr[1:]:
        # If within cutoff (including overlaps), extend the current interval
        gap = s - cur_end
        if gap <= cutoff:
            # merge
            if e > cur_end:
                cur_end = e
        else:
            # finalize current and start new
            merged.append([cur_start, cur_end])
            cur_start, cur_end = s, e

    merged.append([cur_start, cur_end])
    out = np.array(merged, dtype=arr.dtype)
    return out


def merge_intervals(
    intervals: NDArray,
    *,
    assume_sorted: bool = False,
    merge_touches: bool = True,
    normalize: str = "swap",          # "swap" | "drop" | "error" for start>end
    keep_empty: bool = True,          # keep intervals where start==end
) -> NDArray:
    """
    Merge overlapping (and optionally touching) intervals using NumPy.

    Parameters
    ----------
    intervals : (N, 2) array-like
        Each row is [start, end]. May contain start >= end.
    assume_sorted : bool, default False
        If False, intervals are stably sorted by (start, end) before merging.
    merge_touches : bool, default True
        If True, [a,b] and [b,c] are merged; otherwise kept separate.
    normalize : {"swap","drop","error"}, default "swap"
        Handling for rows with start > end:
          - "swap": swap endpoints (treat as [min, max]).
          - "drop": remove those rows.
          - "error": raise ValueError.
    keep_empty : bool, default True
        Keep zero-length intervals (start == end) after normalization.
        If False, they are removed before merging.

    Returns
    -------
    (M, 2) np.ndarray of merged, disjoint intervals.
    """
    arr = np.asarray(intervals)
    if arr.ndim != 2 or arr.shape[1] != 2:
        raise ValueError("`intervals` must be a (N, 2) array-like.")
    if arr.size == 0:
        return np.empty((0, 2), dtype=arr.dtype)

    if np.issubdtype(arr.dtype, np.floating) and np.isnan(arr).any():
        raise ValueError("`intervals` contains NaN values.")

    starts, ends = arr[:, 0], arr[:, 1]
    rev_mask = starts > ends
    if rev_mask.any():
        if normalize == "error":
            raise ValueError("Found rows with start > end.")
        elif normalize == "drop":
            keep = ~rev_mask
            starts, ends = starts[keep], ends[keep]
        elif normalize == "swap":
            s, e = np.minimum(starts, ends), np.maximum(starts, ends)
            starts, ends = s, e
        else:
            raise ValueError("normalize must be one of {'swap','drop','error'}")

    # zero-length handling
    if not keep_empty:
        keep = starts != ends
        starts, ends = starts[keep], ends[keep]

    if starts.size == 0:
        return np.empty((0, 2), dtype=arr.dtype)

    # sort if needed
    if not assume_sorted:
        idx = np.lexsort((ends, starts))
        starts, ends = starts[idx], ends[idx]

    # cumulative max trick over ends to detect group boundaries
    cummax_ends = np.maximum.accumulate(ends)

    # boundary condition depends on whether touches merge
    if merge_touches:
        separators = starts[1:] > cummax_ends[:-1]
    else:
        separators = starts[1:] >= cummax_ends[:-1]

    # valid mask (length N+1) with True at group boundaries
    n = starts.size
    valid = np.empty(n + 1, dtype=bool)
    valid[0] = True
    valid[-1] = True
    valid[1:-1] = separators

    merged_starts = starts[valid[:-1]]
    merged_ends   = cummax_ends[valid[1:]]  # cummax at group ends

    return np.column_stack((merged_starts, merged_ends))


def compute_cov(alns: List[pl.Series]) -> float:
    """Compute query coverage using polars and numpy"""

    t = np.stack((alns[0].to_numpy(), alns[1].to_numpy()), axis=1)
    merged = merge_intervals(np.sort(t, axis=1)[np.argsort(t, axis=0)[:, 0]])
    return np.round(np.sum(np.diff(merged)) / alns[2][0], 2)


# ** Process the grouped data **
def ani_summary(infile: Union[str, StringIO], 
                all: bool, 
                header: list,
                dbdir: str | None = None,
                level: str | None = None) -> Union[Any, Exception]:
    """
    Robust per-group ANI summary that avoids list[f64] issues in Polars by
    grouping via pandas and computing scalar aggregates per (query,target).
    """
   
    try:
        mmseqs_nuc = pl.read_csv(
            infile,
            has_header=False, 
            separator="\t", 
            new_columns=header
        )

        # cast important columns explicitly
        mmseqs_nuc = mmseqs_nuc.with_columns(
            pl.col("fident").cast(pl.Float64),
            pl.col("alnlen").cast(pl.Float64),
            pl.col("qlen").cast(pl.Int64),
            pl.col("tlen").cast(pl.Int64),
            pl.col("qstart").cast(pl.Int64),
            pl.col("qend").cast(pl.Int64),
        )

        # drop rows missing numeric values
        mmseqs_nuc = mmseqs_nuc.filter(
            pl.col("alnlen").is_not_null() & pl.col("fident").is_not_null()
        )

        if (dbdir and level):
            taxdb = TaxDb(
                    nodes_dmp=f"{dbdir}/ictv-taxdump/nodes.dmp",
                    names_dmp=f"{dbdir}/ictv-taxdump/names.dmp",
                    merged_dmp=f"{dbdir}/ictv-taxdump/merged.dmp",
                )
            # accession	accession.version	taxid	gi
            qtaxinfo = pl.read_csv(f"{dbdir}/VMR_latest/virus_genome.accession2taxid",
                                                has_header=True, 
                                                separator="\t",
                                                new_columns=["query", "query_version", "taxid", "gi" ])
            # function(Taxon, level)
            mmseqs_nuc = mmseqs_nuc.with_columns(
                pl.col("taxid").map_elements(lambda x : get_rank_taxid(taxon=Taxon(x, taxdb=taxdb), rank=level), return_dtype=pl.Int64).alias("ttaxrank")
            )
            # qtaxinfo = qtaxinfo.with_columns(
            #     pl.col("taxid").map_elements(lambda x : Taxon(x, taxdb=taxdb).alias("qtaxlineage"))
            # )
            qtaxinfo = qtaxinfo.with_columns(
                pl.col("taxid").map_elements(lambda x : get_rank_taxid(taxon=Taxon(x, taxdb=taxdb), rank=level), return_dtype=pl.Int64).alias("qtaxrank")
            )
            # add qtaxonomic info
            
            mmseqs_nuc = qtaxinfo.select(cs.by_name(["query","qtaxrank"])).join(
                mmseqs_nuc, on="query", how="right"
            )
            # filter level
            mmseqs_nuc = mmseqs_nuc.filter(
                pl.col("qtaxrank") != pl.col("ttaxrank")
            )

        # drop 
        # convert to pandas for reliable group iteration
        # pdf = mmseqs_nuc.to_pandas()
        rows = []
        # for (_, _), g in pdf.groupby(["query", "target"]):
        for _, gpl in  mmseqs_nuc.partition_by("query", "target", as_dict=True).items():
            #gpl = pl.from_pandas(g)

            # basic scalars
            q = gpl["query"][0]
            t = gpl["target"][0]
            qlen = int(gpl["qlen"][0]) if "qlen" in gpl.columns else 0
            tlen = int(gpl["tlen"][0]) if "tlen" in gpl.columns else 0
            taxlineage = gpl["taxlineage"][0] if "taxlineage" in gpl.columns else ""
            taxid = int(gpl["taxid"][0]) if "taxid" in gpl.columns and gpl["taxid"][0] is not None else None

            # numeric arrays
            alnlen = np.asarray(gpl["alnlen"].to_numpy(), dtype=float)
            fident = np.asarray(gpl["fident"].to_numpy(), dtype=float)

            alnlen_sum = alnlen.sum() if alnlen.size else 0.0
            if alnlen_sum == 0:
                ani = 0.0
            else:
                ani = float((alnlen * fident).sum() / alnlen_sum)
            ani = round(ani, 3)

            # compute qcov using existing compute_cov (expects list of Series)
            try:
                qcov = compute_cov([gpl["qstart"], gpl["qend"], gpl["qlen"]])
            except Exception:
                qcov = 0.0

            tani = round(ani * qcov, 4)

            rows.append(
                {
                    "query": q,
                    "target": t,
                    "qlen": qlen,
                    "tlen": tlen,
                    "taxlineage": taxlineage,
                    "taxid": taxid,
                    "qcov": qcov,
                    "ani": ani,
                    "tani": tani,
                }
            )

        out = pl.DataFrame(rows) if rows else pl.DataFrame(
            {"query":[],"target":[],"qlen":[],"tlen":[],"taxlineage":[],"taxid":[],"qcov":[],"ani":[],"tani":[]}
        )

        if not all:
            # keep best hit per query (highest tani)
            out_pdf = out.to_pandas()
            best_pdf = out_pdf.sort_values("tani", ascending=False).drop_duplicates("query", keep="first")
            return pl.from_pandas(best_pdf)
        else:
            return out

    except Exception as e:
        return e


# aai calculation code

def trim_lineage(taxid: int,
                 taxdb: TaxDb,
                 taxomic_level: str = "species"
                 ) -> int:
    # trims lineages to the level
    if taxdb.taxid2rank[taxid] == taxomic_level:
        return taxdb.taxid2parent[taxid]
    return taxid


def get_taxid2taxon_map(
    df: pl.DataFrame, 
    taxdb: TaxDb
) -> Dict[str, OrderedDict]:
    tmap = dict()
    for x in list(df["taxid"].unique()):
        tmap[x] = Taxon(x, taxdb=taxdb).rank_taxid_dictionary
        tmap[x] = Taxon(x, taxdb=taxdb).rank_taxid_dictionary
    return tmap


def cal_axi(
    df: pl.DataFrame, 
    rank: str, 
    threshold: float, 
    taxdb: TaxDb, 
    kind: str, 
    all:bool=False,
) -> tuple[pl.DataFrame, List]:
    tkind = f"t{kind}"
    # print("cal axi prefilter", df['seqid'].unique().len())
    # print("cal axi postfilter == -1", df.filter(pl.col(rank) != -1)["seqid"].unique().len())
    # print("cal axi postfilter == -1 ", df.filter(pl.col(rank) == -1)["seqid"].unique().len())
    # print("cal axi postfilter", df.filter(pl.col(rank) == -1 ).select(["seqid", rank ]))
    # one seq id can have proteins mapping to different taxids
    # only filterout targets that do not map to trank not qrank
    tmp = (
        df.filter(pl.col(rank) != -1)
        .group_by(["seqid", rank])
        .agg(
            ((pl.col("fident") * pl.col("alnlen")).sum() / pl.col("alnlen").sum())
            .round(3)
            .alias(kind),
            (pl.col("qlen").unique().len() / pl.first("qlens").list.len())
            .round(3)
            .alias("qcov"),
            pl.first("qseqlen"),
            pl.col("qlen").unique().len().alias("mgenes"),
            (pl.first("qlens").list.len()).round(3).alias("qgenes"),
        )
        .with_columns(
            (pl.col(kind) * pl.col("qcov")).round(3).alias(tkind),
        )
        # .group_by("seqid")
        # .agg(
        #     pl.all().sort_by(tkind).last(),
        # )
    )
    if not all:
        tmp = tmp.group_by("seqid").agg(
                pl.all().sort_by(tkind).last(),
            )
    prefilt = set(df["seqid"].unique())
    postfilt = set(tmp.filter(pl.col(tkind) >= threshold)["seqid"].unique())
    tmp2 = (
        tmp.sort(tkind)
        .filter(pl.col(tkind) >= threshold)
        .with_columns(
            pl.col(rank)
            .map_elements(lambda x: str(Taxon(x, taxdb=taxdb)), return_dtype=pl.String)
            .alias("taxlineage")
        )
        .with_columns(level=pl.lit(rank))
    )

    # return taxids with taxi < threshold and pl.col(rank) == -1
    # tmp.filter(pl.col(tkind) < threshold)["seqid"]
    return (tmp2.rename({rank: "taxid"}), list(prefilt - postfilt))


def axi_summary(
    input: Union[str, StringIO],
    gff: str,
    dbdir: str,
    header: list,
    thresholds: Dict[str, float],
    top_k: int,
    kind: str,
    all: bool,
    level: str | None = None,
) -> Union[pl.DataFrame, Exception]:
    # ps: optimized for speed not for readability

    taxdb = TaxDb(
        nodes_dmp=f"{dbdir}/ictv-taxdump/nodes.dmp",
        names_dmp=f"{dbdir}/ictv-taxdump/names.dmp",
        merged_dmp=f"{dbdir}/ictv-taxdump/merged.dmp",
    )
    # extract query lengths (genomic seqs)
    leninf = []
    h = ""
    s = 0
    with open(gff) as fh:
        for i in fh:
            if i.startswith("# Sequence Data:"):
                for j in i.rstrip("\n").strip("# Sequence Data: ").split(";", 2):
                    if j.startswith("seqlen="):
                        s = j.lstrip("seqlen=")
                    elif j.startswith("seqhdr="):
                        h = j.lstrip("seqhdr=").strip('"').split()[0]
                        h = j.lstrip("seqhdr=").strip('"').split()[0]
                leninf.append({"seqid": h, "qseqlen": int(s)})

    seqleninfo = pl.DataFrame(leninf)

    col_names = [
        "seqid",
        "source",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attributes",
    ]
    # extract ORF stats
    gff_df = pl.read_csv(
        gff,
        separator="\t",
        comment_prefix="#",
        has_header=False,
        new_columns=col_names,
        schema_overrides={"start": pl.Int64, "end": pl.Int64, "score": pl.Float64},
    )

    gff_df = gff_df.with_columns(
        (
            pl.col("seqid")
            + pl.col("attributes").map_elements(
                lambda x: f'_{x.split(";")[0].strip("ID=").split("_")[-1]}',
                return_dtype=pl.String,
            )
        ).alias("query")
    )
    gff_df = gff_df.join(seqleninfo, on="seqid", how="inner")

    mmseqs_axi = pl.read_csv(
        input, has_header=False, separator="\t", new_columns=header
    )
    # pick the top-k hits per query protein
    explode_cols = [c for c in mmseqs_axi.columns if c != "query"]
    mmseqs_axi = (
        mmseqs_axi.group_by("query")
        .agg(
            pl.all().top_k_by("fident", top_k),
        )
        .explode(explode_cols)
    )

    mmseqs_axi = gff_df.select(
        ~cs.by_name(["source", "type", "phase", "attributes"])
    ).join(mmseqs_axi, on="query", how="right")
    mmseqs_axi = mmseqs_axi.with_columns(pl.col("fident").fill_null(0))
    mmseqs_axi = mmseqs_axi.with_columns(pl.col("taxid").fill_null(1))
    mmseqs_axi = mmseqs_axi.with_columns(pl.col("alnlen").fill_null(0))

    mmseqs_axi = mmseqs_axi.with_columns(
        pl.col("taxid").map_elements(
            lambda x: trim_lineage(x, taxdb=taxdb), return_dtype=pl.Int64
        )
    )
    taxonmap = get_taxid2taxon_map(mmseqs_axi, taxdb=taxdb)

    mmseqs_axi = mmseqs_axi.with_columns(
        pl.col("taxid")
        .map_elements(lambda x: taxonmap[x].get("species", -1), return_dtype=pl.Int64)
        .alias("species")
    )
    mmseqs_axi = mmseqs_axi.with_columns(
        pl.col("taxid")
        .map_elements(lambda x: taxonmap[x].get("genus", -1), return_dtype=pl.Int64)
        .alias("genus")
    )
    mmseqs_axi = mmseqs_axi.with_columns(
        pl.col("taxid")
        .map_elements(lambda x: taxonmap[x].get("family", -1), return_dtype=pl.Int64)
        .alias("family")
    )
    mmseqs_axi = mmseqs_axi.with_columns(
        pl.col("taxid")
        .map_elements(lambda x: taxonmap[x].get("order", -1), return_dtype=pl.Int64)
        .alias("order")
    )
    mmseqs_axi = mmseqs_axi.with_columns(
        pl.col("taxid")
        .map_elements(lambda x: taxonmap[x].get("class", -1), return_dtype=pl.Int64)
        .alias("class")
    )
    mmseqs_axi = mmseqs_axi.with_columns(
        pl.col("taxid")
        .map_elements(lambda x: taxonmap[x].get("phylum", -1), return_dtype=pl.Int64)
        .alias("phylum")
    )
    mmseqs_axi = mmseqs_axi.with_columns(
        pl.col("taxid")
        .map_elements(lambda x: taxonmap[x].get("kingdom", -1), return_dtype=pl.Int64)
        .alias("kingdom")
    )
    seqid2qlens = mmseqs_axi.group_by("seqid").agg(
        pl.col("qlen").unique().alias("qlens")
    )
    mmseqs_axi = mmseqs_axi.join(seqid2qlens, on="seqid", how="left")
    if (dbdir and level):
            # add qtaxonomic info
            # filter level
            # accession	accession.version	taxid	gi
            qtaxinfo = pl.read_csv(f"{dbdir}/VMR_latest/virus_genome.accession2taxid",
                                    has_header=True, 
                                    separator="\t",
                                    new_columns=["seqid", "seqid_version", "taxid", "gi" ])
            # qtaxinfo = qtaxinfo.with_columns(
            #     pl.col("taxid").map_elements(lambda x : Taxon(x, taxdb=taxdb).alias("qlineage"))
            # )
            # if a taxrank is un-available, get the next higher rank 
            qtaxinfo = qtaxinfo.with_columns(
                pl.col("taxid").map_elements(lambda x : get_rank_taxid(taxon=Taxon(x, taxdb=taxdb), rank=level), 
                                             return_dtype=pl.Int64).alias("qtaxrank")
            )
            mmseqs_axi = mmseqs_axi.with_columns(
                pl.col("taxid").map_elements(lambda x : get_rank_taxid(taxon=Taxon(x, taxdb=taxdb), rank=level),
                                              return_dtype=pl.Int64).alias("ttaxrank")
            )
            # add qtaxonomic info
            # filter level
            mmseqs_axi = qtaxinfo.select(cs.by_name(["seqid","qtaxrank"])).join(
                mmseqs_axi, on="seqid", how="right"
            )
            mmseqs_axi = mmseqs_axi.filter(
                pl.col("qtaxrank") != pl.col("ttaxrank")
            )


    try:
        results = []
        unmatched = mmseqs_axi["seqid"].to_list()
        for i in list(thresholds.keys()):
            mmseqs_nuc_f = mmseqs_axi.filter(pl.col("seqid").is_in(unmatched))
            tmp, unmatched = cal_axi(
                mmseqs_nuc_f, rank=i, threshold=thresholds[i], taxdb=taxdb, kind=kind, all=all
            )
            results.append(tmp)
        return pl.concat(results)

    except Exception as e:
        return e


def index_m8(input: Union[str, Path],
             kind: str) -> DefaultDict[int, List]:
    """
    m8 files with blast results can be insanely large and may require a large amount of
    memory to load.
    to over come this, we index the m8 file and load it chunk wise.
    This approch is slow but effective :)
    """
    index = defaultdict(list)

    with open(input, "r+b") as f:
        mmapped_file = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        total_size = len(mmapped_file)

        pbar = tqdm(unit="lines")

        pos = 0
        while pos < total_size:
            end = mmapped_file.find(b"\n", pos)
            if end == -1:  # End of file
                break

            line = mmapped_file[pos:end]  # Process as bytes
            key = line.split(b"\t", 1)[0]
            if kind == "axi":
                key = key.rsplit(b"_", 1)[0]  # Extract key as bytes
            index[key.decode()].append(int(pos))  # Store file offset

            pos = end + 1  # Move to the next line
            pbar.update(1)

    pbar.close()

    pbar = tqdm(unit="records")
    index2 = defaultdict(list)
    for k, (_, byte_offsets) in enumerate(index.items()):
        pbar.update(1)
        index2[k] = byte_offsets

    return index2


def load_chunk(
    input: str, 
    index: DefaultDict[int, List], 
    recstart: int, 
    recend: int
) -> StringIO:
    """
    load a single chunck of indexed blast results file to memory
    """
    string_dump = ""
    with open(input, "r+b") as fh:
        mmapped_file = mmap.mmap(fh.fileno(), 0, access=mmap.ACCESS_READ)

        records = list(range(recstart, recend))
        # Retrieve lines at given offsets
        for rec in records:
            for pos in index[rec]:
                end = mmapped_file.find(b"\n", pos)  # Find the end of the line
                if end == -1:  # Handle last line case
                    end = len(mmapped_file)

                string_dump += mmapped_file[pos:end].decode() + "\n"

    return StringIO(string_dump)



# bbmap.sh in=ERR2368798_1.fastq in2=ERR2368798_2.fastq path=/media/ssd/ICTV-TaxonomyChallenge/vcat/databases/VMR_latest/bbmap_index build=1 out=ERR2368798.sam covstats=ERR2368798.stats idfilter=0.7
# samtools sort ERR2368798.sam -o ERR2368798.bam
# samtools index -M ERR2368798.bam
# less ERR2368798.stats  | awk -F "\t" '$2 > 1'
# samtools coverage ERR2368798.bam -r 'MT249221 [Mourilyan virus] VRL' --histogram
