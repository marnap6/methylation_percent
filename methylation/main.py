import pandas as pd

# Wczytaj sekwencje
seqs = pd.read_csv("non_recombined_filtered_without_per.bed", sep="\t", header=None, names=["chrom", "start", "end", "count", "len", "strand"])
recombined = pd.read_csv("recombined_filtered_without_per.bed", sep="\t", header=None, names=["chrom", "start", "end"])
# Wczytaj metylację
meth = pd.read_csv("hmm_refined_meth_segments.bed", sep="\t", header=None, names=["chrom", "meth_start", "meth_end", "meth_pct", "state", "color"])

# Funkcja do obliczenia metylacji
def calculate_methylation(row):
    overlaps = meth[
        (meth["chrom"] == row["chrom"]) &
        (meth["meth_end"] > row["start"]) &
        (meth["meth_start"] < row["end"])
    ]

    total_overlap = 0
    methylated_bp = 0

    for _, mrow in overlaps.iterrows():
        overlap_start = max(row["start"], mrow["meth_start"])
        overlap_end = min(row["end"], mrow["meth_end"])
        overlap_length = max(0, overlap_end - overlap_start)
        total_overlap += overlap_length
        methylated_bp += overlap_length * (mrow["meth_pct"] / 100)

    if total_overlap > 0:
#        print(overlap_start, overlap_end, overlap_length, total_overlap, methylated_bp, round((methylated_bp / total_overlap) * 100, 3))
        return round((methylated_bp / total_overlap) * 100, 3)

    else:
        return 0.0

seqs["percent_methylation"] = seqs.apply(calculate_methylation, axis=1)
recombined["percent_methylation"] = recombined.apply(calculate_methylation, axis=1)


seqs[["chrom", "start", "end", "percent_methylation"]].to_csv("non_recombined_methylation_percent.bed", sep="\t", index=False)
recombined[["chrom", "start", "end", "percent_methylation"]].to_csv("recombined_methylation_percent.bed", sep="\t", index=False)