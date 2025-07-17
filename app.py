import streamlit as st
import pandas as pd
import gzip

# Read GTF (support .gz)
def read_gtf(file):
    cols = ["seqname", "source", "feature", "start", "end",
            "score", "strand", "frame", "attribute"]
    if file.name.endswith(".gz"):
        with gzip.open(file, "rt") as f:
            return pd.read_csv(f, sep="\t", comment="#", names=cols, low_memory=False)
    else:
        return pd.read_csv(file, sep="\t", comment="#", names=cols, low_memory=False)

# Parse attributes from GTF
def parse_attributes(attr_str):
    attrs = {}
    for item in attr_str.strip().split(";"):
        if item.strip():
            key, val = item.strip().split(" ", 1)
            attrs[key] = val.strip('"')
    return attrs

# Extract introns
def extract_introns(exons):
    exons["attributes"] = exons["attribute"].apply(parse_attributes)
    exons["transcript_id"] = exons["attributes"].apply(lambda d: d.get("transcript_id", "unknown"))
    intron_records = []
    for tid, group in exons.groupby("transcript_id"):
        group = group.sort_values("start")
        coords = list(zip(group.start, group.end))
        for i in range(len(coords) - 1):
            intron_start = coords[i][1] + 1
            intron_end = coords[i + 1][0] - 1
            if intron_end >= intron_start:
                intron_records.append({
                    "transcript_id": tid,
                    "seqname": group.iloc[0]["seqname"],
                    "start": intron_start,
                    "end": intron_end,
                    "strand": group.iloc[0]["strand"]
                })
    return pd.DataFrame(intron_records)

# Check overlap
def overlaps(a_start, a_end, b_start, b_end):
    return max(a_start, b_start) <= min(a_end, b_end)

# Find introns overlapping brain regions
def introns_in_brain_regions(introns, brain_regions):
    result = []
    for idx, intron in introns.iterrows():
        matching = brain_regions[
            (brain_regions["seqname"] == intron["seqname"]) &
            (brain_regions.apply(lambda row: overlaps(intron["start"], intron["end"], row["start"], row["end"]), axis=1))
        ]
        for _, region in matching.iterrows():
            result.append({
                **intron.to_dict(),
                "brain_region_source": region["source"],
                "brain_region_start": region["start"],
                "brain_region_end": region["end"]
            })
    return pd.DataFrame(result)

# Main app
def main():
    st.title("🧠 Brain Region Intron Filter (GTF-based)")

    gtf_file = st.file_uploader("Upload full genome GTF or GTF.GZ file", type=["gtf", "gtf.gz"])
    brain_gtf_file = st.file_uploader("Upload brain region GTF or GTF.GZ file", type=["gtf", "gtf.gz"])
    intron_list_file = st.file_uploader("Upload known introns list (TSV/CSV)", type=["tsv", "csv"])

    if gtf_file and brain_gtf_file and intron_list_file:
        with st.spinner("Processing..."):
            # Parse genome GTF
            gtf = read_gtf(gtf_file)
            exons = gtf[gtf["feature"] == "exon"]
            introns = extract_introns(exons)

            # Parse brain GTF and get regions
            brain_gtf = read_gtf(brain_gtf_file)
            brain_regions = brain_gtf[brain_gtf["feature"] != "exon"][["seqname", "source", "start", "end"]]

            # Filter introns overlapping brain regions
            introns_with_brain = introns_in_brain_regions(introns, brain_regions)

            # Load known introns list
            if intron_list_file.name.endswith(".tsv"):
                known_introns = pd.read_csv(intron_list_file, sep="\t")
            else:
                known_introns = pd.read_csv(intron_list_file)

            # Match introns by coordinate
            filtered = pd.merge(
                known_introns,
                introns_with_brain,
                on=["seqname", "start", "end"],
                how="inner"
            )

        st.success(f"Found {len(filtered)} known introns overlapping brain regions.")
        st.dataframe(filtered)

        out_tsv = filtered.to_csv(sep="\t", index=False)
        st.download_button("Download brain-region introns list", out_tsv, "brain_region_introns.tsv")
    else:
        st.info("Please upload all three required files to begin.")

if __name__ == "__main__":
    main()
