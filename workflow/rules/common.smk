import csv


with open(config["samples"]) as fh:
    SAMPLE_TABLE = {
        row["sample"]: row
        for row in csv.DictReader(
            (line for line in fh if not line.startswith("#")), delimiter="\t"
        )
    }
SAMPLES = list(SAMPLE_TABLE)


wildcard_constraints:
    sample="|".join(SAMPLES),


def ena_url(run, mate):
    """ENA FASTQ URL for a run accession (see ENA's directory layout docs)."""
    base = f"https://ftp.sra.ebi.ac.uk/vol1/fastq/{run[:6]}"
    digits = len(run) - 3
    if digits > 6:
        base += "/" + run[9:].zfill(3)
    return f"{base}/{run}/{run}_{mate}.fastq.gz"
