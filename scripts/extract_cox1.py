#!/usr/bin/env python3
"""
Extract COX1 from MitoFinder results.

Rule: for each sample, find every contig whose raw.gff annotates COX1.
  - Exactly one contig  -> extract COX1 (from whichever contig) and write it.
  - More than one contig -> skip; report the sample and contigs for manual review
                            (use check_cox1.py to decide which copy is the barcode).
  - Zero                 -> skip; report as no COX1.

After processing, all successfully written single-contig sequences are
concatenated into <outdir>/all_COX1.fasta.

Usage:
  extract_cox1.py <results_dir> <outdir>
"""
import sys, glob, os
from Bio import SeqIO

GENE = "COX1"          # change to extract a different gene
TABLE = 5              # invertebrate mitochondrial genetic code (kept for clarity)

def find_gene_hits(final_results_dir, gene):
    """Return (gff, contig_fasta, start, end, strand) for each contig annotating `gene`."""
    hits = []
    for gff in sorted(glob.glob(os.path.join(final_results_dir, "*_mtDNA_contig*_raw.gff"))):
        for line in open(gff):
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene" or f[8].strip() != gene:
                continue
            con = gff.replace("_raw.gff", ".fasta")
            hits.append((gff, con, int(f[3]) - 1, int(f[4]), f[6]))
            break
    return hits

def extract_seq(contig_fasta, start, end, strand):
    rec = next(SeqIO.parse(contig_fasta, "fasta"))
    seq = rec.seq[start:end]
    if strand == "-":
        seq = seq.reverse_complement()
    return seq

def main():
    results_dir, outdir = sys.argv[1], sys.argv[2]
    os.makedirs(outdir, exist_ok=True)

    written_samples = []
    multi, none, err = [], [], []

    for sample_dir in sorted(glob.glob(os.path.join(results_dir, "*.mitogenome"))):
        sample = os.path.basename(sample_dir).replace(".mitogenome", "")
        fr = glob.glob(os.path.join(sample_dir, "*_Final_Results"))
        if not fr:
            err.append((sample, "no Final_Results dir")); continue

        hits = find_gene_hits(fr[0], GENE)

        if len(hits) == 0:
            none.append(sample)
        elif len(hits) > 1:
            multi.append((sample, [os.path.basename(h[1]) for h in hits]))
        else:
            gff, con, start, end, strand = hits[0]
            if not os.path.exists(con):
                err.append((sample, "contig fasta missing: %s" % con)); continue
            seq = extract_seq(con, start, end, strand)
            out = os.path.join(outdir, "%s_%s.fasta" % (sample, GENE))
            open(out, "w").write(">%s_%s\n%s\n" % (sample, GENE, seq))
            written_samples.append(sample)

    combined = os.path.join(outdir, "all_%s.fasta" % GENE)
    with open(combined, "w") as out:
        for sample in written_samples:
            out.write(open(os.path.join(outdir, "%s_%s.fasta" % (sample, GENE))).read())

    print("=== %s extraction summary ===" % GENE, file=sys.stderr)
    print("written (single contig): %d" % len(written_samples), file=sys.stderr)
    print("skipped (multiple contigs): %d" % len(multi), file=sys.stderr)
    print("skipped (no %s): %d" % (GENE, len(none)), file=sys.stderr)
    print("errors: %d" % len(err), file=sys.stderr)
    print("combined file: %s (%d sequences)" % (combined, len(written_samples)), file=sys.stderr)
    for s, contigs in multi:
        print("  MULTI\t%s\t%s" % (s, ",".join(contigs)), file=sys.stderr)
    for s in none:
        print("  NONE\t%s" % s, file=sys.stderr)
    for s, msg in err:
        print("  ERROR\t%s\t%s" % (s, msg), file=sys.stderr)

if __name__ == "__main__":
    main()
