#!/usr/bin/env python3
"""
Resolve multi-contig COX1 samples flagged by extract_cox1.py.

Reads a decisions file (tab- or whitespace-separated), one line per sample:
    PSY_0173    contig_1
    PSY_0180    contig_2
    PSY_0181    contig_1
Lines starting with # and blank lines are ignored.

For each sample it extracts COX1 from the chosen contig, writes
<outdir>/<SAMPLE>_COX1.fasta, then REBUILDS <outdir>/all_COX1.fasta from every
*_COX1.fasta present in <outdir>. The rebuild (not append) makes this safe to
run any number of times: the combined file always reflects exactly the
per-sample files on disk.

Usage:
  resolve_cox1.py <results_dir> <outdir> <decisions_file>
"""
import sys, glob, os
from Bio import SeqIO

GENE = "COX1"

def rebuild_combined(outdir, gene):
    """Rebuild all_<GENE>.fasta from every per-sample file in outdir. Idempotent."""
    combined = os.path.join(outdir, "all_%s.fasta" % gene)
    parts = sorted(f for f in glob.glob(os.path.join(outdir, "*_%s.fasta" % gene))
                   if os.path.basename(f) != "all_%s.fasta" % gene)
    with open(combined, "w") as out:
        for p in parts:
            out.write(open(p).read())
    return combined, len(parts)

def extract_from_contig(final_results_dir, contig_label, gene):
    """Find `gene` on the contig whose fasta name contains contig_label
    (e.g. 'contig_2'). Return the extracted Seq, or None if not found."""
    pattern = os.path.join(final_results_dir, "*_mtDNA_%s_raw.gff" % contig_label)
    for gff in sorted(glob.glob(pattern)):
        for line in open(gff):
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene" or f[8].strip() != gene:
                continue
            con = gff.replace("_raw.gff", ".fasta")
            if not os.path.exists(con):
                return None
            rec = next(SeqIO.parse(con, "fasta"))
            start, end, strand = int(f[3]) - 1, int(f[4]), f[6]
            seq = rec.seq[start:end]
            if strand == "-":
                seq = seq.reverse_complement()
            return seq
    return None

def main():
    results_dir, outdir, decisions = sys.argv[1], sys.argv[2], sys.argv[3]
    os.makedirs(outdir, exist_ok=True)

    resolved, failed = [], []
    for raw in open(decisions):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 2:
            failed.append((raw.rstrip("\n"), "malformed line")); continue
        sample, contig_label = parts[0], parts[1]

        fr = glob.glob(os.path.join(results_dir, sample + ".mitogenome", "*_Final_Results"))
        if not fr:
            failed.append((sample, "no Final_Results")); continue

        seq = extract_from_contig(fr[0], contig_label, GENE)
        if seq is None:
            failed.append((sample, "%s not found on %s" % (GENE, contig_label))); continue

        out = os.path.join(outdir, "%s_%s.fasta" % (sample, GENE))
        open(out, "w").write(">%s_%s\n%s\n" % (sample, GENE, seq))
        resolved.append((sample, contig_label, len(seq)))

    combined, n = rebuild_combined(outdir, GENE)

    print("=== %s resolution summary ===" % GENE, file=sys.stderr)
    print("resolved: %d" % len(resolved), file=sys.stderr)
    print("failed: %d" % len(failed), file=sys.stderr)
    print("combined rebuilt: %s (%d sequences)" % (combined, n), file=sys.stderr)
    for s, c, L in resolved:
        print("  OK\t%s\t%s\t%d bp" % (s, c, L), file=sys.stderr)
    for s, msg in failed:
        print("  FAIL\t%s\t%s" % (s, msg), file=sys.stderr)

if __name__ == "__main__":
    main()
