#!/usr/bin/env python3
"""
Inspect COX1 candidates for samples where COX1 is annotated on >1 contig.

For each named sample, finds every contig carrying COX1, extracts it, and reports:
  - length (bp)             -> expected ~1536-1539 for a full COX1
  - internal stop codons    -> 0 for a genuine CDS (table 5)
  - % amino-acid identity to a reference COX1 protein (local alignment)

Usage:
  check_cox1.py <results_dir> <reference_COX1.fasta> <SAMPLE> [SAMPLE ...]
"""
import sys, glob, os
from Bio import SeqIO, Align

TABLE = 5

_aligner = Align.PairwiseAligner()
_aligner.mode = "local"
_aligner.open_gap_score = -10
_aligner.extend_gap_score = -0.5
_aligner.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")

def identity_to_ref(prot, ref_prot):
    """% identity over aligned columns from a local protein alignment.
    Offset-independent: a true COX1 scores high regardless of where it starts."""
    prot = prot.rstrip("*").replace("*", "X")
    if not prot or not ref_prot:
        return 0.0
    aln = _aligner.align(prot, ref_prot)[0]
    ta, tb = aln[0], aln[1]
    same    = sum(1 for x, y in zip(ta, tb) if x == y and x != "-")
    aligned = sum(1 for x, y in zip(ta, tb) if x != "-" and y != "-")
    return 100.0 * same / aligned if aligned else 0.0

def main():
    results_dir = sys.argv[1]
    ref_path    = sys.argv[2]
    samples     = sys.argv[3:]

    ref_prot = str(next(SeqIO.parse(ref_path, "fasta")).seq.translate(table=TABLE)).rstrip("*")

    for sample in samples:
        fr = glob.glob(os.path.join(results_dir, sample + ".mitogenome", "*_Final_Results"))
        if not fr:
            print("%s: NO Final_Results" % sample); continue
        print("=== %s ===" % sample)
        print("  %-10s %8s %6s %8s" % ("contig", "len(bp)", "stops", "%id"))

        for gff in sorted(glob.glob(os.path.join(fr[0], "*_mtDNA_contig*_raw.gff"))):
            for line in open(gff):
                f = line.rstrip("\n").split("\t")
                if len(f) < 9 or f[2] != "gene" or f[8].strip() != "COX1":
                    continue
                con = gff.replace("_raw.gff", ".fasta")
                rec = next(SeqIO.parse(con, "fasta"))
                start, end, strand = int(f[3]) - 1, int(f[4]), f[6]
                seq = rec.seq[start:end]
                if strand == "-":
                    seq = seq.reverse_complement()
                prot  = str(seq.translate(table=TABLE))
                stops = prot[:-1].count("*")
                ident = identity_to_ref(prot, ref_prot)
                cname = os.path.basename(con).replace(".fasta", "").split("_mtDNA_")[-1]
                print("  %-10s %8d %6d %7.1f" % (cname, len(seq), stops, ident))
                break
        print()

if __name__ == "__main__":
    main()
