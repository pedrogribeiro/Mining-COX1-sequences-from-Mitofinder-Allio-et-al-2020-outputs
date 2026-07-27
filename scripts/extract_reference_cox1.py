#!/usr/bin/env python3
"""
Extract COX1 from an annotated GenBank reference mitogenome, for use as the
reference in check_cox1.py.

Matches the CDS whose /gene qualifier is COX1 (also accepts COI / COX I / COXI).

Usage:
  extract_reference_cox1.py <reference.gb> <out.fasta>
"""
import sys
from Bio import SeqIO

GENE_ALIASES = {"COX1", "COI", "COXI", "COX I"}

def main():
    gb_path, out_path = sys.argv[1], sys.argv[2]
    rec = SeqIO.read(gb_path, "genbank")
    for feat in rec.features:
        if feat.type != "CDS":
            continue
        name = feat.qualifiers.get("gene", [""])[0].upper().replace("_", "").strip()
        if name in {a.replace(" ", "") for a in GENE_ALIASES}:
            seq = feat.extract(rec.seq)
            header = "%s_COX1" % rec.id
            open(out_path, "w").write(">%s\n%s\n" % (header, seq))
            print("COX1: %d bp -> %s" % (len(seq), out_path), file=sys.stderr)
            return
    sys.exit("No COX1/COI CDS found in %s" % gb_path)

if __name__ == "__main__":
    main()
