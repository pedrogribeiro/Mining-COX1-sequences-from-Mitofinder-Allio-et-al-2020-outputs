# Mining-COX1-sequences-from-Mitofinder-Allio-et-al-2020-outputs

Extract the COX1 barcode from a set of [MitoFinder](https://github.com/RemiAllio/MitoFinder)
runs, reading the per-contig `raw.gff` annotations and contig FASTAs directly (I had to use these
because Mitofinder did not generate its final *_NT.fasta files). The workflow has four main scripts
further explained below.

**Requires:** Python 3, [Biopython](https://biopython.org/).

## Scripts

**`extract_reference_cox1.py <reference.gb> <out.fasta>`**
Pulls COX1 from an annotated GenBank reference mitogenome, producing the reference
FASTA used by `check_cox1.py`.

**`extract_cox1.py <results_dir> <outdir>`**
Extracts COX1 for every sample with COX1 on a single contig; writes per-sample
files and a combined `all_COX1.fasta`. Samples with COX1 on multiple contigs are
skipped and reported for manual review.

**`check_cox1.py <results_dir> <reference_COX1.fasta> SAMPLE [...]`**
For the flagged multi-contig samples, reports each candidate contig's length,
internal stop codons, and % identity to a reference COX1. A genuine barcode is
full length (~1536–1539 bp), has no internal stops, and is highly similar to the
reference; NUMTs and fragments deviate.

**`resolve_cox1.py <results_dir> <outdir> <decisions_file>`**
Extracts COX1 from the contig you chose per sample and rebuilds `all_COX1.fasta`.
Decisions file, one line per sample:

example:
```
PSY_0173    contig_1
PSY_0180    contig_2
```

The combined file is rebuilt (not appended), so this is safe to re-run.

## Workflow

Note: results is an actual folder from Mitofinder

```bash
python3 extract_reference_cox1.py <reference.gb> <reference_COX1.fasta>
python3 extract_cox1.py <results> <output>
python3 check_cox1.py <results> <reference_COX1.fasta> <sample1 sample2 sampleN> ## here you will name your samples accordingly
python3 resolve_cox1.py <results> <output> cox1_decisions.tsv
```

Set `GENE` at the top of the scripts to target a different gene. Translation uses
genetic code table 5 (invertebrate mitochondrial).
