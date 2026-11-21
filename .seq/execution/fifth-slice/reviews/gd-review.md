# The Genomes dataset (GD): a review (2026-11-21)

Recommendation 12 of [whole-designs.md](whole-designs.md). Read against the five drafts and state.md's decisions (three
formats built, the twelve recommendations). Findings ranked; each applied to the draft.

## Blocking

1. **Four views name dropped runtimes.** Synteny ribbons, both genomes on one ring, the dot plot and the liftover read
   PAF and chain; the decision moves the alignment to the `seq` format, written from the alignment truth in code, read
   whole and indexed by both ends, sided (RR, GF). The liftover is generated only (fetch 9 is dropped, and the row
   cited "fetch 8", which is the GTF).
2. **The letters view needs bgzipped FASTA, `.gzi` and 2bit**, all designed only. The toys' letters come from an
   in-memory adapter on the seed (as the heavy fixture's), lengths inline or `chrom-sizes`; one small plain FASTA (chrM
   and the HLA contig) exercises the built FASTA reader and the quoted name. This also drops about 330 kB of fixtures.
3. **Status "served" on views the decisions defer**: the Manhattan family (EC-214, 216, 217, 219), `split` and the
   keyed-results join (EC-224) go to the close-out slice; their generator steps (GWAS, results) and the heavy
   fixture's GWAS go with them.

## To fix

4. Open by name: GFF3 names only (no bigBed `extraIndex`).
5. MAF: cut, its row removed; EC-213 cut.
6. EC-221 to RR's matrix (values), EC-226 to themes (Ranges, deferred): the GAP row becomes two assignments.
7. Fetches: 1, 2, 7, 8 kept; 3, 4, 5, 6, 9, 10 cut, with their rows.
8. The methylation long form has a path: RR's records transform at open (the wide form in the `seq` format's matrix).
   "GF has no long-form reader" is answered.
9. The ring's stems reduce bins to bins (`coarsenBins`, RR), not `binPairs`; its square needs `navigator.switch` and
   `section.include` (GA).
10. The ClinVar slice keeps its `.tbi` whole: the per-view byte bound (`maxBytesPerLayer`, RR) is then tested at the fit
    from the real index without its data.
11. Gaps 1 to 4, 6 and 8 are settled by the decisions; 5 and 7 stay. "For the gate" and "Review, unresolved" go.

## Checked and sound

Seeded streams per concern, the byte-identical fixture test, the 1 MB and 2 MB bounds (now about 0.5 MB), links
generated at load, the 3 Mb quick representative, an edge case named by each view's test.
