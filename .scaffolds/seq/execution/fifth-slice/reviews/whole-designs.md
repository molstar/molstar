# The genome designs as a whole: a review (2026-11-21)

After the three-format decision (state.md, Decided). The author's questions first.

## Can the Rondo ring be built on BED, GFF and VCF?

Yes, with no new runtime. Sizes: inline `assembly.sequences` (or `chrom-sizes` in `runtimes/bed`); cytobands: BED read
whole (`fields: { stain }`, the `bands` layer); links: BEDPE read whole or held pairs, then RR's `binPairs` and
`PairsReader`; genes and labels: GFF by region through tabix (RN still says bigBed: to change); variants: VCF by region
(no draft puts them on the ring). Lost: 2bit or `.fai` as an assembly source (inline sizes instead), bigBed genes (GFF),
bigBed whole-genome density (a host's bigWig).

## With streaming?

For the 1D layers, yes; links are never streamed at whole-genome zoom. Links are read whole and summarised in time
slices within `summaries.maxBytes` (16 MB, about 600,000 pairs): 12,000 is trivial; 1M passes the open bounds but not the
pyramid's and falls to binning tile by tile (held rows outside the store's cache, the text parse not sliced); real Hi-C
has no path. BEDPE through tabix is indexed by one end, so GF ignores the index; pairix (`.px2`) or a doubled BEDPE
would give either-end reads, but still read every tile at the fit. GFF and VCF at the fit read every tile of a listed
level (exact or absent): ClinVar at the fit is the whole file (90 to 170 MB) — tiles bounded, the view not.

## Seams

- GF still designs the dropped runtimes throughout, and additions only they need (`measure: 'coverage'`, `block_index`,
  `.gzi`, `pageSize`, `assembly.sizes`).
- Sided pairs and `map.source` (synteny, the dot plot, the liftover) now need the `seq` runtime to serve sided systems
  and pairs found from either end; LX assumes the chain runtime; GD still names PAF and chain.
- No owner builds BEDPE's index of both ends after a whole read.
- GF's lane order loses PF; LA waits on the `seq` runtime's sided reads.
- GA and RR differ on the square (a 2D summary vs reads by block); RN's `section.include` replaces GA's tension.
- RN's stems call `binPairs` (items to bins) for a bins-to-bins reduction.
- Open by name reads GFF3 whole within 1M rows: a human GENCODE GFF3 (about 3M lines) finds nothing.
- GD: the Manhattan family has no design; `split` and EC-224's join have no owner; EC-221 and EC-226 unassigned; the
  toys' 2bit and bgzipped FASTA dropped; GD never reviewed.
- io: `readUncompressed`, `parseGzi`, `GziIndex` serve only the dropped bgzipped FASTA.

## Questions (recommendation first)

1. Strip the six dropped runtimes from GF to a short "designed only" appendix; drop `measure`, `block_index`,
   `pageSize`, `assembly.sizes`? Yes.
2. Remove `gzi` from io's public surface? Yes.
3. The `seq` runtime serves sided systems and pairs found from either end (synteny, dot plot, liftover); LA after RG
   and FR? Yes.
4. BEDPE's index of both ends after a whole read explicit in `runtimes/bed`; held rows count against a byte budget? Yes.
5. A doubled BEDPE (PIF-like, through tabix) in `runtimes/bed`? No: 12,000 links is the milestone; pairix and a levelled
   pairs format noted as the Hi-C route for later.
6. A per-view byte bound for tabix levels (the tiles in view's `regionBytes` summed; past it "zoom in")? Yes.
7. The whole-genome remedy for 1D layers is a host's bigWig, not bigBed? Yes.
8. The Manhattan family (EC-214, 216, 217, 219), `split` and EC-224's join deferred to the close-out slice? Yes.
9. MAF cut; EC-221 to the matrix, EC-226 to themes? Yes.
10. Fetch 3 (the 1000 Genomes panel) cut, fetches 1, 2, 7, 8 kept? Yes.
11. RN's gene detail reads GFF; the stems' reduction a bins-to-bins `coarsenBins`? Yes.
12. GD reviewed before promotion? Yes.
