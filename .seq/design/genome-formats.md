# The genome formats: a runtime each, BGZF and tabix in the base library

Lane GF's design (the fifth slice, stage 0; agreed by the author 2026-11-21 and promoted from its draft in `execution/fifth-slice/drafts/`). It amends
[files](files.md) (BED, GFF3's tension closed, `fields` by column number), [sources](sources.md)
(the genome runtimes), [base-library](base-library.md) (BGZF and tabix in `@molstar/io`),
[protocol](protocol.md) (an item's extent, a tile's read cost, names), [security](security.md)
(the formats' bounds), [performance](performance.md) (their rows) and RR's design
([region-reads](region-reads.md): sided pairs sources, items at a tile size for ranges and matrices) (below,
Consequences). Built on Gate 0's answers as recommended (the fifth slice's plan): BGZF and tabix in `@molstar/io`, beside its
readers, a public surface agreed here (2); generated fixtures first and one real slice a format, each fetch asked (7);
cytobands synthetic in tests, UCSC's asked (6). Built on GA's design ([genome-coordinates](genome-coordinates.md):
assemblies, the names rule, a map between assemblies as one map over sequence pairs, the locus form's fallback to a name
index) and RR's ([region-reads](region-reads.md): the `pairs`, `matrix` and `rows` shapes, items read at a tile size,
bounds in bytes, the one-end index rule). BAM and CRAM are out of scope (plan, Out of scope).

## Intent

sources.md: "binary and indexed genome formats are runtimes an application adds"; the plan: "each a runtime an
application adds, opt-in entries like `runtimes/mmcif`, each with a conformance suite on recorded or generated
fixtures: BED (`base: 0`, BED12 blocks as pieces), bigBed, 2bit, FASTA with its index by region, GFF3 and GTF (the
`Parent` hierarchy into transcripts, files.md's tension closed), VCF with tabix, PAF, chain files". Files.md left three
things for here by name: "Files read by region (tabix, BGZF) wait for the Genomes milestone", "**BED** waits for a
dataset that has one … BED12's blocks as pieces", and its Tension "GFF3 without its hierarchy … Open until a dataset
draws transcripts from GFF3". The JBrowse reply's first ask is the same list ("bigWig, bigBed, 2bit, VCF with a tabix
index, FASTA with fai, GFF3, … PAF", EC-185 to EC-190). What a user does with these files: open a genome by region, see
its genes as transcripts, its variants and genotypes, its gaps and letters, find a gene by name (EC-228), and compare two
genomes through an alignment's pairs (LX's figures). Cost follows pixels ([performance](performance.md)):
a region read costs the region, a view's reads are bounded together, every claim a file makes is bounded before it
allocates ([security](security.md)).

**Three runtimes are built**: `runtimes/bed`, `runtimes/gff` and `runtimes/vcf`, with BGZF and tabix in `@molstar/io`;
the built `seq` runtime gains sided pairs read from either end. bigBed, 2bit, indexed and bgzipped FASTA, PAF, PIF and
chain files are designed only (the appendix), no runtime until a use asks. Why: a runtime is maintenance and fixtures,
and a host brings its own through `withAdapter` (the author, 2026-11-21: "Yes, those three").

## Precedents

| Precedent | Taken | Not taken |
|---|---|---|
| htslib's specifications: BGZF (SAM spec section 4.1: gzip members of at most 64 KiB with a `BC` extra field, virtual offsets, the 28-byte EOF block), tabix (`.tbi`: bins, chunks, a 16 kb linear index, the column header, `0x10000` zero-based), CSI (`min_shift`, `depth`, a bin's `loffset`), `tbx_parse1`'s extent of a VCF line (`POS` to `POS + len(REF)`, `INFO/END` when present) | every rule as specified; CSI beside TBI (chromosomes past 512 Mb: wheat, axolotl); a line's extent by the index's own rule, so a region read agrees with the index at its edges | writing (bgzip, tabix): tests write their own fixtures with a test-only writer |
| htslib's VCF 4.4: header lines, `##contig` lengths, `##INFO`/`##FORMAT` `Number` and `Type`, `GT` with `/` and `|`, symbolic alleles, `*`, `<NON_REF>` | typed INFO fields by the header, a row per ALT allele (EC-196), genotypes as allele dosage | phasing (a Tension), BCF |
| UCSC: bigBed (the BBI container shared with bigWig; autoSql; `extraIndex` B+ trees on a field), 2bit (versions 0 and 1, N and mask blocks), BED (`bed N+M`, BED12 blocks, `thickStart`, `itemRgb`), chain (`chain` headers and `size dt dq` blocks, `qStrand -` counted on the reverse strand), liftOver | every container rule; a field's name by autoSql; `extraIndex` as the name index; BED's `N+M` as `fields` by column number; chains as the runs of one map (GA) | `itemRgb` as a colour without a theme (a field here, Tensions); bigBed's untyped extras guessed (autoSql or `fields` say) |
| minimap2's PAF (12 columns, 0-based, `cg:Z` CIGAR, `tp:A`, `cs:Z`) and JBrowse's PIF (a PAF indexed by both ends: each line twice, names prefixed `q` and `t`, bgzipped and tabix-indexed) | PAF read whole; PIF read by region from either end through tabix (RR's "pairs found from either end") | `cs:Z` (the CIGAR suffices for runs) |
| JBrowse 2's readers (`@gmod/bgzf-filehandle`, `@gmod/tabix`, `@gmod/vcf`, `@gmod/gff`, `@gmod/bbi`, `@gmod/twobit`, `@gmod/indexedfasta`; the PAF, PIF and chain adapters; trix name search) | an adapter a format; GFF's `Parent` assembled into nested features, discontinuous features by `ID`; genotypes parsed only when asked; bigBed's `searchExtraIndex`; bgzipped FASTA with `.fai` and `.gzi` | the byte-size estimate before a read ("force load", EC-171; RR refuses it too); trix (a later format, Gaps) |
| igv.js: tabix and BGZF of its own, a feature's `visibilityWindow` (VCF and BED zoomed out show "zoom in"), GFF and GTF combined into transcripts by `Parent` and `transcript_id`, an indel drawn without its anchor base | the visibility window as each tabix source's coarsest level; GTF's `transcript_id` and `gene_id` hierarchy | trimming the anchor base (here the record's REF extent, the index's rule: Decision, VCF) |
| The current Seq\* (`stream/gmod`): 2bit, bigWig, bigBed and tabix VCF over range reads, its bounds (a BGZF block at most 64 KB with its ISIZE checked, a `.tbi` run at most 256 MB, 2bit's mask blocks found by binary search on disk), ids `name@start-end` and `pos:ref>alt`, a VCF's `estimate` from the index's bytes, recorded ranges of hg38 files | the readers' rules and bounds, the ids, recording real files by the ranges read | the estimate drawn as density (what is drawn is checked against the data); a range recording keyed by its old reads (here core's pages) |
| This codebase: the bigWig runtime (BBI header, chromosome tree, R-tree, `walk`'s bounds, `BigWigLimits`, `DecompressionStream`), core's range reads (`openRanges`, 64 KiB pages, `AssetRangeError`), core's `inflateStream` and `ungzip`, the FASTA runtime (`.fai` checked), the records runtime and its GFF3 preset, `withAdapter`, the `seq` runtime's join, `delayedRuntime` | one adapter a format through `withAdapter`; the records' issue rules (EC-86) and `limits`; the `seq` runtime for sided pairs | a second inflater, a second page cache |

The bigBed, 2bit, FASTA, PAF, PIF and chain rows serve the appendix's designs.

### Interactions in the precedents

| What a user does there | Here |
|---|---|
| open a file by URL with its index (all) | built: a `source` with `uri` and `index` (`.tbi`, `.csi`), never `<uri>.tbi` assumed (EC-186) |
| see a gene's transcripts, exons, CDS and UTRs (JBrowse, IGV, UCSC) | built: GFF3 and GTF as transcripts (`features: transcripts`), drawn by the ranges layer's `line: gene` (GG) |
| collapse transcripts to genes, or expand (IGV's collapsed and expanded, JBrowse's display modes) | built: `features: genes` or `transcripts` by patch; `rows` for the file's lines |
| click a variant: its INFO in a popup (all) | built: the hover card reads the fields (built) |
| filter variants by an INFO value (`CLNSIG`, `AF > 0.01`) (JBrowse) | built: the layer's data filter (the query sublanguage, EC-225) |
| genotypes of many samples, sorted by a field (IGV) | built: RR's `matrix` layer over the VCF's matrix; `samples` narrows them (new) |
| "zoom in to see features" past a visibility window (IGV, JBrowse) | built: a tabix source's coarsest level (listed per file from its index's bytes) and a view's read cost (RR's `maxBytesPerLayer`), past either the layer's too-many-tiles state, never an estimate |
| search a gene by name (IGV, JBrowse's trix, UCSC) | built in data: `findNames` over a GFF3 or GTF read whole; the box and the jump are GA's locus fallback (GC3) |
| show N gaps and soft-masked repeats (UCSC) | designed only: a 2bit's N and mask blocks as `ranges` (appendix) |
| lift a region from one assembly to another (UCSC liftOver) | built in data: a pairs source's equal-length runs as a map (`map.source`: the `seq` runtime, a host's chain adapter); the figure LX's |
| colour items by `itemRgb` (UCSC) | refused here: a field; a theme reading an RGB text is themes' (Tensions) |
| export the features in view (JBrowse) | refused: data export is a host's |

## Decision

### Whose each choice is

By [axis](axis.md)'s rule. What a file means is the spec's (a source's params: its `index`, its
`fields`, its `features`, its `samples`); how much a read may cost is the runtime's (an option of
`<format>Runtime(over)`, as `bigWigRuntime`) or the store's (`limits`, `maxTileBytes`); a format's conventions (BED's
base 0, GFF's base 1) are fixed by the format, never a param a writer can get wrong (EC-195). No new port: each format
is a registry entry of the existing kind ([sources](sources.md), [seq-data-ports](seq-data-ports.md)).

### One runtime a format family, each an opt-in entry

| Entry (`@molstar/seqx-core/…`) | Exports | Formats | Shapes | Read |
|---|---|---|---|---|
| `runtimes/bed` | `BedRuntime`, `bedRuntime` | `bed`, `bedgraph`, `bedpe`, `chrom-sizes` | ranges; values; pairs; coverage only | whole, or by region with a tabix `index` (not `bedpe`, not `chrom-sizes`) |
| `runtimes/gff` | `GffRuntime`, `gffRuntime`, `assembleTranscripts` | `gff3` (in place of the records preset by `withAdapter`), `gtf` | ranges | whole, or by region with a tabix `index` |
| `runtimes/vcf` | `VcfRuntime`, `vcfRuntime`, `readVcfHeader` | `vcf` | ranges (a row an ALT allele); matrix (genotypes, FORMAT values) | by region with `index`; whole without |
| `seq` (built, amended) | as built | `seq` | pairs, sided, found from either end; grid; matrix | whole, indexed by both ends (below) |

- **Why an entry a family, opt in**: a reader of a spec, the minimal viewer and the default viewer pay for none
  ([performance](performance.md), Size: "a format runtime … loaded only when composed"; the plan: "the
  minimal viewer not grown"). A host lists what it reads, as `BigWigRuntime` today. Why families: BED, bedGraph, BEDPE
  and chrom.sizes share one tab-line reader and one row of size.
- **A format the default viewer already reads** (`gff3` flat, through the records preset) is taken over by the genome
  runtime when a host composes it with `withAdapter` (built: one adapter a format, the default narrowed, never
  order-dependent). Without it the built behavior stands. Why not change the defaults: the default viewer is at its
  size bound, and a plasmid's GFF3 needs no region read or hierarchy.
- **Internals shared, not published**: `runtimes/internal/lines.ts` (a tab-line source read whole or by region through
  tabix, below), `runtimes/internal/text.ts` (bytes to text: plain, gzip or BGZF by magic, bounded). The package
  exports `"./runtimes/internal/*": null`, as core's `util/internal`. Why: core's `./*` rule publishes every file
  otherwise.
- **Each with its synthetic twin** (`delayedRuntime(<Runtime>)` on a fake clock), as sources.md requires.

### BGZF and tabix in `@molstar/io` (Gate 0, 2): the public surface

Two modules beside io's containers (`src/common/`, where `binary-cif`, `msgpack` and `netcdf` are), types over core's
`AssetRangeReader` (io depends on core; a type import), inflating each block with core's `inflateStream(…, 'gzip', buf)`
(the platform's decoder checks CRC-32 and ISIZE; core's JS inflater where the platform has none). Nothing Seq\*-shaped:
lines and offsets, not tables, so a Mol\* reader (a bgzipped trajectory index, a tabix-indexed annotation) uses them too.

```ts
// @molstar/io/common/bgzf
/** A block's file offset × 2^16 + an offset into its inflated bytes; a number (files under 2^37 bytes: 128 GiB). */
export type VirtualOffset = number;
export const VirtualOffset: {
  of(block: number, within: number): VirtualOffset;
  block(v: VirtualOffset): number;
  within(v: VirtualOffset): number;
  /** A little-endian uint64 at `at` (a `.tbi`'s), refused past 2^53. */
  read(view: DataView, at: number): VirtualOffset;
};
export interface BgzfLimits {
  /** Inflated bytes one read gives, its blocks together (default 64 MiB: the store's answer bound). */
  readonly maxInflatedBytes: number;
  /** Inflated blocks a reader keeps for a pan back (default 16). */
  readonly cacheBlocks: number;
}
export const DefaultBgzfLimits: BgzfLimits;
/** Not BGZF, a block's `BC` field, size, CRC-32 or ISIZE wrong, or a bound passed; by `name` across copies. */
export class BgzfError extends Error { readonly kind: 'malformed' | 'too-large' }
export function isBgzf(bytes: Uint8Array): boolean;                    // gzip, FEXTRA, a `BC` subfield
/** A BGZF file read whole (a `.tbi`, an unindexed `.vcf.gz`): every block inflated in turn, joined. */
export function inflateBgzf(bytes: Uint8Array, limits?: Partial<BgzfLimits>, signal?: AbortSignal): Promise<Uint8Array>;
export interface BgzfReader {
  /** Inflated bytes from `begin` to `end` (the end block read for its first `within(end)` bytes). */
  read(begin: VirtualOffset, end: VirtualOffset, signal?: AbortSignal): Promise<Uint8Array>;
  /** Whether the file ends with the 28-byte EOF block (one read of its last bytes; false: perhaps cut short). */
  hasEof(signal?: AbortSignal): Promise<boolean>;
  dispose(): void;
}
export function bgzfReader(ranges: AssetRangeReader, limits?: Partial<BgzfLimits>): BgzfReader;
```

```ts
// @molstar/io/common/tabix
export interface TabixIndex {
  readonly kind: 'tbi' | 'csi';
  readonly format: 'generic' | 'sam' | 'vcf';
  readonly zeroBased: boolean;                 // the 0x10000 flag (BED)
  readonly columns: { readonly seq: number; readonly begin: number; readonly end: number }; // 1-based; end 0: none
  readonly meta: string;                       // a header line's first character (`#`)
  readonly skip: number;                       // lines skipped at the top
  readonly names: readonly string[];           // the sequences, in index order
  readonly minShift: number;                   // 14 for a .tbi
  readonly depth: number;                      // 5 for a .tbi
}
export interface TabixLimits extends BgzfLimits {
  readonly maxIndexBytes: number;              // an index inflated (64 MiB)
  readonly maxRegionBytes: number;             // a region's chunks, compressed, together (8 MiB)
  readonly maxLineBytes: number;               // one line (16 MiB: a VCF line of 100,000 samples at 160 bytes each)
  readonly maxHeaderBytes: number;             // the header lines (16 MiB)
}
export const DefaultTabixLimits: TabixLimits;
export class TabixError extends Error { readonly kind: 'malformed' | 'too-large' }
/** A `.tbi` or `.csi` (BGZF, told by its magic), every count and offset checked against `limits` and its own bytes. */
export function readTabixIndex(bytes: Uint8Array, limits?: Partial<TabixLimits>, signal?: AbortSignal): Promise<TabixIndex>;
export interface Chunk { readonly begin: VirtualOffset; readonly end: VirtualOffset }
/** The chunks that may hold lines meeting `[start, end)` (0-based, half-open): reg2bins, the linear index, merged. */
export function regionChunks(index: TabixIndex, name: string, start: number, end: number): readonly Chunk[];
/** The compressed bytes of the chunks `regionChunks` gives, without reading them (a level's cost, VCF below). */
export function regionBytes(index: TabixIndex, name: string, start: number, end: number): number;
/** Lines meeting a region, as bytes: each line's span in `bytes` and its extent by the index's columns (tbx_parse1). */
export interface TabixLines {
  readonly bytes: Uint8Array;
  readonly lineStart: Uint32Array; readonly lineEnd: Uint32Array;
  readonly begin: Int32Array; readonly end: Int32Array;   // 0-based, half-open
}
export function readRegion(reader: BgzfReader, index: TabixIndex, name: string, start: number, end: number,
  limits?: Partial<TabixLimits>, signal?: AbortSignal): Promise<TabixLines>;
/** The header lines (those starting with `meta`, and `skip` lines) at the top of the file. */
export function readHeaderLines(reader: BgzfReader, index: TabixIndex, limits?: Partial<TabixLimits>,
  signal?: AbortSignal): Promise<TabixLines>;
```

- **No `.gzi` reader**: a `.gzi` (bgzipped FASTA by uncompressed offset) serves only the designed-only FASTA runtime,
  so `parseGzi`, `GziIndex` and `readUncompressed` are not in io's surface (appendix). Why: a public name is a promise.
- **Internal to io**: the block walk (`BSIZE`) and `reg2bins` are not exported; tests reach them through the module's
  `_test` folder. Why: a public name is a promise, and no caller outside the readers needs them.
- **Lines as bytes**, not strings: io's rule (text is tokenized as bytes, io's spec, Decisions), and a region's lines
  are read once into the format's columns. **A line's extent is the index's**: `tbx_parse1`'s rule per format (generic:
  `begin`, `end` or `begin + 1`; VCF: `POS − 1` to `INFO/END` when the line's INFO has `END=` (whatever its alleles),
  else to `POS − 1 + len(REF)`; `0x10000` 0-based), so
  `readRegion` returns exactly what the index would; a line outside the region from a chunk's edge is dropped there.
- **CSI with TBI**: a `.csi` has no linear index (a bin's `loffset` instead) and any `min_shift` and `depth`; one reader
  over both. Why: a `.tbi` cannot index past 2^29 (512 Mb), and plant chromosomes pass it. **Positions stay int32**
  (the protocol's): a sequence longer than 2^31 − 1 (an axolotl or lungfish chromosome) is refused at open, said
  (`seq.source.too_large`, naming it), the others read. Why: widening every position column for a handful of genomes
  is the protocol's change, not a format's.
- **The EOF block**: `hasEof` is read once at open; a missing one is a warning (`seq.source.bgzf.eof`: "the file may be
  cut short"), the file still read (htslib's rule).
- **Why io, and why not core's `util/zip`** (base-library.md placed BGZF in core's zip "when one is built"): Gate 0, 2
  put it beside io's readers; it reads files (offsets, indexes, lines), which is io's, and core's zip is a codec. It
  reuses core's `inflateStream` and `ungzip`, never a second inflater.
- **What it amends elsewhere**: base-library.md's BGZF row moves to `@molstar/io/common/bgzf` and gains tabix; io's
  `_spec` gains `bgzf.md` (BGZF, TBI, CSI: their rules and bounds); `seqx-core`'s restricted imports
  (`biome.json`) gain the two entries beside the text CIF parser: an architecture rule widened by the gate's decision,
  not to make an import pass. A v6 changelog entry (an addition, no break).
- **Core's `ungzip` fails on a multi-member stream** (verified): `gzipDeclaredSize` reads the last
  member's ISIZE (BGZF's EOF block: 0), so the platform path throws on its first chunk; the platform decoders refuse
  concatenated members anyway (Node: "trailing junk"); the fallback `ungzipBlocks` goes through `_inflate`'s raw path
  and throws too. Only where `DecompressionStream` is missing does it come back short (the first member). **The fix
  splits members, never resizes a buffer**: each member found (by its `BC` field's BSIZE when present, else where the
  JS inflater's deflate stream of the member ends), inflated as one member with its own ISIZE, the outputs joined within
  the caller's bound. BZ lands the test first (BGZF, two plain members, a member cut short), then the fix. The genome
  runtimes read BGZF raw (no `decompress`) through io either way. Decided, 2.

### What every runtime shares

- **Bytes to text** (`internal/text.ts`): a file is plain, gzip (one member, `DecompressionStream`) or BGZF (by `isBgzf`;
  `inflateBgzf`), told by its magic, never its extension; bounded by `limits.maxBytes` before and after inflating (as
  mmCIF). UTF-8 (the records' rule; U+FFFD said once).
- **Whole or by region**: a tab-line format (`bed`, `bedgraph`, `gff3`, `gtf`, `vcf`) with an `index` reads by region
  through tabix (`internal/lines.ts`); without one it is read whole at open within `limits` and answered by region by
  the `seq` runtime's join (`openSeq`, as records). One line parser either way. **Why one parser**: the conformance
  property "a region read equals the whole read clipped" is then a property of the line source alone. An index that
  names no line the file has, or a `.tbi` for another file, shows as lines outside their chunk's region: `malformed`
  for that region, never rows in the wrong place.
- **Levels of a tab-line source by region, per source, from its index's bytes**: candidate levels `1, 16, 256, 4096`
  (up to the runtime option `maxLevel`, 4,096: 4 Mb tiles), each a tile width for items (RR's items read at a tile
  size, extended to `ranges` and `matrix`: RR's), and above 1 a `density` (ranges) or `bins`
  (matrix) **counted exactly from the records of that tile**. At open the runtime lists a level only when **every tile
  of it on every system fits `maxRegionBytes`** (`regionBytes`: the chunks' compressed bytes, read from the index,
  never from the data). So ClinVar (about 30 bytes a kb) lists all four; a 1000 Genomes chromosome (tens of MB a Mb,
  mostly genotypes) lists `1` and `16` (16 kb tiles, views up to about 4 Mb under `maxTilesPerWant`), and a wider view
  is the layer's built too-many-tiles state ("zoom in"). Why per source and from bytes: one runtime reads sparse and
  dense files, and the bytes a tile costs are known exactly before a read, where its rows are not. **Why counts exact
  or absent**: a density estimated from the index's bytes (the current Seq\*'s `estimate`, JBrowse's gate) would draw a
  guess as data; the listed coarsest level is igv's visibility window, made per file. A dense file's sites wider than
  that are a host's bigWig of their density (built; Tensions).
- **A view's tiles together are bounded too**: the source gives each tile's cost before reading it
  (`readCost(request)`: `regionBytes`, the chunks' compressed bytes from the index), and the store sums a layer's tiles
  within `maxBytesPerLayer` (16 MiB, RR); past it the layer reads coarser, then shows "zoom in". So ClinVar at the fit,
  every tile of whose four levels fits `maxRegionBytes` but whose tiles together are its whole file (about 90 MB), is
  "zoom in" until the view is a few hundred Mb wide. Why: a bound per tile does not bound a view. `samples` narrows the states held, not the bytes read (genotypes are
  inside the compressed lines); a tile's states stay within RR's `maxTileBytes`.
- **Issues** (`seq.source.<format>.*`, one a cause with a count and the first ten places, EC-86; severities EC-87): a
  line that does not parse is left out (`.line`, warning, with the line number or virtual offset); a structure that
  cannot be read refuses the source (`seq.source.malformed`); a bound passed is `too-large` for a region and
  `seq.source.too_large` at open, hinting the bound and the remedy (an index, a coarser level, a smaller file).
- **Range reads**: `AssetReader.openRange` (built); a host without it reads the file whole within `limits.maxBytes`
  (bigWig's rule). A failed range read (`AssetRangeError`: the file changed under the reader) fails the region
  (`failed`, EC-174), never wrong rows.
- **A format's conventions are fixed** (`bed`, `bedgraph`, `bedpe`, `vcf`, `gff3`, `gtf`): the
  source's built `base` and `ends` (files' defaults `1` and `inclusive`) are not read for them, said once if written
  (`seq.source.<format>.transform`, info), as the GFF3 preset already does. Why: a writer setting `base: 1` on a BED
  would move every row by one (EC-195).
- **Names are the file's**: a source's systems are the names its file writes (`1` or `chr1`); GA's names rule maps
  them onto an assembly (EC-188, EC-189). A system's length is the file's where it says (`##contig`,
  `##sequence-region`, a chrom.sizes), else absent: the assembly's stands.

### BED, bedGraph, BEDPE, chrom.sizes (`runtimes/bed`)

- **BED**: tab-separated (a line without a tab split on whitespace, said once), `track`, `browser` and `#` lines read
  past; base 0, ends exclusive, fixed. Columns past the third are BED's by place (`name`, `score`, `strand`,
  `thick_start`, `thick_end`, `item_rgb`, `block_count`, `block_sizes`, `block_starts`), typed by the format (`score`
  float32, `thick_*` int32, `strand` the role). **A `bed N+M` file names its own through the built `fields`, a column
  by its 0-based number** (the records' path rule: an all-digit segment is an index; a headerless line is an array of
  cells): `fields: { stain: '4' }` for UCSC's cytobands (the navigator's `bands`), `{ signal: '6', p: '7' }` for
  `narrowPeak`. A column a mapping reads is not read under BED's name (the records' rule), so `gneg` is never parsed
  as a score; columns past twelve that no mapping names are `col<k>` category fields. Types by `types` (built). Why no
  new param: `fields` already says which column is which (EC-190), and files.md's refusal of `names` stands.
- **BED12's blocks are pieces** of one item: a row a block (`[start + blockStart, + blockSize)`), the item's whole
  extent on every piece (`itemStart`, `itemEnd`, below), `thick_start` and `thick_end` fields; blocks sorted, inside the
  item, counts equal to `block_count`, else the line left out (`.line`). **Identity** (EC-196): `name@start-end` (BED
  names repeat), `id` the category; the name is the label.
- **bedGraph**: `chrom start end value`, served as `values` (sparse positions, NaN between); coarse levels the store's
  summaries (read whole) or bins counted from the tile (by region), `mean`, `min`, `max`, `count` as `ValuesBins`.
- **BEDPE**: `chrom1 start1 end1 chrom2 start2 end2 name score strand1 strand2` and `fields` by number past them, as RR's
  `pairs` list; **read whole** within `limits`, an `index` ignored and said (a `.tbi` finds a pair by its first end
  only), then **indexed by both ends** once: each row's position in two orders, by its first end and by its second,
  sorted by system and start (two `Int32Array` permutations), so a region read finds every pair with either end in it
  (RR's contract). Its rows, both orders and RR's pairs summaries count together against `summaries.maxBytes` (16 MB,
  about 600,000 pairs); past it `too-large`, hinting pairix or a levelled pairs format (the Hi-C route, later). No
  doubled BEDPE (each line twice under tabix) is read: 12,000 links is the milestone. `id` the name unless `.` (RR:
  `.` is no identity), else `row`.
- **chrom.sizes**: `name length` lines; a coverage of systems and lengths, no shapes. Why here: an assembly's lengths
  come inline or from a chrom.sizes (GA).

### GFF3 and GTF: transcripts (`runtimes/gff`; files.md's tension closed)

- **The hierarchy, assembled by the line's type** (`assembleTranscripts(rows)`, exported: a pure function over parsed
  rows, tested alone), so a tile sees the same thing as the whole file: a line is a **transcript** when its `type` is
  one of `transcriptTypes` (a runtime option; default the Sequence Ontology's transcript types written in Ensembl,
  GENCODE, NCBI and Prokka files: `mRNA`, `transcript`, `primary_transcript`, `lnc_RNA`, `ncRNA`, `rRNA`, `tRNA`,
  `snRNA`, `snoRNA`, `miRNA`, `scRNA`, `misc_RNA`, `pseudogenic_transcript`, `unconfirmed_transcript`,
  `V_gene_segment`, `D_gene_segment`, `J_gene_segment`, `C_gene_segment`, and any type ending in `_transcript`), a
  **piece** when its type is `exon`, `CDS`, `five_prime_UTR`, `three_prime_UTR` or `UTR` (GTF adds `start_codon`,
  `stop_codon`), a **gene** when a transcript names it as `Parent`, else of its own. **A piece whose parent is not a
  transcript type** (a prokaryote's `gene` → `CDS`, NCBI's) makes that parent a transcript-level item too: exact by
  region only where its pieces cover it (a prokaryote's CDS spans its gene), and said once by type when read whole
  (`seq.source.gff.type`, info: "`X` read as a transcript; add it to `transcriptTypes`"). Never by children seen in a
  tile: a tile inside an intron holds the gene and transcript lines and no piece (EC-227, EC-229).
- **Each transcript an item** (`id` its `ID`): its pieces its `exon` lines, else its CDS and UTR lines merged;
  `thick_start`, `thick_end` the hull of its CDS (in GTF with `start_codon` and `stop_codon`, which GTF writes outside
  the CDS), `cds_phase` the phase of its first CDS piece in the transcript's direction (what a codon map later needs;
  `rows` keeps every line's `phase` as an int8 field); its strand its own; its fields its line's attributes, with
  `gene` (the gene's `ID`), `gene_name` (the gene's `Name`, GTF's `gene_name`) and `type`; **its extent its own
  line's** (`itemStart`, `itemEnd` on every row). A child naming several parents (`Parent=tx1,tx2`: a shared exon) is a
  piece of each. Every other feature is an item of its own by `ID` (rows sharing an `ID` its pieces, as built).
- **A span row**: a transcript whose pieces all fall outside a tile (the tile inside an intron, or between its UTRs)
  gives one row in that tile, `start = end` at the tile's clipped start of the item, with `itemStart`, `itemEnd`: the
  item, no piece here (protocol, below). Readers draw its hairline across the tile and hit it as the item. **The whole
  read clipped to a tile gives the same row** (`sliceTable` emits a span row for an item with pieces elsewhere), so
  region and whole agree by construction, and the property test holds them to it.
- **GTF**: attributes `key "value";` (unquoted numbers too; a key repeated, GENCODE's `tag`, joined by `,`); the
  hierarchy by `transcript_id` and `gene_id`. **With `transcript` lines** (GENCODE's and Ensembl's GTFs) the extent is
  exact, as GFF3's. **Without them** (GTF 2.2): read whole, a transcript's extent is the hull of its rows, exact; by
  region, a tile only knows the rows it read, so `itemStart`, `itemEnd` are that tile's hull and readers take the union
  over the tiles held (it grows as tiles arrive), a tile inside an intron draws nothing for it, and the source says so
  once (`seq.source.gtf.extent`, warning: "no transcript lines: read it whole, or extents are partial by region").
- **`features`** (new source param): `transcripts` (the default), `genes` (one item a gene: its transcripts' exons
  merged, its CDS hull, the gene's attributes), `rows` (the file's lines, as the records preset). Why a param: IGV's
  expanded and collapsed are two figures over one file, and a patch switches them (one undo step, the source reopened
  in place: files.md's key rule).
- **`##sequence-region`** gives a system's length (EC-188), and a `region`, `chromosome` or `supercontig` feature
  covering a whole sequence describes it, not an item (said once): Ensembl writes one a chromosome, which would be a
  bar under every view.
- **By region** (a bgzipped, sorted GFF with a tabix `index`): a tile's lines are assembled alone, by type; a
  transcript's line spans its children, so every tile meeting a transcript holds its line (and its gene's), and the
  pieces of one item join across tiles by `id` (protocol, Tables). A child outside its parent's extent (malformed GFF3)
  keeps its parent's `id` and lacks its fields in that tile, said once (`seq.source.gff.extent`).
- **`##FASTA`**: everything after it is ignored (the records' rule; Prokka and NCBI append the sequences), its bytes
  still counted against `limits.maxBytes`; letters are a FASTA source's.
- **Bounds**: a `Parent` chain deeper than 16 or a cycle leaves those rows as items of their own (`.parent`, warning); a
  row's attributes at most 256; read whole within `limits`.
- **Why transcripts by default**: EC-227 (gene, transcript, exon) and every precedent draw a gene model; a flat GFF3 of
  exons is rows with no reading. **Why not dynamic transcript nodes** with codon maps (dynamic-nodes.md's p53
  transcripts): a protein's coordinate system from a GFF's CDS is a figure no Genomes view needs yet (Decided, 7).
- **The picture**: the ranges layer's `line: 'gene'` (new `LineStyles` member, GG): pieces as boxes, thick between
  `thick_start` and `thick_end`, thin outside, a hairline between pieces across the item's extent with strand
  chevrons. Why a line style: the ranges layer already lays out items with pieces; a gene model is how they are joined.

### VCF (`runtimes/vcf`)

- **Header** (`readVcfHeader`, within `maxHeaderBytes`): `##contig` (names, lengths: EC-188), `##INFO` and `##FORMAT`
  (`Number`, `Type`), the sample names. The fixed columns and INFO are parsed per line; **genotypes only when a `matrix`
  is asked** (most reads are sites).
- **A row per ALT allele** (EC-196), **every allele of a line with the line's extent, the index's**: `POS − 1` to
  `INFO/END` when the line has `END=` (whatever its alleles), else to `POS − 1 + len(REF)`; the anchor base kept. Why
  one rule equal to `tbx_parse1`'s: a region read then agrees with the index at its edges (igv trims the anchor; a
  Tension). **`id`** `pos:ref>alt` (the VCF `ID` repeats over alleles), with `ref` and `alt` replaced by
  `<length>:<FNV-1a hex>` when longer than 32 characters (a sequence-resolved SV's), and `#2`, `#3` … on a second and
  later identical key in file order (a duplicated record: always read together, being at one `POS`, so the suffix is
  stable across tiles), and `@<end>` when `INFO/END` extends the record past its REF (two records differing only in `END` are two items, whole and by region; VC's review); the `ID` a field `name`, the label. Fields: `ref`, `alt`, `qual` (float32), `filter` (category), `category` (`SNV`, `MNV`,
  `insertion`, `deletion`, `indel`, `structural` for `<…>`, `breakend` for `]`/`[` alleles), and each declared INFO key
  under its name (`Integer` int32, `Float` float32, `Flag` and `String` categories; `Number=A` split by allele,
  `Number=R` the reference's and the allele's, other numbers their text), up to `maxInfoFields` (256): the keys the
  source's `fields` or `types` name first, then the header's order (gnomAD v4 declares about 900), the rest said once
  with their count. A key clashing
  with a fixed field is said once (`.cell`, the records' rule). ALT `.`, `*` and `<NON_REF>` are no item (counted, info).
- **The matrix** (RR's `matrix` shape, `MatrixItems`): items the ALT alleles in the tile × states the samples
  (`ShapeInfo.states`, in the header's order, narrowed by **`samples`**, a new source param: a list of names);
  `codes` the **genotype class** of that allele in `GT`, whatever the ploidy: `0` no copy, `1` some but not every copy
  (`0/1`, `1/2` for either allele), `2` every copy (`1/1`, and a haploid `1`: chrY, chrM, a male's chrX), `3` partly
  called (`./1`, `./0`: a `.` beside called alleles), `255` uncalled (`.`, `./.`). Why classes and not a dosage count:
  a haploid `1` and a heterozygote's `1` are different genotypes a theme must not paint alike; the exact copy count is
  `value` float32 by `RegionRequest.field: 'GT'` (copies of the allele, NaN uncalled), as any FORMAT field (`DS`, `GQ`;
  `Number=1` or `A`). Coarse levels: `bins` counted
  from the tile (majority and share, `mean`), as ranges. Bounded by RR's `maxTileBytes`. A long deletion is one item
  across the SNPs it covers (EC-222).
- **Without an index**: a plain or bgzipped VCF read whole within `limits`.
- **Dense files** (1000 Genomes, a biobank): the levels the index's bytes allow and the view's read cost (above), so a
  dense chromosome is read in 16 kb tiles of a few hundred kB; a wider view is "zoom in", never an 80 MB read. Genotypes are parsed only for a
  `matrix` request, and only for the `samples` held; the line's bytes are read either way (they are in the compressed
  block), which the level rule already counts.
- **One file a chromosome** (1000 Genomes, gnomAD) is one source a file this slice; a source over several files is
  open (Tensions).

### The `seq` runtime: sided pairs found from either end, a map's runs

The dataset's synteny, dot plot and liftover views read their alignment from the `seq` format, written from the
alignment truth in code (GD); PAF and chain runtimes are designed only (appendix).

- **Side `a` is the side the alignment maps from** (RR's sided sources): a map read from the source is `from` side `a`
  `to` side `b` with no per-format flag; a host's PAF adapter puts the query on `a`, a chain adapter its `t`.
- **Sided pairs, found from either end**: a pairs table in the `seq` format may mark each system's side
  (`SystemInfo.side`, RR); the runtime reads it whole and indexes it by both ends as BEDPE (above), each side's
  systems apart, within RR's `summaries.maxBytes`, and answers a request by its `side` and region from either end.
  Rows sharing `row` are the blocks of one alignment (its pieces: a hit names the alignment). Fields as written
  (`matches`, `block_length`, `score`).
- **As a map** (`map.source`, new language param): a `map` with `method: segments` may name a pairs source instead of
  inline `segments`; its runs are the source's rows of equal-length ends, read by region of either side
  (`RegionRequest.side`): the `from` side for projection forward, the `to` side for a layer drawn on the `to` section
  (a liftover's old genes on the new assembly); the runtime indexes both ends once, so it is one read
  (building-blocks.md: "a genome-scale map is held by region"; GA: one map over every sequence pair). A row of unequal
  ends is refused as a run (counted, `seq.source.seq.run`). Why the pairs shape: a run is a pair of equal length, so one
  source serves the ribbons (LX) and the liftover (EC-212) from one read.

### A name index (EC-228)

- **A capability, not a shape**: `OpenedSource.findNames?(text, { max, prefix }, signal): Promise<NameHit[]>`, with
  `Coverage.names` (the fields searched) and `NameHit { system, start, end, id?, label, field }`. Case-insensitive,
  exact first, then by prefix when asked; at most `max` (20) hits.
- **Who serves it**: a GFF3 or GTF read whole (GFF3's `Name`, `ID` and `Alias`; GTF's `gene_name`, `gene_id`,
  `transcript_name`, `transcript_id`), from a map built once at open (the author, 2026-11-21). BED and VCF names are
  not indexed this slice, nor a tabix source (a text index such as trix is a later format), so
  a GFF3 past `limits.maxRows` (a human GENCODE GFF3 is about 3M lines) finds nothing by name (Tensions).
- **Who calls it**: `findNames(sources, text, options, signal)` (`data/names.ts`, a function over the opened sources
  a viewer's source sync holds: not in `store.ts`, so no store lane is touched), merging hits in source order; GA's
  `(locus :text …)` falls back to it when no sequence matches (GC3's: an asynchronous focus, cancelled by the next),
  and a hit focuses the item's address (`(item :id ENST…)` on its source) or, without an `id`, its range.

### Protocol additions (`packages/seq`, public API, recorded in the lanes' changelog entries)

Each with its owning lane (below, Lanes):

- `RangesList.itemStart?`, `itemEnd?` (`Int32Array`): each row's item's whole extent, when rows are pieces (BED12,
  transcripts); and **span rows** (`start = end` with an extent: the item, no piece in this tile), which `sliceTable`
  emits for an item whose pieces lie outside the slice. Why: a reader of region data knows an item's ends from any
  tile (today `extents` exist for inline data only, so region items get no arrowhead), and an intron-only tile still
  draws its transcript (EC-229: flanks of streamed genes need the extent). RR's text; RG (tables and `checkAnswer`).
- `SystemInfo.side?`, `RegionRequest.side?`, the key's sided system part (RR's). RG.
- RR's `RegionRequest.items` for `ranges` and `matrix` lists too (keys change: RR's `i<level>` part). RG.
- `OpenedSource.readCost?(request)` and the store's `maxBytesPerLayer` (RR's). RG; the tabix sources give it in FR.
- `OpenedSource.findNames?`, `Coverage.names?`, `NameHit` (types in `protocol.ts`, no store change) and
  `data/names.ts`. BZ (types only, before RG; RG merges over them).
- Field names fixed by convention and documented in the Data spec: `thick_start`, `thick_end`, `cds_phase`, `phase`,
  `item_rgb`, `gene`, `gene_name`, `category`, `name`.
- Language (GC1's bump, the only lane on the language in stage 1): a source's `features`, `samples`; `map.source`; the
  ranges layer's `line: 'gene'`. Declared, documented and checked by GC1; read by the format lanes after GC1 merges.

### Bounds (security.md: a file is untrusted input)

Each refused before it allocates; `malformed` for a structure, `too-large` for a region (a coarser level, a smaller
region), `seq.source.too_large` at open.

| Format | Bound | Default | Over it |
|---|---|---|---|
| BGZF | a block's compressed and inflated size; its ISIZE and CRC-32 | 65,536; the platform's check | malformed |
| BGZF | a read's inflated bytes, its blocks together | `maxInflatedBytes` 64 MiB | too-large |
| tabix | the index inflated; its references; a reference's bins (TBI: 37,450 and the metadata pseudo-bin), chunks, a chunk `begin ≤ end` and inside the file | `maxIndexBytes` 64 MiB; `limits.maxRows` | malformed (the source refused) |
| tabix | a region's chunks, compressed, together (a level not listed when any of its tiles passes it); one line; the lines of an answer | `maxRegionBytes` 8 MiB; `maxLineBytes` 16 MiB; `limits.maxRows` | too-large; malformed; too-large |
| tabix | a layer's tiles in view, their chunks' compressed bytes together (`readCost`) | RR's `maxBytesPerLayer` 16 MiB | a coarser level, else "zoom in" |
| all | a sequence's length | 2^31 − 1 (int32 positions) | that sequence refused at open, said |
| VCF | the header; samples; ALT alleles a record; INFO fields | `maxHeaderBytes` 16 MiB; 100,000 (a line of them within `maxLineBytes`); 256; 256 | malformed; malformed; the record left out; the rest not fields (said) |
| GFF, BED, VCF read whole | bytes (after inflating), lines | `limits.maxBytes` 64 MB, `limits.maxRows` 1,000,000 | `seq.source.too_large`, hinting an index |
| BEDPE, the `seq` format's pairs | rows held with both orders and their pyramid | RR's `summaries.maxBytes` 16 MB (about 600,000 pairs) | `seq.source.too_large`, hinting pairix or a levelled pairs format |
| GFF | a `Parent` chain's depth; attributes a row | 16; 256 | rows kept flat (said); the rest dropped (said) |

## Edge cases this design takes or assigns

- **Taken here**: EC-186 (an index named, never `<uri>.tbi`), EC-188 (lengths from `##contig`, `##sequence-region`,
  chrom.sizes), EC-190 (`fields` by name or column number), EC-195 (a format's base fixed, converted at the edge,
  EC-1), EC-196 (ids `name@start-end`, `pos:ref>alt`, a row a BED12 block and an ALT allele), EC-206 (strand on every
  pairs row: BEDPE's, the `seq` format's), EC-208 (pairs read by region from either end: BEDPE and the `seq` format's
  pairs read whole and indexed by both ends within a byte budget), EC-220 and EC-222 (the VCF's matrix: missing as 255,
  a long deletion one item), EC-227 and EC-229 (gene, transcript, exon; a span row in an intron-only tile), EC-40 and
  EC-92 (every reader fuzzed without a throw), EC-51 (a row past a known length clipped, counted), EC-86 and EC-87
  (issues), EC-169 and EC-174 (pending, failed and `AssetRangeError`), EC-171 (its outcome for tabix: exact counts
  within bytes, a view's read cost bounded, else "zoom in").
- **Its data part here, the rest assigned**: EC-228 (`findNames` over a GFF3 or GTF here; the box and the focus GA's locus
  fallback, GC3), EC-189 (a file's names here; aliases GA's names rule), EC-207 (a pairs source as one map's runs,
  `map.source`, here; paths GA's, GC1), EC-212 (the runs read by either side here; the projection GA's, the figure
  LX's), EC-209 and EC-210 (data here; ribbons LX's).
- **Assigned**: EC-213 (MAF: cut), EC-223 (`.hic`, `.mcool`: not in the formats, RR's Decided 10), EC-224 (keyed results:
  RR's `rows`; the join deferred to the close-out slice), EC-216 and EC-219 (GWAS points and LD: the Manhattan family,
  deferred to the close-out slice).

## Defaults, with a picture

What a reader sees with a source of the format and one layer over it, nothing else set.

1. **GFF3 or GTF** (`features: transcripts`, `line: gene`): a transcript a lane, boxes for exons, thick CDS, thin UTRs,
   a hairline with chevrons between, labelled by `Name` (else `ID`).
   ```
   TP53-201 ◂━━┳━━━━━━▮▮▮━━━▮▮━━━━━━━▮▮▮▮━━━━━━━▮▮▮▮▮▬▬   (thin UTR ▬, thick CDS ▮, intron ━ with ◂)
   TP53-202 ◂━━━━━━━━━▮▮▮━━━▮▮━━━━━━━━━━━━━━━━━━━▮▮▮▮▮▬▬
   ```
   `features: genes`: one lane a gene, its exons merged.
2. **VCF zoomed in**: a variant a mark at its REF extent, coloured by `category` (the next categorical colours);
   a long deletion a bar. With a matrix layer, RR's picture 6 (a row a sample).
   ```
   ▏ ▌  ▌▌    ███████   ▌  ▏   SNV ▌, deletion ███
   ```
3. **VCF past its coarsest level** (`[1, 16, 256]`: a view wider than about 67 Mb) **or past the view's read cost**
   (ClinVar wider than a few hundred Mb): the layer's too-many-tiles state (built: the band hatched, "zoom in to read
   this track"), never an estimate.
4. **BED with `fields: { stain: '4' }`** under a `bands` layer: the ideogram (navigators.md's picture).
5. **Pending**: as built (the band's pending hatch); **failed** (a range read refused, a changed file): the placeholder
   with its issue; **malformed** region: the placeholder, the issue naming the virtual offset or line.
6. **Dark**: tokens only; nothing a format adds is a colour.

## The interaction matrix's cells

Formats add no gesture of their own: what a file holds is drawn by the layers that read it, as files.md says. The cells
they change:

| | select | extend | add / remove | focus | back | zoom | pan | hide | edit |
|---|---|---|---|---|---|---|---|---|---|
| **an element: a transcript** | click: the transcript (every piece, `(item :id …)` on its source) | Shift+click: the hull to another | Cmd/Ctrl+click | double-click: its extent (`itemStart`, `itemEnd`) | Alt+Left | — | — | no: a filter | no: data |
| **an element: a variant allele** | click: the allele's REF extent | Shift+click | Cmd/Ctrl+click | double-click: its extent | Alt+Left | — | — | no: a filter on `category`, `filter`, INFO | no |
| **an element: an alignment's block** (the `seq` format's pairs) | click: the block as drawn (blocks merged under a pixel are one), both ends lit | Shift+click | Cmd/Ctrl+click | double-click: the alignment's hull on each side | Alt+Left | — | — | no | no |
| **a track or layer** | as built | as built | `track.add` as built | as built | as built | — | — | `layer.hidden`; the data filter | the source's `features`, `samples`, `fields`: a patch, one undo step, the source reopened in place (files.md) |
| **a bookmark or link** | — | — | as built | an item's address (`(item :id ENST…)`) waits for its region (dynamic-nodes' rule) and survives a `features` switch only for `transcripts` ↔ `rows` ids it keeps (a Tension) | as built | — | — | — | — |
| **the view** | — | — | — | a name typed (a host's box, the locus text) finds through `findNames`, then focuses (GC3) | Alt+Left | — | — | — | — |

## Configuration

| Item | Whose | Default | Why the default; the alternatives |
|---|---|---|---|
| each runtime | registry entry | not listed | opt in: a viewer pays for what it reads; listed by a host |
| `gff3` by the genome runtime | `withAdapter` | the default (flat) | the built behavior; the genome runtime by composition |
| `source.index` (built) | spec param | none | where the index is (EC-186); `.tbi`, `.csi` |
| `source.fields` by column number (built param, new reading) | spec param | BED's names by place | a `bed N+M` file names its own columns |
| `source.features` | spec param | `transcripts` | a gene model (EC-227); `genes`, `rows` |
| `source.samples` | spec param | every sample | a figure's subset; a list of names |
| `source.fields`, `types` (built) | spec params | the format's names and types | renames (EC-190) |
| `map.source` | spec param | — (inline `segments`) | a pairs source's equal-length runs, read by region of either side |
| `maxLevel` of a tab-line source by region | runtime option | 4,096 | the coarsest candidate (4 Mb tiles); the levels listed are those the file's index says fit `maxRegionBytes` |
| `maxBytesPerLayer` (RR's) | `StoreSettings` | 16 MiB of read cost | a view's tiles together, from `readCost`; past it coarser, then "zoom in" |
| `transcriptTypes` (GFF, GTF) | runtime option | the SO transcript types above, and `*_transcript` | a line's role by its type, so a tile reads as the whole; a list |
| `maxRegionBytes` | runtime option (`TabixLimits`) | 8 MiB compressed | a tile's read stays a few hundred ms remote; any |
| `maxInflatedBytes` | runtime option | 64 MiB | the store's answer bound |
| `maxIndexBytes` | runtime option | 64 MiB | a dbSNP `.tbi` fits; any |
| `maxLineBytes`, `maxHeaderBytes` | runtime option | 16 MiB; 16 MiB | a VCF line of 100,000 samples at about 160 bytes each (`GT:AD:DP:GQ:PL` of a multi-allelic site); gnomAD's header a few hundred kB |
| `maxInfoFields`, `maxSamples`, `maxAlleles` (VCF) | runtime option | 256; 100,000; 256 | bounded columns and states |
| `cacheBlocks` (BGZF) | runtime option | 16 | a pan back without a read |
| `findNames`' `max`, `prefix` | option of the call | 20; false | a list a box shows; prefix for typing |
| a format's base and ends; side `a` the map's from side; the anchor base kept; a line's extent by the index's rule | fixed | — | the reasons in the Decision |

Not built: FR (the formats by region past level 1, VCF's matrix); BZ, BE, GG and VC built (GFF3 and GTF read whole; VCF sites whole and by region at level 1).

## Tensions

- **The anchor base**: a VCF indel spans its REF with the anchor (the index's extent), where igv draws the deletion
  without it and an insertion between bases. Recommendation: the index's extent (a region read agrees with the
  index); a `trimmed` drawing is the variants layer's later, if the use gate asks.
- **Phasing**: `0|1` and `1|0` give one dosage; haplotypes as two rows a sample need a states system twice the
  samples. Open until a view asks (a phased panel's figure).
- **No density past the coarsest level of a tabix source or its view's read cost**: ClinVar at the whole genome draws
  "zoom in". Accepted: honest and bounded; the whole-genome remedy for a 1D layer is a host's bigWig of its density
  (the built runtime), drawn beside or instead.
- **Breakends** (`N[chr2:321682[`) are one-position items with their mate in `alt`; served as `pairs` they would draw
  as links. Open until an SV view asks (a later lane, the pairs shape exists).
- **`itemRgb`** is a field, not a colour: a theme reading `r,g,b` text is the themes design's (a `token:value` form
  over a field). Open.
- **GTF without transcript lines by region** has partial extents (the tiles' hull) and draws nothing in an intron-only
  tile. Accepted, said; read whole for exact extents.
- **A type outside `transcriptTypes`** with exon children reads as a transcript only where its pieces cover it by
  region. Accepted, said once by type; the list is a runtime option.
- **A source over several files** (one VCF a chromosome: 1000 Genomes, gnomAD): one source a file this slice. Open: a
  `uri` template by system, when a view asks (GD's do not).
- **Genotype classes** lose the copy count past diploid (`1/1/1/0` and `0/0/0/1` are both `1`); the count is `field: 'GT'`.
  Accepted.
- **A bookmark of an item across `features`**: a gene's id in `genes` is not a transcript's; switching modes leaves a
  bookmark of the other mode waiting. Accepted: bookmarks of positions survive (files.md's rule).
- **Name search reads a GFF3 whole**: within `limits.maxRows` (1M lines); a human GENCODE GFF3 (about 3M) is
  `too-large` and finds nothing, and tabix sources have no name index. Open: a host's adapter serving `findNames`, or a
  text index (trix) later.
- **A GFF3 child outside its parent's extent** differs by region from the whole read (its fields absent in the tile).
  Accepted: malformed GFF3, said.

## Self-check

The claim: **this design covers what a user does with genome files.**

1. **Assume it holds.** A host composes the runtimes it reads; each opens its file by range, reads a region within
   bounds into the protocol's tables, transcripts and variants and genotypes and pairs; a name finds its
   item; the files' conventions are converted at the edge.
2. **What it rests on, negated.** (a) Files are indexed; negated: a plain VCF, a BEDPE, a GFF without `.tbi`. (b) Files
   are well formed; negated: crafted headers, stale indexes, a file replaced under the reader, a truncated BGZF. (c) A
   file is read by position; inverse: found by name. (d) Features are flat; negated: a three-level hierarchy, a shared
   exon, a GTF without transcript lines, a prokaryote's gene with a CDS. (e) One sample; plural: thousands. (f) One
   genome; plural: two sharing `chr1`. (g) Data is dense enough to read at the fit; negated: a sparse VCF over a
   genome, or one so dense a tile is megabytes. (h) Names match the assembly; negated: `1` against `chr1`, contigs the
   assembly lacks. (i) The platform inflates gzip as one stream; negated: BGZF's members, no `DecompressionStream`. (j)
   What is drawn is counted; negated: an estimate. (k) A format's conventions are one; negated: BED 0, GFF 1, VCF's
   REF extent. (l) A bounded tile is a bounded view; negated: every tile fits and the view is the whole file.
3. **The negated world.** (a) whole reads within bounds, a hint to index, BEDPE and the `seq` format's pairs indexed by
   both ends after the read; (b) bounds before allocation, EOF said, `AssetRangeError` failing the region; (c)
   `findNames`; (d) `assembleTranscripts` by type, a piece of each parent, span rows, hulls; (e) `samples`, genotypes
   parsed only when asked, RR's byte bound; (f) sided systems (the `seq` runtime); (g) levels listed per file from its
   index's bytes, exact counts, `maxRegionBytes`; (h) the file's names kept, GA's rule maps them; (i) io's block reader
   and core's `ungzip` fixed; (j) no estimate; (k) conventions fixed per format and converted once; (l) a view's read
   cost bounded (`readCost`, `maxBytesPerLayer`).
4. **Both held.** Integrated: every item of 3 in the Decision. Standing: phasing (e); breakends as pairs; name search
   past a GFF3 read whole (c); no density past a tabix source's coarsest level or its view's cost (g, l: a host's
   bigWig); the anchor base (k); GTF without transcript lines by region (d).

## Consequences

- **`@molstar/io`** (BZ): `common/bgzf.ts`, `common/tabix.ts` (the surface above), tests in `common/_test/` (blocks
  across a read's edge, the EOF block, TBI and CSI, `reg2bins` against the specification's table, a region's lines
  equal to a whole inflate filtered, `regionBytes` equal to the chunks' sum, every bound refused before it allocates,
  fuzzed blocks and indexes without a throw); `_spec/bgzf.md`; a v6 changelog entry.
- **`@molstar/core`** (BZ): `ungzip` over a multi-member stream (every member, within its bound), its test; a v6
  changelog entry (a fix).
- **`biome.json`** (the orchestrator's): `seqx-core`'s allowlist of `@molstar/io` gains `common/bgzf` and
  `common/tabix`.
- **`seqx-core`**:
  - `runtimes/bed`, `runtimes/gff`, `runtimes/vcf` and `runtimes/internal/` (text, lines); the package's exports
    `"./runtimes/internal/*": null` (hot file); the `seq` runtime's sided pairs, its index of both ends, `map.source`.
  - protocol (`data/protocol.ts`): `findNames`, `names`, `NameHit`, `data/names.ts` (BZ); RR's part (span rows and item
    extents, `side`, items for ranges and matrices, `readCost`, `checkAnswer`, `sliceTable`) in RG.
  - language (GC1's bump): `source.features`, `samples`; `map.source`; `line: 'gene'`; issues
    `seq.source.{bgzf,tabix,bed,gff,vcf,seq}.*`.
  - render: the ranges layer's `gene` line (GG), item extents from `itemStart`/`itemEnd` for region data.
- **`seqx-viewer`**: nothing but a host composing the runtimes (`withAdapter` for `gff3`); the find box is GC3's.
- **`seqx-examples`**: the test-only BGZF and tabix writer on core's `PCG`, which GD's generator uses for its bgzipped
  files; the recording script for the real slices (below); `sources.json` rows with licences.
- **Storybook**: `Features/Genome formats/<format>` (the file's head or its index's summary and the source node beside
  the canvas and the SVG, a misread file and its issues), each with its headless twin on recorded or generated files;
  a Design page "Genome formats".
- **Tests**, per runtime: the conformance suite (`data/_test/source-suite.ts`: coverage, region answers, empty and edge
  regions, cancellation) on a generated fixture and on its recorded real slice; **a region read equals the whole read
  clipped** at every level (a property over generated files and random regions, items deduped across tiles); fuzzed
  and truncated files read without a throw, every bound refused before it allocates; per format: BED12 pieces and
  `fields` by column number; BEDPE and the `seq` format's sided pairs found from either end and either side, deduped
  once, past the byte budget `too-large`; `assembleTranscripts` (three levels, a shared exon, CDS-only, GTF with and
  without transcript lines, a cycle, a prokaryote's gene, a type outside the list) and its region twin: **a tile inside
  an intron gives the transcript's span row**, every tiling of a generated GFF equal to the whole read clipped (EC-227,
  EC-229); VCF's extent equal to `tbx_parse1`'s per line (an `END=` on a SNV line too), genotype classes (haploid,
  partial, multi-allelic), ids of duplicates and long SVs, alleles, INFO types and `Number`s, a long deletion; levels
  listed from a dense and a sparse generated index; the view's read cost past `maxBytesPerLayer` at the fit of the
  recorded ClinVar index ("zoom in", no data read); a map's runs read by either side; EC-186, 188, 190, 195, 196, 206,
  208, 212, 220, 222, 227, 228, 229 named by tests.
- **Performance** (each with its path): a BGZF region of 64 blocks inflated, quick; a `.tbi` of ClinVar's size parsed
  under 50 ms (heavy, `itWithin`); **a VCF region of 10,000 records decoded under 100 ms**, sites, and with 100
  samples' genotypes under 200 ms (heavy; 1,000 records quick); 100,000 GFF3 lines assembled under 500 ms (heavy;
  10,000 quick); 600,000 BEDPE rows indexed by both ends under 300 ms of slices (heavy; 12,000 quick); the whole genome
  by region (M5).
- **Size**, each its own row under 30 kB and in no bundle that does not list it: `bed` about 12 kB, `gff` about 12 kB,
  `vcf` about 15 kB, each with io's BGZF and tabix (about 8 kB) where it reads them; measured at the gate. The minimal
  and default viewers unchanged.
- **`_spec` sections that change**: core Data (a section a runtime family: BGZF and tabix, BED, GFF and GTF, VCF; the
  `seq` runtime's sided pairs; Protocol's additions; names), Language (the params), Render (the ranges layer's `gene`
  line), Configuration; io's `_spec` (`bgzf.md`).
- **Other designs amended**: [files](files.md) (BED built; GFF3's Tension closed by `features`; region
  reads; `fields` by column number, `names` still refused); [sources](sources.md) (the runtimes listed,
  the designed-only ones named); [base-library](base-library.md) (BGZF to io with tabix; core's
  `ungzip`); [protocol](protocol.md) (the additions); [security](security.md) (the
  bounds table's formats rows); [performance](performance.md) (the time and size rows);
  [examples](examples.md) (the fixture rule: real slices as recorded pages).
- **Lanes** (the plan's rules: the format lanes touch runtimes, not the language; RG runs alone in the store; one lane
  at a time on the language):
  - **GC1** (stage 1, its language bump) declares GF's params (`features`, `samples`, `map.source`, `line: 'gene'`),
    with their docs and checks; no reader. The format lanes reading a new param start after GC1 merges.
  - **Stage 2, runtimes only**: **BZ** (strong, first) io's BGZF and tabix, core's `ungzip`, the allowlist,
    `runtimes/internal/` (text, lines), `findNames`' types and `data/names.ts`; then beside each other **BE** the BED
    family (whole and by region at level 1, the built `list` and `density` forms; BEDPE read whole and indexed by both
    ends); **GG** GFF3 and GTF read whole (`assembleTranscripts`, `features`, `findNames`) and `line: 'gene'` (render)
    after GC1; **VC** VCF sites at level 1. bigBed, 2bit, FASTA by region, PAF and chain have no lane (appendix).
  - **RG** (stage 3, alone in the store) takes the protocol and keys part (RR's): span rows and item extents in
    `tables.ts`, `side` in keys, items at a tile size for every list shape, `readCost` and `maxBytesPerLayer`.
  - **FR** (strong, after RG): the format lanes' by-region rest on RG's protocol: tabix levels from the index's bytes
    with exact coarse counts and `readCost` (BED, GFF, VCF), GFF's span rows by region, VCF's matrix and genotype
    classes; the `seq` runtime's sided pairs found from either end and `map.source` read by either side.
  - **LA after RG and FR** (its synteny and liftover read the `seq` runtime's sided pairs), so stage 3 is RG → FR → LA;
    the GATE of stage 2 checks BZ to VC, stage 3's checks FR.

## Real slices to record (each fetch asked first, Gate 0, 7)

Recorded as **pages** of the published file (core's 64 KiB pages the runtime read, keyed by file and page; a read the
recording did not make fails naming its page), so the container written by the real tools is what the test reads; a
change to a runtime's reads is a re-recording. Each with its `sources.json` row and licence, checked when asked. Around
TP53 (GRCh38 chr17:7,668,402-7,687,550), the p53 dataset's gene. Four fetches, numbered as GD's.

| # | Format | File to fetch | Pages | Licence to check |
|---|---|---|---|---|
| 1 | BED | UCSC hg38 `cytoBandIdeo.txt.gz` (Gate 0, 6's ask) | whole (a few kB) | UCSC |
| 2 | VCF sites | NCBI ClinVar `clinvar_20240107.vcf.gz` and `.tbi` (GRCh38 archive) | the `.tbi` whole (about 520 kB), the header, TP53's chunks | NCBI, public domain |
| 7 | GFF3 | Ensembl GRCh38 `.chr.gff3.gz` | TP53's gene lines cut by the script (a small whole file, `##sequence-region` kept) | Ensembl |
| 8 | GTF | GENCODE `annotation.gtf.gz` | TP53's lines cut | GENCODE |

The repository's current Seq\* recording of ClinVar (by its own reads) cannot serve page reads as recorded, so it is
not reused.

## Decided by the author

2026-11-21: every gap as recommended (the author on the drafts: "generally reasonable and we can fold in the reviews";
on the whole-designs review's twelve recommendations: "yes with the recommendations"; on the fold-in: "yes"). An item
the author decided by name carries its own date.

1. **BGZF and tabix (TBI and CSI) in `@molstar/io/common/bgzf` and `common/tabix`** with the surface above (lines as
   bytes, extents by the index's rule, no `.gzi`), and `seqx-core`'s io allowlist widened to them: **decided yes**.
2. **Core's `ungzip` fixed for multi-member streams by splitting members** (BZ's test first): **decided yes**.
3. **One opt-in entry a format family**, `gff3` taken from the default by `withAdapter` (the default unchanged): **decided yes**.
4. **BED's extra columns named by the built `fields` with a column number** (`{ stain: '4' }`), no new param: **decided yes**.
5. **`source.features`** (`transcripts` by default, `genes`, `rows`), a line's role by its type (`transcriptTypes`), an
   intron-only tile giving a span row; a GTF without transcript lines exact only when read whole: **decided yes**.
6. **The ranges layer's `line: 'gene'`** (thick CDS, thin UTR, intron hairline with chevrons), built in GG: **decided yes**.
7. **No dynamic transcript nodes with codon maps from a GFF** this slice (transcripts are items): **decided yes**.
8. **Coarse counts of a tabix source exact or absent**, its levels listed per file from its index's bytes (a level
   only when every tile fits `maxRegionBytes`), no estimate drawn: **decided yes**.
9. **GF's language additions ride GC1's bump** (`features`, `samples`, `map.source`, `line: 'gene'`): **decided yes**.
10. **`findNames` as an optional capability**, served by a GFF3 or GTF read whole (GTF's `gene_name`, `gene_id`,
    `transcript_name`, `transcript_id` beside GFF3's `Name`, `ID`, `Alias`: the fold-in, the author 2026-11-21: "yes"); tabix
    sources none, trix later: **decided yes**.
11. **Real slices as recorded pages** of the four files listed, each fetch asked: **decided yes**.
12. **The anchor base kept** in a VCF allele's extent (the index's rule; igv trims): **decided yes**.
13. **Genotype classes** (`0`, `1`, `2` with a haploid call, `3` partly called, `255`) in RR's matrix codes, the copy
    count a `field: 'GT'` value: **decided yes**.

Taken (the author, 2026-11-21): the three runtimes built and the rest designed only; no `.gzi` in io; BEDPE indexed by
both ends after a whole read, its held rows within a byte budget, no doubled BEDPE; the `seq` runtime's sided pairs
from either end; a view's read cost bounded; a host's bigWig as the whole-genome remedy; the lanes BZ → BE, GG, VC →
RG → FR → LA. Cut (the author, 2026-11-21): the Manhattan family, `split` and EC-224's join to the close-out slice;
MAF; fetches 3 to 6, 9 and 10.

## Appendix: designed only

No runtime this slice; a host brings its own through `withAdapter`, or a later slice builds one when a use asks. What
each would take, in short:

- **bigBed** (`runtimes/bigbed`): the BBI container shared with bigWig (moved to `runtimes/internal/bbi.ts` then);
  items `chromId, start, end, rest`, fields typed by the autoSql (else BED's names and `fields` by number), BED12's
  lists as pieces; items at any listed level as a tile width; `extraIndex` B+ trees as `findNames`. Its zoom levels
  count covered bases, not items, so a density from them needs a `measure` on `RangesDensity` (`items`, `coverage`),
  added with the runtime.
- **2bit** (`runtimes/twobit`): versions 0 and 1, either byte order; lengths live in each record's header, so every
  length at open costs a read a header (4 KiB pages through `openRange`'s page size, bounded; past the bound none
  listed, the assembly's lengths inline or `chrom-sizes`); letters by region with N blocks applied; mask blocks found
  by binary search on disk; N and mask blocks served as `ranges` (`gap`, `masked`).
- **FASTA by region** (`runtimes/fasta`, in place of the built one by `withAdapter`): a `.fai`'s offsets, each read
  checking the line structure it spans (a stale index `malformed` for that region); bgzipped with a `.gzi`, which
  needs io's `parseGzi` and `BgzfReader.readUncompressed` (added with the runtime, through a source's `block_index`).
- **PAF and PIF** (`runtimes/paf`): side `a` the query; an alignment an item, its ungapped blocks its pieces (rows
  sharing `row`); CIGAR `M`, `=`, `X` runs, `I`, `D` gaps, `N` a gap on `b`, `S`, `H` outside the range (EC-212), at
  most `maxCigarOps` a line; on strand `-` the query walked backwards from `qEnd`; a line without `cg:Z` one row of
  unequal ends, never a run. PIF (each line twice under tabix, `q` and `t` prefixed, `t` lines with `I` and `D`
  swapped: to verify on a PIF made from a PAF) reads by region from either end.
- **Chain** (`runtimes/chain`): side `a` the `t`; blocks `size dt dq` walked from `tStart`, `qStart`; `qStrand -`
  counted on the reverse strand; every block inside its ranges and sizes, else the chain left out; read whole and
  indexed by both ends, as BEDPE.
