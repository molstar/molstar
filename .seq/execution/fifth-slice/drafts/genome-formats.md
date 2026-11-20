# The genome formats: a runtime each, BGZF and tabix in the base library

Draft of lane GF (the fifth slice, stage 0) for Gate 0's agreement; on agreement it moves to `design/` and amends
[files](../../../design/files.md) (BED, GFF3's tension closed, `fields` by column number), [sources](../../../design/sources.md)
(the genome runtimes), [base-library](../../../design/base-library.md) (BGZF and tabix in `@molstar/io`),
[protocol](../../../design/protocol.md) (an item's extent, a density's measure, names), [security](../../../design/security.md)
(the formats' bounds), [performance](../../../design/performance.md) (their rows) and RR's draft
([region-reads](region-reads.md): sided pairs sources, items at a tile size for ranges and matrices) (below,
Consequences). Built on Gate 0's answers as recommended ([plan](../plan.md)): BGZF and tabix in `@molstar/io`, beside its
readers, a public surface agreed here (2); generated fixtures first and one real slice a format, each fetch asked (7);
cytobands synthetic in tests, UCSC's asked (6). Built on GA's draft ([genome-coordinates](genome-coordinates.md):
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
genomes through a PAF or a chain file (LX's figures). Cost follows pixels ([performance](../../../design/performance.md)):
a region read costs the region, every claim a file makes is bounded before it allocates
([security](../../../design/security.md)).

## Precedents

| Precedent | Taken | Not taken |
|---|---|---|
| htslib's specifications: BGZF (SAM spec §4.1: gzip members of at most 64 KiB with a `BC` extra field, virtual offsets, the 28-byte EOF block), tabix (`.tbi`: bins, chunks, a 16 kb linear index, the column header, `0x10000` zero-based), CSI (`min_shift`, `depth`, a bin's `loffset`), `tbx_parse1`'s extent of a VCF line (`POS` to `POS + len(REF)`, `INFO/END` when present) | every rule as specified; CSI beside TBI (chromosomes past 512 Mb: wheat, axolotl); a line's extent by the index's own rule, so a region read agrees with the index at its edges | writing (bgzip, tabix): tests write their own fixtures with a test-only writer |
| htslib's VCF 4.4: header lines, `##contig` lengths, `##INFO`/`##FORMAT` `Number` and `Type`, `GT` with `/` and `|`, symbolic alleles, `*`, `<NON_REF>` | typed INFO fields by the header, a row per ALT allele (EC-196), genotypes as allele dosage | phasing (a Tension), BCF |
| UCSC: bigBed (the BBI container shared with bigWig; autoSql; `extraIndex` B+ trees on a field), 2bit (versions 0 and 1, N and mask blocks), BED (`bed N+M`, BED12 blocks, `thickStart`, `itemRgb`), chain (`chain` headers and `size dt dq` blocks, `qStrand -` counted on the reverse strand), liftOver | every container rule; a field's name by autoSql; `extraIndex` as the name index; BED's `N+M` as `fields` by column number; chains as the runs of one map (GA) | `itemRgb` as a colour without a theme (a field here, Tensions); bigBed's untyped extras guessed (autoSql or `fields` say) |
| minimap2's PAF (12 columns, 0-based, `cg:Z` CIGAR, `tp:A`, `cs:Z`) and JBrowse's PIF (a PAF indexed by both ends: each line twice, names prefixed `q` and `t`, bgzipped and tabix-indexed) | PAF read whole; PIF read by region from either end through tabix (RR's "pairs found from either end") | `cs:Z` (the CIGAR suffices for runs) |
| JBrowse 2's readers (`@gmod/bgzf-filehandle`, `@gmod/tabix`, `@gmod/vcf`, `@gmod/gff`, `@gmod/bbi`, `@gmod/twobit`, `@gmod/indexedfasta`; the PAF, PIF and chain adapters; trix name search) | an adapter a format; GFF's `Parent` assembled into nested features, discontinuous features by `ID`; genotypes parsed only when asked; bigBed's `searchExtraIndex`; bgzipped FASTA with `.fai` and `.gzi` | the byte-size estimate before a read ("force load", EC-171; RR refuses it too); trix (a later format, Gaps) |
| igv.js: tabix and BGZF of its own, a feature's `visibilityWindow` (VCF and BED zoomed out show "zoom in"), GFF and GTF combined into transcripts by `Parent` and `transcript_id`, an indel drawn without its anchor base | the visibility window as each tabix source's coarsest level; GTF's `transcript_id` and `gene_id` hierarchy | trimming the anchor base (here the record's REF extent, the index's rule: Decision, VCF) |
| The current Seq\* (`stream/gmod`): 2bit, bigWig, bigBed and tabix VCF over range reads, its bounds (a BGZF block at most 64 KB with its ISIZE checked, a `.tbi` run at most 256 MB, 2bit's mask blocks found by binary search on disk), ids `name@start-end` and `pos:ref>alt`, a VCF's `estimate` from the index's bytes, recorded ranges of hg38 files | the readers' rules and bounds, the ids, recording real files by the ranges read | the estimate drawn as density (what is drawn is checked against the data); a range recording keyed by its old reads (here core's pages) |
| This codebase: the bigWig runtime (BBI header, chromosome tree, R-tree, `walk`'s bounds, `BigWigLimits`, `DecompressionStream`), core's range reads (`openRanges`, 64 KiB pages, `AssetRangeError`), core's `inflateStream` and `ungzip`, the FASTA runtime (`.fai` checked), the records runtime and its GFF3 preset, `withAdapter`, the `seq` runtime's join, `delayedRuntime` | the BBI code shared by bigWig and bigBed; one adapter a format through `withAdapter`; the records' issue rules (EC-86) and `limits` | a second inflater, a second page cache |

### Interactions in the precedents

| What a user does there | Here |
|---|---|
| open a file by URL with its index (all) | built: a `source` with `uri` and `index` (`.tbi`, `.csi`, `.fai`), never `<uri>.tbi` assumed (EC-186); a `.gzi` by `block_index` (new) |
| see a gene's transcripts, exons, CDS and UTRs (JBrowse, IGV, UCSC) | built: GFF3 and GTF as transcripts (`features: transcripts`), drawn by the ranges layer's `line: gene` (GG) |
| collapse transcripts to genes, or expand (IGV's collapsed and expanded, JBrowse's display modes) | built: `features: genes` or `transcripts` by patch; `rows` for the file's lines |
| click a variant: its INFO in a popup (all) | built: the hover card reads the fields (built) |
| filter variants by an INFO value (`CLNSIG`, `AF > 0.01`) (JBrowse) | built: the layer's data filter (the query sublanguage, EC-225) |
| genotypes of many samples, sorted by a field (IGV) | built: RR's `matrix` layer over the VCF's matrix; `samples` narrows them (new) |
| "zoom in to see features" past a visibility window (IGV, JBrowse) | built: a tabix source's coarsest level (listed per file from its index's bytes), past it the layer's too-many-tiles state, never an estimate |
| search a gene by name (IGV, JBrowse's trix, UCSC) | built in data: `findNames` on sources that index names (read whole, bigBed's `extraIndex`); the box and the jump are GA's locus fallback (GC3) |
| show N gaps and soft-masked repeats (UCSC) | built: a 2bit serves its N and mask blocks as `ranges` |
| lift a region from one assembly to another (UCSC liftOver) | built in data: a chain source as a map's runs (`map.source`); the figure LX's |
| colour items by `itemRgb` (UCSC) | refused here: a field; a theme reading an RGB text is themes' (Tensions) |
| export the features in view (JBrowse) | refused: data export is a host's |

## Decision

### Whose each choice is

By [axis](../../../design/axis.md)'s rule. What a file means is the spec's (a source's params: its `index`, its
`fields`, its `features`, its `samples`); how much a read may cost is the runtime's (an option of
`<format>Runtime(over)`, as `bigWigRuntime`) or the store's (`limits`, `maxTileBytes`); a format's conventions (BED's
base 0, GFF's base 1) are fixed by the format, never a param a writer can get wrong (EC-195). No new port: each format
is a registry entry of the existing kind ([sources](../../../design/sources.md), [seq-data-ports](../../../design/seq-data-ports.md)).

### One runtime a format family, each an opt-in entry

| Entry (`@molstar/seqx-core/…`) | Exports | Formats | Shapes | Read |
|---|---|---|---|---|
| `runtimes/bed` | `BedRuntime`, `bedRuntime` | `bed`, `bedgraph`, `bedpe`, `chrom-sizes` | ranges; values; pairs; coverage only | whole, or by region with a tabix `index` (not `bedpe`, not `chrom-sizes`) |
| `runtimes/bigbed` | `BigBedRuntime`, `bigBedRuntime`, `openBigBed` | `bigbed` | ranges (items, coverage density) | by region (range reads) |
| `runtimes/twobit` | `TwoBitRuntime`, `twoBitRuntime`, `openTwoBit` | `2bit` | letters; ranges (N and mask blocks) | by region |
| `runtimes/fasta` | `IndexedFastaRuntime`, `indexedFastaRuntime` | `fasta` (in place of the default by `withAdapter`) | letters | by region with `index` (`.fai`), bgzipped with `block_index` (`.gzi`) |
| `runtimes/gff` | `GffRuntime`, `gffRuntime`, `assembleTranscripts` | `gff3` (in place of the records preset by `withAdapter`), `gtf` | ranges | whole, or by region with a tabix `index` |
| `runtimes/vcf` | `VcfRuntime`, `vcfRuntime`, `readVcfHeader` | `vcf` | ranges (a row an ALT allele); matrix (genotypes, FORMAT values) | by region with `index`; whole without |
| `runtimes/paf` | `PafRuntime`, `pafRuntime` | `paf`, `pif` | pairs (a map's runs when equal-length) | `paf` whole; `pif` by region from either end |
| `runtimes/chain` | `ChainRuntime`, `chainRuntime` | `chain` | pairs (a map's runs) | whole |

- **Why an entry a family, opt in**: a reader of a spec, the minimal viewer and the default viewer pay for none
  ([performance](../../../design/performance.md), Size: "a format runtime … loaded only when composed"; the plan: "the
  minimal viewer not grown"). A host lists what it reads, as `BigWigRuntime` today. Why families: BED, bedGraph, BEDPE
  and chrom.sizes share one tab-line reader and one row of size; PAF and chain do not share one (two small entries).
- **Two formats the default viewer already reads** (`fasta` whole, `gff3` flat through the records preset) are taken
  over by the genome runtime when a host composes it with `withAdapter` (built: one adapter a format, the default
  narrowed, never order-dependent). Without it the built behavior stands. Why not change the defaults: the default
  viewer is at its size bound, and a protein's FASTA or a plasmid's GFF3 needs no region read or hierarchy.
- **Internals shared, not published**: `runtimes/internal/bbi.ts` (the BBI container: header, zoom headers, chromosome
  B+ tree, R-tree `walk`, block inflation; bigWig is moved onto it, no change to its surface),
  `runtimes/internal/lines.ts` (a tab-line source read whole or by region through tabix, below),
  `runtimes/internal/text.ts` (bytes to text: plain, gzip or BGZF by magic, bounded). The package exports
  `"./runtimes/internal/*": null`, as core's `util/internal`. Why: core's `./*` rule publishes every file otherwise.
- **Each with its synthetic twin** (`delayedRuntime(<Runtime>)` on a fake clock), as sources.md requires.

### BGZF and tabix in `@molstar/io` (Gate 0, 2): the public surface proposed for the gate

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
export interface GziIndex { readonly compressed: Float64Array; readonly uncompressed: Float64Array }
export function parseGzi(bytes: Uint8Array, maxEntries: number): GziIndex;
export interface BgzfReader {
  /** Inflated bytes from `begin` to `end` (the end block read for its first `within(end)` bytes). */
  read(begin: VirtualOffset, end: VirtualOffset, signal?: AbortSignal): Promise<Uint8Array>;
  /** Inflated bytes `[start, end)` of the uncompressed stream, by a `.gzi` (bgzipped FASTA). */
  readUncompressed(gzi: GziIndex, start: number, end: number, signal?: AbortSignal): Promise<Uint8Array>;
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
  readonly maxLineBytes: number;               // one line (16 MiB: a VCF line of 100,000 samples' GT:AD:DP:GQ:PL)
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
  `_spec` gains `bgzf.md` (BGZF, `.gzi`, TBI, CSI: their rules and bounds); `seqx-core`'s restricted imports
  (`biome.json`) gain the two entries beside the text CIF parser: an architecture rule widened by the gate's decision,
  not to make an import pass. A v6 changelog entry (an addition, no break).
- **Core's `ungzip` fails on a multi-member stream** (the review verified it): `gzipDeclaredSize` reads the last
  member's ISIZE (BGZF's EOF block: 0), so the platform path throws on its first chunk; the platform decoders refuse
  concatenated members anyway (Node: "trailing junk"); the fallback `ungzipBlocks` goes through `_inflate`'s raw path
  and throws too. Only where `DecompressionStream` is missing does it come back short (the first member). **The fix
  splits members, never resizes a buffer**: each member found (by its `BC` field's BSIZE when present, else where the
  JS inflater's deflate stream of the member ends), inflated as one member with its own ISIZE, the outputs joined within
  the caller's bound. BZ lands the test first (BGZF, two plain members, a member cut short), then the fix. The genome
  runtimes read BGZF raw (no `decompress`) through io either way. Gaps, 2.

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
  size, extended to `ranges` and `matrix`: below, for RR word for word), and above 1 a `density` (ranges) or `bins`
  (matrix) **counted exactly from the records of that tile**. At open the runtime lists a level only when **every tile
  of it on every system fits `maxRegionBytes`** (`regionBytes`: the chunks' compressed bytes, read from the index,
  never from the data). So ClinVar (about 30 bytes a kb) lists all four; a 1000 Genomes chromosome (tens of MB a Mb,
  mostly genotypes) lists `1` and `16` (16 kb tiles, views up to about 4 Mb under `maxTilesPerWant`), and a wider view
  is the layer's built too-many-tiles state ("zoom in"). Why per source and from bytes: one runtime reads sparse and
  dense files, and the bytes a tile costs are known exactly before a read, where its rows are not. **Why counts exact
  or absent**: a density estimated from the index's bytes (the current Seq\*'s `estimate`, JBrowse's gate) would draw a
  guess as data; the listed coarsest level is igv's visibility window, made per file. A dense file's sites wider than
  that are a bigBed's (a host converts; Tensions). `samples` narrows the states held, not the bytes read (genotypes are
  inside the compressed lines); a tile's states stay within RR's `maxTileBytes`.
- **Issues** (`seq.source.<format>.*`, one a cause with a count and the first ten places, EC-86; severities EC-87): a
  line that does not parse is left out (`.line`, warning, with the line number or virtual offset); a structure that
  cannot be read refuses the source (`seq.source.malformed`); a bound passed is `too-large` for a region and
  `seq.source.too_large` at open, hinting the bound and the remedy (an index, a coarser level, a smaller file).
- **Range reads**: `AssetReader.openRange` (built); a host without it reads the file whole within `limits.maxBytes`
  (bigWig's rule). A failed range read (`AssetRangeError`: the file changed under the reader) fails the region
  (`failed`, EC-174), never wrong rows.
- **A format's conventions are fixed** (`bed`, `bedgraph`, `bedpe`, `vcf`, `gff3`, `gtf`, `paf`, `pif`, `chain`): the
  source's built `base` and `ends` (files' defaults `1` and `inclusive`) are not read for them, said once if written
  (`seq.source.<format>.transform`, info), as the GFF3 preset already does. Why: a writer setting `base: 1` on a BED
  would move every row by one (EC-195).
- **Names are the file's**: a source's systems are the names its file writes (`1` or `chr1`); GA's names rule maps
  them onto an assembly (EC-188, EC-189). A system's length is the file's where it says (`##contig`,
  `##sequence-region`, 2bit, `.fai`, BBI), else absent: the assembly's stands.

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
  `pairs` list; **read whole**, an `index` ignored and said (RR's one-end rule: a `.tbi` finds a pair by its first end
  only); within `limits`, then RR's pairs summaries. `id` the name unless `.` (RR: `.` is no identity), else `row`.
- **chrom.sizes**: `name length` lines; a coverage of systems and lengths, no shapes. Why here: GA's assembly takes "a
  chrom.sizes table" as its `source`, and a draft genome's 2bit is too many record headers to read (2bit, below).

### bigBed (`runtimes/bigbed`)

- **The BBI container** shared with bigWig (`internal/bbi.ts`): magic `0x8789F2EB` either way round, the chromosome
  tree, the R-tree, blocks inflated one after another within `maxInflatedBytes` (bigWig's bounds and reasons, unchanged).
- **Items**: `chromId, start, end, rest` with `rest` a NUL-terminated tab-separated string; fields named and typed by
  the file's autoSql (`uint`, `int` int32; `float`, `double` float32; `string`, `lstring`, `char[n]` categories; a
  list type (`int[blockCount]`) its text, BED12's lists read as pieces), else by BED's names and `fields` by number. A `fields`
  map renames as for any source (EC-190: `{ label: 'geneName' }`). Pieces, ids and `thick_*` as BED12.
- **Levels**: `1` and items at any listed level as a tile width (RR's items rule, amended); the file's zoom levels at
  their reductions as `density` with **`measure: 'coverage'`** (new, below): a zoom record covers bases (its summed
  depth prorated over the bins it meets, bigWig's rule), not items. Why say so: a density strip of coverage read as a
  count of items would lie by the items' lengths; the hover says "coverage", the strip draws either.
- **Name index**: each `extraIndex` B+ tree (`bedToBigBed -extraIndex=name`) is a `findNames` field (below).

### 2bit (`runtimes/twobit`)

- **Header and index**: the signature read either way round; versions 0 (32-bit offsets) and 1 (64-bit, files over
  4 GiB) both read; the index (`nameSize`, name, offset) read in one call within `maxIndexBytes`.
- **Lengths** live in each record's header (its first 4 bytes), not the index: **every sequence's is read at open**,
  in offset order, in batches of 1,024 headers a call, through a range reader of 4 KiB pages (`AssetReader.openRange`
  gains `{ pageSize }`, a `packages/seq` addition passed to core's `RangeReadOptions`, built), so a header costs a 4 KiB
  page, and neighbouring small scaffolds share one. Bounded by `maxLengthBytes` (16 MiB of pages: 4,096 pages; hg38's
  455 sequences about 1.8 MB, a draft of 100,000 small scaffolds packed in a few hundred pages). **Past the bound no
  length is listed** (never half a list): the coverage names every sequence without a length, said
  (`seq.source.twobit.lengths`, warning: "give its chrom.sizes as the assembly's `sizes`"), and its letters are still
  served by region. **An assembly never draws without lengths**: GA gains `sizes` and an error when a sequence has none
  (below, word for word). Why every length: GA's remedies for a draft (`include: { minLength }`, `order: length`, the
  closing run) need them all; a 2bit for letters and a chrom.sizes for lengths is UCSC's own pair.
- **Letters by region**: the packed bases of `[start, end)` (four a byte, `T C A G`, the first in the high bits), N
  blocks applied (`N`), upper case. A record's N blocks are read with its header (hundreds); its **mask blocks are
  found by binary search on disk** (chr1 of hg38 has 331,893: 2.7 MB never read whole).
- **Ranges**: `N` blocks (category `gap`) and mask blocks (category `masked`) served as `ranges` on request, so a
  genome's gaps and repeats are a track of their own. Why not lower case in the letters: the letters are one alphabet
  for every reader (schemes, translation, search); the mask is data.

### FASTA by region (`runtimes/fasta`)

- **With an `index`** (`.fai`) and range reads: a region is `offset + ⌊s / lineBases⌋ · lineWidth + s mod lineBases`
  to the same for `e`, line breaks removed. The index is parsed as built (`parseFai`); since the file is not read whole,
  **each read checks the line structure it spans** (a line break exactly where the index puts one, none elsewhere),
  else that region is `malformed` (a stale index), and at open the first record's header line ends where its offset
  starts (one small read). Coverage from the `.fai` alone, so a FASTA index is an assembly's `source` (GA) without a
  byte of the FASTA.
- **Bgzipped** (`.fa.gz` with `.fai` and `.gzi`, as Ensembl's `dna_index` publishes): the `.fai`'s offsets are of the
  uncompressed stream; `block_index` (new source param) names the `.gzi`; read through `readUncompressed`. A bgzipped
  FASTA without a `.gzi` is read whole within `limits`.
- **Without an index** the built runtime's behavior: read whole, the index built (`faidx`).

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
  transcripts): a protein's coordinate system from a GFF's CDS is a figure no Genomes view needs yet (Gaps, 7).
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
  stable across tiles); the `ID` a field `name`, the label. Fields: `ref`, `alt`, `qual` (float32), `filter` (category), `category` (`SNV`, `MNV`,
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
- **Dense files** (1000 Genomes, a biobank): the levels the index's bytes allow (above), so a dense chromosome is read
  in 16 kb tiles of a few hundred kB; a wider view is "zoom in", never an 80 MB read. Genotypes are parsed only for a
  `matrix` request, and only for the `samples` held; the line's bytes are read either way (they are in the compressed
  block), which the level rule already counts.
- **One file a chromosome** (1000 Genomes, gnomAD) is one source a file this slice; a source over several files is
  open (Tensions).

### PAF and chain files: pairs, and a map's runs (`runtimes/paf`, `runtimes/chain`)

- **Side `a` is the side the alignment maps from**: PAF's query (`minimap2 target query` maps query onto target),
  chain's `t` (liftOver's chain files are written old-to-new, `tName` the old assembly). Why one rule: a map read from
  either is `from` side `a` `to` side `b` with no per-format flag.
- **An alignment is an item, its ungapped blocks its pieces**: a pairs row a block (`aStart`, `aEnd`, `bStart`, `bEnd`
  of equal length, `aStrand` 1, `bStrand` ±1), all of one alignment sharing `row` (amends RR: rows sharing `row` with
  other starts are blocks of one alignment, as ranges' pieces); fields: PAF's `matches`, `block_length`, `mapq`, `tp`;
  chain's `score`, `chain_id`. A PAF line without `cg:Z` is one row of unequal ends (a pair, never a run: a map refuses
  it, `seq.source.paf.cigar`). CIGAR operations `M`, `=`, `X` are runs, `I`, `D` gaps, `N` a gap on `b`, `S` and `H`
  outside the aligned range (EC-212); at most `maxCigarOps` (1,048,576) a line, which fits `maxLineBytes` (16 MiB).
  **On strand `-`** the CIGAR walks the target forward from `tStart` and the query backwards from `qEnd` (minimap2's
  rule; PAF's query coordinates are of the forward strand), so each run is `[qEnd − q − len, qEnd − q)` on `a` against
  `[tStart + t, + len)` on `b`, strand −1.
- **PIF's `t` lines** carry the line with query and target swapped, and their CIGAR with `I` and `D` swapped (JBrowse's
  PIF writer; PF verifies it on a PIF made from a PAF before relying on it); the runtime swaps both back, so either
  copy gives the same row with side `a` the query, and a pair met from both ends is one row (RR's dedupe).
- **Chain arithmetic**: blocks `size dt dq` walked from `tStart`, `qStart`; `qStrand -` counts `q` on the reverse
  strand (`qSize − end`, `qSize − start`); every block inside `[tStart, tEnd)` and `[qStart, qEnd)` and within the sizes,
  else the chain left out (`.line`). Each chain header's `tSize` and `qSize` are its systems' lengths.
- **Read whole** within `limits` (a chain file is small: hg19ToHg38's is about a megabyte gzipped), indexed by both
  ends once (RR's one-end rule), then RR's pairs summaries (its pyramid within `summaries.maxBytes`). **`pif`** reads by
  region from either end (each line twice, `q` and `t` prefixed, under tabix): the large PAF's path.
- **Sided systems** (amends RR): a pairs source between two genomes may name `chr1` on both sides, distinct
  sequences. Its coverage marks each system's side (`SystemInfo.side: 'a' | 'b'`), a request names the side it reads
  (`RegionRequest.side`), and the store's key carries it in its system part. A `pairs` node's `a` and `b` (RR's spec
  params, an assembly each) map each side's names through GA's names rule. A source over one genome (BEDPE, self
  synteny, Hi-C) is unsided, as RR has it.
- **As a map** (`map.source`, new language param): a `map` with `method: segments` may name a pairs source instead of
  inline `segments`; its runs are the source's rows of equal-length ends, read by region of the `from` side as any
  large map (building-blocks.md: "a genome-scale map is held by region"; GA: one map over every sequence pair). A row
  of unequal ends is refused as a run (counted). Why the pairs shape: a run is a pair of equal length, so one runtime
  serves the ribbons (LX) and the liftover (EC-212) from one read.

### A name index (EC-228)

- **A capability, not a shape**: `OpenedSource.findNames?(text, { max, prefix }, signal): Promise<NameHit[]>`, with
  `Coverage.names` (the fields searched) and `NameHit { system, start, end, id?, label, field }`. Case-insensitive,
  exact first, then by prefix when asked; at most `max` (20) hits.
- **Who serves it**: a source read whole (BED's `name`, GFF's `Name`, `ID`, `Alias`, `gene_name`, `transcript_name`,
  `gene_id`, `transcript_id`; VCF's `ID`), from a map built once at open; bigBed by its `extraIndex` B+ trees (a key
  read and its block filtered). A tabix source serves none (a text index such as trix is a later format: Gaps, 12).
- **Who calls it**: `findNames(sources, text, options, signal)` (`data/names.ts`, a function over the opened sources
  a viewer's source sync holds: not in `store.ts`, so no store lane is touched), merging hits in source order; GA's
  `(locus :text …)` falls back to it when no sequence matches (GC3's: an asynchronous focus, cancelled by the next),
  and a hit focuses the item's address (`(item :id ENST…)` on its source) or, without an `id`, its range.

### Protocol additions (`packages/seq`, public API, recorded in the lanes' changelog entries)

Each with its owning lane (below, Lanes and the plan's amendment):

- `RangesList.itemStart?`, `itemEnd?` (`Int32Array`): each row's item's whole extent, when rows are pieces (BED12,
  transcripts); and **span rows** (`start = end` with an extent: the item, no piece in this tile), which `sliceTable`
  emits for an item whose pieces lie outside the slice. Why: a reader of region data knows an item's ends from any
  tile (today `extents` exist for inline data only, so region items get no arrowhead), and an intron-only tile still
  draws its transcript (EC-229: flanks of streamed genes need the extent). RG (tables and `checkAnswer`).
- `RangesDensity.measure?: 'items' | 'coverage'` (absent: items): what `count` counts. RG.
- `SystemInfo.side?`, `RegionRequest.side?`, the key's sided system part (RR amended). RG.
- RR's `RegionRequest.items` for `ranges` and `matrix` lists too (keys change: RR's `i<level>` part). RG.
- `OpenedSource.findNames?`, `Coverage.names?`, `NameHit` (types in `protocol.ts`, no store change) and
  `data/names.ts`. BZ (types only, before RG; RG merges over them).
- `AssetReader.openRange(uri, signal, options?: { pageSize })`. BZ.
- Field names fixed by convention and documented in the Data spec: `thick_start`, `thick_end`, `cds_phase`, `phase`,
  `item_rgb`, `gene`, `gene_name`, `category`, `name`.
- Language (GC1's bump, the only lane on the language in stage 1): a source's `features`, `samples`, `block_index`;
  an assembly's `sizes` (GA's); `map.source`; the ranges layer's `line: 'gene'`. Declared, documented and checked by
  GC1; read by the format lanes after GC1 merges.

### Bounds (security.md: a file is untrusted input)

Each refused before it allocates; `malformed` for a structure, `too-large` for a region (a coarser level, a smaller
region), `seq.source.too_large` at open.

| Format | Bound | Default | Over it |
|---|---|---|---|
| BGZF | a block's compressed and inflated size; its ISIZE and CRC-32 | 65,536; the platform's check | malformed |
| BGZF | a read's inflated bytes, its blocks together | `maxInflatedBytes` 64 MiB | too-large |
| tabix | the index inflated; its references; a reference's bins (TBI: 37,450 and the metadata pseudo-bin), chunks, a chunk `begin ≤ end` and inside the file | `maxIndexBytes` 64 MiB; `limits.maxRows` | malformed (the source refused) |
| tabix | a region's chunks, compressed, together (a level not listed when any of its tiles passes it); one line; the lines of an answer | `maxRegionBytes` 8 MiB; `maxLineBytes` 16 MiB; `limits.maxRows` | too-large; malformed; too-large |
| all | a sequence's length | 2^31 − 1 (int32 positions) | that sequence refused at open, said |
| VCF | the header; samples; ALT alleles a record; INFO fields | `maxHeaderBytes` 16 MiB; 100,000 (a line of them within `maxLineBytes`); 256; 256 | malformed; malformed; the record left out; the rest not fields (said) |
| bigBed | bigWig's (`BigWigLimits`); the autoSql; fields; a BED12 item's blocks | as built; 64 KiB; 256; 10,000 | as built; malformed; malformed; the item left out |
| 2bit | the index; sequences; the pages read for lengths; a record's N blocks; a region's packed bytes | `maxIndexBytes` 8 MiB; `limits.maxRows`; `maxLengthBytes` 16 MiB; `limits.maxRows`; `limits.maxBytes` | malformed; malformed; no lengths listed (said; GA's `sizes`); malformed; too-large |
| FASTA | `.fai` and `.gzi` lines; a region's bytes | `limits.maxRows`; `limits.maxBytes` | malformed; too-large |
| GFF, BED, PAF, chain read whole | bytes (after inflating), lines | `limits.maxBytes` 64 MB, `limits.maxRows` 1,000,000 | `seq.source.too_large`, hinting an index (`pif` for a PAF) |
| GFF | a `Parent` chain's depth; attributes a row | 16; 256 | rows kept flat (said); the rest dropped (said) |
| PAF | CIGAR operations a line | 1,048,576 | the line left out |

## Edge cases this design takes or assigns

- **Taken here**: EC-186 (an index named, never `<uri>.tbi`; `.gzi` by `block_index`), EC-188 (lengths from
  `##contig`, `##sequence-region`, 2bit, `.fai`, chain headers, chrom.sizes), EC-190 (`fields` by name or column number, autoSql),
  EC-195 (a format's base fixed, converted at the edge, EC-1), EC-196 (ids `name@start-end`, `pos:ref>alt`, a row a
  BED12 block and an ALT allele), EC-206 (strand on every pairs row: PAF's, chain's `qStrand`), EC-208 (pairs read by
  region: `pif`; PAF and BEDPE whole within bounds), EC-212 (chains and CIGARs with `S`, `H`, `N`), EC-220 and EC-222
  (the VCF's matrix: missing as 255, a long deletion one item), EC-227 (gene, transcript, exon), EC-40 and EC-92
  (every reader fuzzed without a throw), EC-51 (a row past a known length clipped, counted), EC-86 and EC-87 (issues),
  EC-169 and EC-174 (pending, failed and `AssetRangeError`), EC-171 (its outcome for tabix: exact counts within bytes,
  else "zoom in").
- **Its data part here, the rest assigned**: EC-228 (`findNames` here; the box and the focus GA's locus fallback,
  GC3), EC-189 (a file's names here; aliases GA's names rule), EC-207 (a chain as one map's runs here; paths GA's,
  GC1), EC-209 and EC-210 (data here; ribbons LX's).
- **Assigned**: EC-213 (MAF: cut first, Gate 0, 1), EC-223 (`.hic`, `.mcool`: not in the formats, RR's Gaps 10),
  EC-224 (keyed results: RR's `rows`), EC-216 and EC-219 (GWAS points and LD: GD's data, a BED or TSV here).

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
3. **VCF past its coarsest level** (`[1, 16, 256]`: a view wider than about 67 Mb): the layer's too-many-tiles state
   (built: the band hatched, "zoom in to read this track"), never an estimate.
4. **bigBed zoomed out**: the zoom levels' coverage as a density strip, its hover "coverage 12.3 kb in 1 Mb";
   zoomed in, items as BED12 genes (picture 1).
5. **2bit**: letters at letter zoom; with a ranges layer over its gaps and masks, `gap` and `masked` bars.
   ```
   ACGTNNNNNNNNNNACGTACGT…
       ████████████            gap
                   ▬▬▬▬▬▬▬     masked
   ```
6. **BED with `fields: { stain: '4' }`** under a `bands` layer: the ideogram (navigators.md's picture).
7. **Pending**: as built (the band's pending hatch); **failed** (a range read refused, a changed file): the placeholder
   with its issue; **malformed** region: the placeholder, the issue naming the virtual offset or line.
8. **Dark**: tokens only; nothing a format adds is a colour.

## The interaction matrix's cells

Formats add no gesture of their own: what a file holds is drawn by the layers that read it, as files.md says. The cells
they change:

| | select | extend | add / remove | focus | back | zoom | pan | hide | edit |
|---|---|---|---|---|---|---|---|---|---|
| **an element: a transcript** | click: the transcript (every piece, `(item :id …)` on its source) | Shift+click: the hull to another | Cmd/Ctrl+click | double-click: its extent (`itemStart`, `itemEnd`) | Alt+Left | — | — | no: a filter | no: data |
| **an element: a variant allele** | click: the allele's REF extent | Shift+click | Cmd/Ctrl+click | double-click: its extent | Alt+Left | — | — | no: a filter on `category`, `filter`, INFO | no |
| **an element: an alignment's block** (PAF, chain) | click: the block, both ends lit | Shift+click | Cmd/Ctrl+click | double-click: the alignment's hull on each side | Alt+Left | — | — | no | no |
| **a track or layer** | as built | as built | `track.add` as built | as built | as built | — | — | `layer.hidden`; the data filter | the source's `features`, `samples`, `fields`: a patch, one undo step, the source reopened in place (files.md) |
| **a bookmark or link** | — | — | as built | an item's address (`(item :id ENST…)`) waits for its region (dynamic-nodes' rule) and survives a `features` switch only for `transcripts` ↔ `rows` ids it keeps (a Tension) | as built | — | — | — | — |
| **the view** | — | — | — | a name typed (a host's box, the locus text) finds through `findNames`, then focuses (GC3) | Alt+Left | — | — | — | — |

## Configuration

| Item | Whose | Default | Why the default; the alternatives |
|---|---|---|---|
| each runtime | registry entry | not listed | opt in: a viewer pays for what it reads; listed by a host |
| `fasta`, `gff3` by the genome runtimes | `withAdapter` | the defaults (whole, flat) | the built behavior; the genome runtimes by composition |
| `source.index` (built), `source.block_index` | spec params | none | where the index is (EC-186); `.tbi`, `.csi`, `.fai`; `.gzi` |
| `source.fields` by column number (built param, new reading) | spec param | BED's names by place | a `bed N+M` file names its own columns |
| `source.features` | spec param | `transcripts` | a gene model (EC-227); `genes`, `rows` |
| `source.samples` | spec param | every sample | a figure's subset; a list of names |
| `source.fields`, `types` (built) | spec params | the format's names and types | renames (EC-190) |
| `map.source` | spec param | — (inline `segments`) | a chain or PAF's runs by region |
| `maxLevel` of a tab-line source by region | runtime option | 4,096 | the coarsest candidate (4 Mb tiles); the levels listed are those the file's index says fit `maxRegionBytes` |
| `transcriptTypes` (GFF, GTF) | runtime option | the SO transcript types above, and `*_transcript` | a line's role by its type, so a tile reads as the whole; a list |
| `maxRegionBytes` | runtime option (`TabixLimits`) | 8 MiB compressed | a tile's read stays a few hundred ms remote; any |
| `maxInflatedBytes` | runtime option | 64 MiB | the store's answer bound |
| `maxIndexBytes` (tabix; 2bit) | runtime option | 64 MiB; 8 MiB | a dbSNP `.tbi` fits; any |
| `maxLineBytes`, `maxHeaderBytes` | runtime option | 16 MiB; 16 MiB | a line of 100,000 samples' `GT:AD:DP:GQ:PL`, a PIF line of 1M CIGAR operations; gnomAD's header a few hundred kB |
| `maxInfoFields`, `maxSamples`, `maxAlleles` (VCF) | runtime option | 256; 100,000; 256 | bounded columns and states |
| `maxLengthBytes` (2bit) | runtime option | 16 MiB of 4 KiB pages | every length read at open; past it none listed and GA's `sizes` |
| `assembly.sizes` (GA's) | spec param | none | a chrom.sizes or `.fai` for lengths beside a 2bit |
| `maxCigarOps` (PAF) | runtime option | 1,048,576 | a chromosome-scale alignment's CIGAR |
| `cacheBlocks` (BGZF) | runtime option | 16 | a pan back without a read |
| `findNames`' `max`, `prefix` | option of the call | 20; false | a list a box shows; prefix for typing |
| a format's base and ends; side `a` the map's from side; the anchor base kept; a line's extent by the index's rule | fixed | — | the reasons in the Decision |

Not built: all of it (BZ, BD, SQ, GG, VC, PF); bigWig's BBI moves onto `internal/bbi.ts` in BD.

## Tensions

- **The anchor base**: a VCF indel spans its REF with the anchor (the index's extent), where igv draws the deletion
  without it and an insertion between bases. Recommendation: the index's extent (a region read agrees with the
  index); a `trimmed` drawing is the variants layer's later, if the use gate asks.
- **Phasing**: `0|1` and `1|0` give one dosage; haplotypes as two rows a sample need a states system twice the
  samples. Open until a view asks (a phased panel's figure).
- **No density past the coarsest level of a tabix source**: ClinVar at the whole genome draws "zoom in", where a
  bigBed of the same data draws its coverage. Accepted: honest and bounded; a host converts a large VCF to bigBed or
  converts a large VCF to bigBed for its genome view.
- **Breakends** (`N[chr2:321682[`) are one-position items with their mate in `alt`; served as `pairs` they would draw
  as links. Open until an SV view asks (a later lane, the pairs shape exists).
- **`itemRgb`** is a field, not a colour: a theme reading `r,g,b` text is the themes design's (a `token:value` form
  over a field). Open.
- **2bit lengths past 16 MiB of pages** need a chrom.sizes (`sizes`): a draft genome's 2bit alone lists names without
  lengths, and its assembly says so. Accepted; the issue names the remedy.
- **GTF without transcript lines by region** has partial extents (the tiles' hull) and draws nothing in an intron-only
  tile. Accepted, said; read whole for exact extents.
- **A type outside `transcriptTypes`** with exon children reads as a transcript only where its pieces cover it by
  region. Accepted, said once by type; the list is a runtime option.
- **A source over several files** (one VCF a chromosome: 1000 Genomes, gnomAD): one source a file this slice. Open: a
  `uri` template by system, if GD's views need it.
- **Genotype classes** lose the copy count past diploid (`1/1/1/0` and `0/0/0/1` are both `1`); the count is `field: 'GT'`.
  Accepted.
- **A bookmark of an item across `features`**: a gene's id in `genes` is not a transcript's; switching modes leaves a
  bookmark of the other mode waiting. Accepted: bookmarks of positions survive (files.md's rule).
- **Name search over tabix sources** has no index: a host's large GFF finds nothing by name. Open until trix (Gaps, 12).
- **A GFF3 child outside its parent's extent** differs by region from the whole read (its fields absent in the tile).
  Accepted: malformed GFF3, said.

## Self-check

The claim: **this design covers what a user does with genome files.**

1. **Assume it holds.** A host composes the runtimes it reads; each opens its file by range, reads a region within
   bounds into the protocol's tables, transcripts and variants and genotypes and pairs and letters; a name finds its
   item; the files' conventions are converted at the edge.
2. **What it rests on, negated.** (a) Files are indexed; negated: a plain VCF, a PAF, a GFF without `.tbi`. (b) Files
   are well formed; negated: crafted headers, stale indexes, a file replaced under the reader, a truncated BGZF. (c) A
   file is read by position; inverse: found by name. (d) Features are flat; negated: a three-level hierarchy, a shared
   exon, a GTF without transcript lines, a prokaryote's gene with a CDS. (e) One sample; plural: thousands. (f) One
   genome; plural: two sharing `chr1`. (g) Data is dense enough to read at the fit; negated: a sparse VCF over a
   genome, or one so dense a tile is megabytes. (h) Names match the assembly; negated: `1` against `chr1`, contigs the
   assembly lacks. (i) The platform inflates gzip as one stream; negated: BGZF's members, no `DecompressionStream`. (j)
   What is drawn is counted; negated: an estimate. (k) A format's conventions are one; negated: BED 0, GFF 1, VCF's
   REF extent, chain's reverse-strand counting.
3. **The negated world.** (a) whole reads within bounds, a hint to index, `pif` for PAFs; (b) bounds before
   allocation, a line check per FASTA read, EOF said, `AssetRangeError` failing the region; (c) `findNames`; (d)
   `assembleTranscripts` by children, not types, a piece of each parent, hulls; (e) `samples`, genotypes parsed only
   when asked, RR's byte bound; (f) sided systems; (g) levels listed per file from its index's bytes, exact counts, `maxRegionBytes`;
   (h) the file's names kept, GA's rule maps them; (i) io's block reader and core's `ungzip` fixed; (j) no estimate;
   (k) conventions fixed per format and converted once.
4. **Both held.** Integrated: every item of 3 in the Decision. Standing: phasing (e); breakends as pairs; name search
   over tabix sources (c); no density past a tabix source's coarsest level (g); the anchor base (k); 2bit lengths past
   the page bound (h); GTF without transcript lines by region (d).

## Consequences

- **`@molstar/io`** (BZ): `common/bgzf.ts`, `common/tabix.ts` (the surface above), tests in `common/_test/` (blocks
  across a read's edge, the EOF block, a `.gzi`, TBI and CSI, `reg2bins` against the specification's table, a region's
  lines equal to a whole inflate filtered, every bound refused before it allocates, fuzzed blocks and indexes without a
  throw); `_spec/bgzf.md`; a v6 changelog entry.
- **`@molstar/core`** (BZ): `ungzip` over a multi-member stream (every member, within its bound), its test; a v6
  changelog entry (a fix).
- **`biome.json`** (the orchestrator's): `seqx-core`'s allowlist of `@molstar/io` gains `common/bgzf` and
  `common/tabix`.
- **`seqx-core`**:
  - `runtimes/` entries (above) and `runtimes/internal/`; the package's exports `"./runtimes/internal/*": null` (hot
    file); bigWig onto `internal/bbi.ts`.
  - protocol (`data/protocol.ts`): `itemStart`, `itemEnd`, `measure`, `findNames`, `names`, `NameHit`, `side`; the
    `data/names.ts`; RG's part (span rows, `measure`, `side`, items for ranges and matrices, `checkAnswer`, `sliceTable`)
    as the lanes below.
  - language (GC1's bump): `source.features`, `samples`, `block_index`; `assembly.sizes`; `map.source`; `line: 'gene'`;
    issues `seq.source.{bgzf,tabix,bed,bigbed,twobit,fasta,gff,vcf,paf,chain}.*`.
  - render: the ranges layer's `gene` line (GG), item extents from `itemStart`/`itemEnd` for region data, the density
    hover's measure.
- **`seqx-viewer`**: nothing but a host composing the runtimes (`withAdapter` for `fasta`, `gff3`); the find box is
  GC3's.
- **`seqx-examples`**: the test-only writers (BGZF and tabix, bigBed on bigWig's writer, 2bit, `.fai` and `.gzi`) on
  core's `PCG`, which GD's generator uses for the toy assemblies' files; the recording script for the real slices
  (below); `sources.json` rows with licences.
- **Storybook**: `Features/Genome formats/<format>` (the file's head or its index's summary and the source node beside
  the canvas and the SVG, a misread file and its issues), each with its headless twin on recorded or generated files;
  a Design page "Genome formats".
- **Tests**, per runtime: the conformance suite (`data/_test/source-suite.ts`: coverage, region answers, empty and edge
  regions, cancellation) on a generated fixture and on its recorded real slice; **a region read equals the whole read
  clipped** at every level (a property over generated files and random regions, items deduped across tiles); fuzzed
  and truncated files read without a throw, every bound refused before it allocates; per format: BED12 pieces and
  `fields` by column number; bigBed's autoSql and `extraIndex`; 2bit's N and mask blocks against a naive decode, both versions and
  byte orders; FASTA's stale index refused per region, bgzipped with `.gzi`; `assembleTranscripts` (three levels, a
  shared exon, CDS-only, GTF with and without transcript lines, a cycle, a prokaryote's gene, a type outside the list)
  and its region twin: **a tile inside an intron gives the transcript's span row**, every tiling of a generated GFF
  equal to the whole read clipped (EC-227, EC-229); VCF's extent equal to `tbx_parse1`'s per line (an `END=` on a SNV
  line too), genotype classes (haploid, partial, multi-allelic), ids of duplicates and long SVs; levels listed from a
  dense and a sparse generated index; 2bit lengths through 4 KiB pages within the bound and none past it; VCF's alleles,
  INFO types and `Number`s, dosage codes, a long deletion; PAF's CIGAR runs, chain's reverse strand, a liftover of the
  recorded chain against UCSC's own coordinates recorded at fetch; EC-186, 188, 190, 195, 196, 206, 208, 212, 220, 222,
  227, 228 named by tests.
- **Performance** (each with its path): a BGZF region of 64 blocks inflated, quick; a `.tbi` of ClinVar's size parsed
  under 50 ms (heavy, `itWithin`); **a VCF region of 10,000 records decoded under 100 ms**, sites, and with 100
  samples' genotypes under 200 ms (heavy; 1,000 records quick); a bigBed region at the level that fits under 50 ms
  (heavy, as bigWig); a 2bit megabyte of letters under 20 ms (heavy; 64 kb quick); 100,000 GFF3 lines assembled under
  500 ms (heavy; 10,000 quick); 100,000 chain blocks read under 300 ms (heavy); the whole genome by region (M5).
- **Size**, each its own row under 30 kB and in no bundle that does not list it: `bed` about 12 kB, `bigbed` about
  10 kB (with the shared BBI, bigWig's row unchanged), `twobit` about 6 kB, `fasta` about 5 kB, `gff` about 12 kB,
  `vcf` about 15 kB, `paf` about 6 kB, `chain` about 4 kB, each with io's BGZF and tabix (about 8 kB) where it reads
  them; measured at the gate. The minimal and default viewers unchanged.
- **`_spec` sections that change**: core Data (a section a runtime family: BGZF and tabix, BED, bigBed, 2bit, FASTA by
  region, GFF and GTF, VCF, PAF and chain; Protocol's additions; names), Language (the params), Render (the ranges
  layer's `gene` line), Configuration; io's `_spec` (`bgzf.md`).
- **Other designs amended**: [files](../../../design/files.md) (BED built; GFF3's Tension closed by `features`; region
  reads; `fields` by column number, `names` still refused); [sources](../../../design/sources.md) (the runtimes listed);
  [base-library](../../../design/base-library.md) (BGZF to io with tabix; core's `ungzip`);
  [protocol](../../../design/protocol.md) (the additions); [security](../../../design/security.md) (the bounds table's
  formats rows); [performance](../../../design/performance.md) (the time and size rows);
  [examples](../../../design/examples.md) (the fixture rule: real slices as recorded pages); RR's draft (sided pairs
  sources; pieces of an alignment; items at a tile size for ranges and matrices); GA's draft (the locus fallback's
  `findNames`).
- **Lanes, and the plan's amendment** (the plan's rules: the format lanes touch runtimes, not the language; RG runs
  alone in the store; one lane at a time on the language):
  - **GC1** (stage 1, its language bump) declares GF's params (`features`, `samples`, `block_index`, `map.source`,
    `line: 'gene'`) and GA's `assembly.sizes`, with their docs and checks; no reader. The format lanes reading a new
    param start after GC1 merges (BD's bigBed and SQ's 2bit need none and start at once).
  - **Stage 2, runtimes only**: **BZ** (strong, first) io's BGZF and tabix, core's `ungzip`, the allowlist,
    `runtimes/internal/` (text, lines, BBI), `findNames`' types and `data/names.ts`, `openRange`'s `pageSize`; then
    beside each other **BD** BED family (whole and by region at level 1, the built `list` and `density` forms) and
    bigBed (bigWig onto BBI, density without `measure` until RG); **SQ** 2bit and FASTA by region; **GG** GFF3 and GTF
    read whole (`assembleTranscripts`, `features`) and `line: 'gene'` (render) after GC1; **VC** VCF sites at level 1.
  - **RG** (stage 3, alone in the store) takes GF's protocol and keys part, as RR's draft (below, word for word): span
    rows and item extents in `tables.ts`, `measure`, `side` in keys, items at a tile size for every list shape.
  - **FR** (new, strong, after RG, before LA): the format lanes' by-region rest on RG's protocol: tabix levels from the
    index's bytes with exact coarse counts (BED, GFF, VCF), GFF's span rows by region, VCF's matrix and genotype
    classes, bigBed's `measure` and items at a tile size; then **PF** PAF, PIF and chain (pairs, `side`, a map's runs).
  - **LA after PF** (its liftover reads the chain runtime), so stage 3 is RG → FR → PF → LA; the GATE of stage 2 checks
    BZ to VC, stage 3's checks FR and PF.

## Real slices to record (each fetch asked first, Gate 0, 7)

Recorded as **pages** of the published file (core's 64 KiB pages the runtime read, keyed by file and page; a read the
recording did not make fails naming its page), so the container written by the real tools is what the test reads; a
change to a runtime's reads is a re-recording. Each with its `sources.json` row and licence, checked when asked. Around
TP53 (GRCh38 chr17:7,668,402-7,687,550), the p53 dataset's gene.

| Format | File to fetch | Pages | Licence to check |
|---|---|---|---|
| VCF sites | NCBI ClinVar `clinvar_20240107.vcf.gz` and `.tbi` (GRCh38 archive) | the `.tbi` whole (about 520 kB), the header, TP53's chunks | NCBI, public domain |
| VCF genotypes | 1000 Genomes 30x on GRCh38, chr17's phased panel `.vcf.gz` and `.tbi` (IGSR) | the header, a 10 kb region's chunks (about 300 records × 2,504 samples) | IGSR's terms |
| bigBed | UCSC hg38 GENCODE `.bb` (the `knownGene` track) | header, chromosome tree, TP53's R-tree path and blocks, the `name` extraIndex path for `TP53` | GENCODE and UCSC's terms |
| 2bit | UCSC `hg38.2bit` | header, index, chr17's record header, the mask pages TP53's search reads, its bases | UCSC (GRC sequence) |
| FASTA, bgzipped | Ensembl `dna_index` GRCh38 toplevel `.fa.gz` with `.fai` and `.gzi` | the indexes whole, TP53's blocks | Ensembl, no restrictions |
| GFF3 | Ensembl GRCh38 `.chr.gff3.gz` | TP53's gene lines cut by the script (a small whole file, `##sequence-region` kept) | Ensembl |
| GTF | GENCODE `annotation.gtf.gz` | TP53's lines cut | GENCODE |
| chain | UCSC `hg19ToHg38.over.chain.gz` | the chains touching chr17 cut, and UCSC liftOver's TP53 coordinates recorded beside | UCSC |
| PAF | a published PAF between two small genomes (JBrowse's grape–peach demo the candidate, GA's peach and grape) | a whole small file or its region | to check; generated if it fails |
| BED | UCSC hg38 `cytoBandIdeo.txt.gz` (already Gate 0, 6's ask) | whole (a few kB) | UCSC |

The repository already holds the current Seq\*'s recordings of hg38's 2bit, GENCODE bigBed and ClinVar (by its own
reads); they cannot serve page reads as recorded, so they are not reused (Gaps, 16).

## Where this differs from GA's and RR's drafts, and the text they take word for word

- **For RR's draft, Pairs (after "Pairs are tiled in 1D, by either end")**:
  > **A pairs source between two genomes is sided.** Its coverage marks each system's side (`SystemInfo.side: 'a' |
  > 'b'`), so `chr1` may be listed once a side; a request on it names the side it reads (`RegionRequest.side`), and the
  > store's key carries it in its system part (`a:chr1`, escaped as built). Side `a` is the side the alignment maps
  > from (PAF's query, a chain's `t`). A `pairs` node's `a` and `b` map each side's names through GA's names rule. A
  > source over one genome (BEDPE, self-synteny, Hi-C) is unsided. Rows sharing `row` with other starts are the blocks
  > of one alignment (its pieces): a hit names the alignment.
- **For RR's draft, "Items are read at a tile size, not a bin"**, its first sentence replaced by:
  > **Items are read at a tile size, not a bin, for every list shape** (`pairs`, `ranges`, `matrix`): the `list` form
  > (and a matrix's `items`) may be asked at any listed level, the level then only the tile's width.
- **For RR's draft, the Tables and the store** (a new bullet under The answer shapes):
  > **A ranges row may carry its item's whole extent** (`RangesList.itemStart`, `itemEnd`) when rows are pieces; a
  > **span row** (`start = end`, with an extent) is an item with no piece in the tile, drawn as its connecting line and
  > hit as the item. `sliceTable` emits a span row for an item whose pieces lie outside the slice, so a region read and
  > the whole read clipped agree. `RangesDensity.measure` (`items`, the default, or `coverage`) says what `count`
  > counts.
- **For RR's draft, the matrix's `codes`** ("genotypes `0`, `1`, `2`, 255 missing") replaced by:
  > `codes` (uint8: a genotype class of the item's allele: `0` no copy, `1` some but not every copy, `2` every copy, a
  > haploid call included, `3` partly called, `255` uncalled; the copy count is a `value` channel, `field: 'GT'`)
- **For GA's draft, An assembly** (after `alias_source`'s bullet):
  > **`sizes`** (a ref: a source whose coverage lists the sequences' lengths, a chrom.sizes or a `.fai`): read before
  > `source`, its lengths winning; absent, `source`'s. A derived sequence with no length from either makes the
  > assembly unavailable (`seq.resolve.assembly.lengths`, error, naming the sequences and `sizes`), never drawn with a
  > guessed length. Why: a large 2bit lists names without lengths past its read bound, and UCSC ships its chrom.sizes
  > beside it.
- **GA's map between assemblies** is read by region; its data here is a pairs source's equal-length rows
  (`map.source`); the projection by region stays GA's (GC1) and the figure LX's. No text change asked.

## Gaps for the author

Each a yes or no; the recommendation first.

1. **BGZF and tabix (TBI and CSI) in `@molstar/io/common/bgzf` and `common/tabix`** with the surface above (lines as
   bytes, extents by the index's rule), and `seqx-core`'s io allowlist widened to them: yes.
2. **Core's `ungzip` fixed for multi-member streams by splitting members** (verified failing: the ISIZE of the EOF
   block, the platforms' refusal of trailing members; BZ's test first): yes.
3. **One opt-in entry a format family**, `fasta` and `gff3` taken from the defaults by `withAdapter` (the defaults
   unchanged): yes.
4. **BED's extra columns named by the built `fields` with a column number** (`{ stain: '4' }`), no new param: yes.
5. **`source.features`** (`transcripts` by default, `genes`, `rows`), a line's role by its type (`transcriptTypes`), an
   intron-only tile giving a span row; a GTF without transcript lines exact only when read whole: yes.
6. **The ranges layer's `line: 'gene'`** (thick CDS, thin UTR, intron hairline with chevrons), built in GG: yes.
7. **No dynamic transcript nodes with codon maps from a GFF** this slice (transcripts are items): yes.
8. **Coarse counts of a tabix source exact or absent**, its levels listed per file from its index's bytes (a level
   only when every tile fits `maxRegionBytes`), no estimate drawn: yes.
9. **bigBed's zoom levels as coverage density** (`measure: 'coverage'`): yes.
10. **Sided pairs sources** (`side` on systems and requests; side `a` the map's from side: PAF's query, chain's `t`),
    amending RR: yes.
11. **GF's language additions ride GC1's bump** (`features`, `samples`, `block_index`, GA's `sizes`, `map.source`,
    `line: 'gene'`), so the language changes once: yes.
12. **`findNames` as an optional capability** (read-whole sources and bigBed's `extraIndex`; tabix sources none, trix
    later): yes.
13. **`pif`** (a PAF indexed from both ends through tabix) in PF: yes; else large PAFs only whole within `limits`.
14. **The plan amended**: GC1 declares GF's and GA's params; stage 2 runtimes only (BZ, BD, SQ, GG, VC sites); RG takes
    the protocol and keys part; a new lane **FR** after RG for the by-region rest; PF after FR; LA after PF: yes.
15. **Every 2bit length read at open through 4 KiB pages within 16 MiB**, else none listed and GA's `assembly.sizes`
    (a chrom.sizes, `chrom-sizes` in the BED entry) required: yes.
16. **Real slices as recorded pages** of the files listed, each fetch asked, the current Seq\*'s recordings not reused:
    yes.
17. **bgzipped FASTA with `.gzi`** (`block_index`) in SQ: yes, cut first if SQ runs long.
18. **The anchor base kept** in a VCF allele's extent (the index's rule; igv trims): yes.
19. **Genotype classes** (`0`, `1`, `2` with a haploid call, `3` partly called, `255`) in RR's matrix codes, the copy
    count a `field: 'GT'` value: yes.
20. **VCF levels per file from its index's bytes** (a dense file read in 16 kb tiles, no `levels` per source to
    write): yes.

## Review, unresolved

The independent review (2026-11-21) found three blocking issues (GFF by region inside an intron, 2bit lengths past a
bound against GA's draft remedies, the lane order against the plan's rules) and the rest to fix; the revision took
them all. Still open after it, for the gate:
- A GTF without transcript lines read by region has partial extents (said; exact when read whole).
- A dense VCF's sites past its listed levels draw "zoom in" (a bigBed is the remedy), and a source over several
  files (one a chromosome) is not designed.
- Phasing, breakends as pairs, name search over tabix sources (trix), the anchor base, `itemRgb` as a colour.
- PIF's `I`/`D` swap on `t` lines is taken from JBrowse's writer, verified in PF before it is relied on.
- RR's and GA's drafts take the texts above word for word, or the gate reconciles them.

## Review, unresolved

The independent review (2026-11-21) found three blocking issues (GFF read by region inside an intron, 2bit lengths past
4,096 sequences against GA's `minLength`, the lane order against the plan's rules) and eight to fix; the revision took
them all. The text RR's and GA's drafts must take (section "Where this differs…") is carried to them at the gate. Still
open after it:
- A GTF without transcript lines has exact extents only when read whole.
- Dense VCFs show "zoom in" past their levels.
- One source over several files (1000 Genomes and gnomAD ship one per chromosome) is not designed.
- Phasing, breakends as pairs, name search over tabix sources (trix), the VCF anchor base, `itemRgb`.
- PIF's swapped I/D is to be verified in PF.
