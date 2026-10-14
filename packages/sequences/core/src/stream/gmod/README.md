# @molstar/seq-core/stream/gmod

Stream providers (seqspec:sources#streamed-sources, seqspec:sources#the-viewers-side) for indexed genome files, over byte-range reads. Four readers,
one per format; each takes a **transport** and returns a `StreamProvider`:

| Reader | File | `data_kind` | Answers |
|---|---|---|---|
| `twoBitProvider(read, { only?, maxSequences? })` | `.2bit` | `letters` | the manifest (a system per sequence), `letters` patches (level -1) |
| `bigWigProvider(read)` | `.bw` | `values` | the manifest (`levels`: the file's zoom reductions; `domain`: its range), `bins` (`min`, `max`, `mean`) per level, `items` (a dense `value` column) at level -1 |
| `bigBedProvider(read, { maxItems?, label?, category?, vocabulary? })` | `.bb` | `items` | the manifest, `items` (a row per BED item, per block of a BED12; `label` and `category` from the fields the host names), `bins` (`count`; coarse levels from the file's zoom records) |
| `vcfProvider(read, tbi, { info?, lengths?, maxItems? })` | `.vcf.gz` + `.tbi` | `items` | as BigBed, a row per ALT allele |

`gmodStream(open)` is a host's `loaders.stream`: `open(source)` maps a source to a `GmodFile` (or `null`: a 2bit's
`only` keeps the sequences a source needs; a VCF's `index` is its `.tbi`, as bytes or a function fetching them), and
`httpFile(format, url)` is one over HTTP ranges. The name mapping of a live source is the host's: `renamed(provider,
name)` gives a file's coordinate systems the document's names (`name(fileName)`, null leaving one out; requests mapped
back, a name the file lacks `StreamNotServed`), and `chrNames` is UCSC's for Ensembl and NCBI names (`1` → `chr1`,
`MT` → `chrM`). Coordinates need no mapping: every reader answers 0-based, half-open (a VCF's 1-based positions and a
BED's are converted where read). A VCF's `lengths` come from the host (ClinVar's header has no `##contig` lengths).

## The transport

`read(offset, length, signal) → Uint8Array`: `length` bytes from `offset`, fewer where the file ends, none past it.
`fetchRanges(url)` sends `Range: bytes=a-b` and requires 206 (a server that ignores `Range` answers 200 with the
whole file: an error here, not a download). Tests read recorded bytes (`_test/slices.ts`). A host may cache ranges or
batch them: the readers keep their own small caches (the opened header, a 2bit sequence's record, 8 kb pages for its
binary searches), nothing per request.

A reader opens its file at its first request, which is the manifest; a failed open is not kept.

**Bounds.** A length or count a header gives is checked before it is read, and a refusal is the request's error, as
an `.svsx`'s bounds are (`seqspec:sources#archives`): a BigWig's or BigBed's chromosome tree and R-tree at most 32
levels deep and 65,536 items a node, its chromosome tree at most 64 MB (about 2M contigs) and its autoSql 1 MB, a
block inflating to at most the header's `uncompressBufSize`; a BGZF block at most 64 KB, its ISIZE checked, and a
run (a `.tbi`) at most 256 MB; a 2bit's index read only as far as the file holds it, and at most a million
sequences kept (after `only`; `maxSequences`).

## Items, levels and the zoom gate (BigBed, VCF)

- **Levels never read every item of a chromosome** (ST4b, David: "don't they have some level of detail?"). A BigBed
  lists one level counted from items (`256` positions a bin, regions of 64 kb, cheap to read), then its file's zoom
  levels: their records' summed coverage over a bin's width is its count (the items' mean depth there; a zoom record
  covers bases, not items), read through the zoom index, never the data. A VCF lists `[256, 4096, 65536]`
  (`ITEM_LEVELS`) and answers every level from its tabix index alone: a 16 kb window's compressed bytes (from the
  linear index) over a record's compressed bytes, measured from the first blocks past the header that are read
  (the header's own read, else a region's; `BYTES_PER_ITEM`, 28, as measured on ClinVar, until then), spread over the bins; the patch says
  `estimate: true`. An items source serves levels through `ItemsSource.bins` (null: count them from the items) and
  lists them through `levels`.
- **The zoom gate:** a data request (level -1) whose region holds more than `maxItems` (2,000) items, or more bytes
  than its `max_bytes` (a VCF decides it from its index before reading anything, `ItemsSource.estimate`, as JBrowse 2
  does; `force: true` reads it whole), is answered with
  `bins` at the finest level whose bins fit in a tile; the patch's `level` is that level, not -1, and its `start`
  and `stop` are the request's. The viewer (ST4b) holds the region at that level and asks for data again only
  for a smaller region.
- **Rows.** `start`, `stop` (0-based, half-open), `id`, `label`, then the format's fields as text. An `id` names an
  item: the viewer joins rows of one `id` into a feature, so a BED12 gives a row per block under one id (exons), and
  an id is `name@start-end` (BED names repeat) or `pos:ref>alt` (a VCF `ID` repeats over ALT alleles).
- An item crossing a region's edge is in the region whole, in every region it meets.

## 2bit

- **Byte order.** The signature `0x1A412743` read little-endian says little-endian; `0x4327411A` big. Every field
  follows. Version 1 (64-bit offsets, files over 4 GB) is refused by name.
- **Layout.** A header (signature, version, count, reserved), an index of (name length, name, offset) per
  sequence, then per sequence at its offset: `dnaSize`, the N blocks (count, starts, sizes), the mask blocks (count,
  starts, sizes), a reserved word, and `ceil(dnaSize/4)` bytes of bases. **A sequence's size is in its record**, not
  in the index: the manifest reads 4 bytes per sequence (455 for hg38: `only` keeps the ones a source needs).
- **Bases** are 2 bits, four to a byte, the first base in the high bits; the codes are `T=0, C=1, A=2, G=3`.
- **N blocks** are runs of `N` the packed bases know nothing of (they hold `T`s there); **mask blocks** are soft
  masking (repeats), lower case. Both lists are sorted by start. chr1 of hg38 has 331,893 mask blocks (2.7 MB):
  a region finds its blocks by binary search over the on-disk list (a handful of 8 kb pages), never reads the list.
  N is applied first and the mask over it (`twoBitToFa` writes `n` for a masked N).
- **0-based, half-open**, as the port's. The viewer holds the letters upper-case; the case is kept in the patch.
- A region past the sequence's end is cut at it (the patch has fewer letters than `stop - start`).
- Levels: none (`levels: []`); a level request is `StreamNotServed`.

## BigWig and BigBed

One container (`bbi.ts`), two data formats.

- **Byte order** by the magic (`0x888FFC26` BigWig, `0x8789F2EB` BigBed), whichever way it reads: files written on
  big-endian machines exist. Offsets are 64-bit (read as numbers; files stay under 2^53).
- **Header** (64 bytes), then the zoom headers (24 bytes each, at most 10 in practice; 16 are read at once): reduction
  level, data offset, index offset. Then optionally the autoSql (BigBed), the total summary (40 bytes: bases covered,
  min, max, sum, sum of squares: the BigWig's `domain`), the **chromosome B+ tree** (name to id and size; the
  keys are NUL-padded, and the whole tree lies between its offset and the data's, so it is read once), the data, the
  data's **R-tree index** and each zoom level's data and index.
- **Finding data.** The R-tree is read top-down: each node is a read (`4 + blockSize * 32` bytes is the most it
  can be; leaves hold 32-byte items (chrom range, base range, offset, size), inner nodes 24 bytes), its matching
  children in parallel. A leaf item's range is conservative: filter what a block holds by position.
- **Blocks** are zlib-compressed when the header's `uncompressBufSize` is not 0 (zlib, not gzip: `DecompressionStream("deflate")`).
- **BigWig data** is a section per block: a 24-byte header (chrom id, start, end, step, span, type, item count),
  then `bedGraph` (start, end, value), `variableStep` (start, value; span from the header) or `fixedStep` (value;
  start from the header and step) items, values 32-bit floats.
- **Zoom level choice.** The manifest's `levels` are the file's own zoom reductions, finest first, so a level
  request is the zoom level of the same index: the viewer asks for the coarsest whose bin is at most a device
  column, and the data (level -1) where finer than all. A file with no zoom levels serves data only.
- **Zoom records** (32 bytes: chrom id, start, end, valid bases, min, max, sum, sum of squares) are **not on a grid**:
  a record starts where its data starts and spans up to the reduction. A request's bins are on its grid, so a
  record is spread over the bins it meets by its share in each (valid bases and sum prorated: `mean` is sum over
  valid bases; `min` and `max` are of every record a bin meets). So a bin's `min` and `max` may be wider than its own
  data, never narrower, and the means keep the total.
- **BigBed items** are (chrom id, start, end, rest) with the rest one NUL-terminated, tab-separated string: the BED
  fields after the first three, named by the autoSql (`bigGenePred` has 12 BED fields and 15 more) or, without it,
  by BED position. The BED12 block fields (count, sizes, starts relative to the item's start) become rows.
- **0-based, half-open** in the file, as the port's. BigBed has no per-item summary; see Levels above.

## Tabix and BGZF

- **BGZF** is a series of gzip members of at most 64 kb each, with an extra field `BC` holding the member's total size
  minus 1: a block is found by its size, so a run of blocks is split by reading each header (the members are
  inflated one at a time with `DecompressionStream("gzip")`; one multi-member stream is not portable). A last empty
  block marks the end. A read that cuts a block loses that block only.
- **Virtual offsets** are 64-bit: the block's file offset shifted left 16 bits, ORed with the offset into the block's
  inflated bytes. A chunk is (begin, end) virtual offsets: it ends `uend` bytes into the block at `cend`, so that
  block is read too (its size is in its 18-byte header: a second, tiny read).
- **The `.tbi`** is itself BGZF, then little-endian: magic `TBI\1`, the reference count, the format (2: VCF; flag
  `0x10000` marks 0-based), the columns (sequence, begin, end), the meta character, skip count, the names (NUL
  separated), then per reference its bins (bin number, chunks) and its **linear index** (the smallest virtual
  offset of any record meeting each 16 kb window). It is read whole, given as bytes or a function fetching them
  under the open's signal (`GmodFile.index`; `httpFile` fetches `${url}.tbi`).
- **A region's chunks.** The bins that can meet `[beg, end)` are the SAM specification's `reg2bins` (the levels
  of 512 Mb, 64 Mb, 8 Mb, 1 Mb, 128 kb and 16 kb bins, numbered from 0, 1, 9, 73, 585 and 4681); take their chunks,
  drop those ending before the linear index's offset for `beg >> 14`, sort and merge. Chunks hold records of other
  positions: filter every line by chromosome and overlap.
- **VCF is 1-based, closed.** A record is `[POS - 1, POS - 1 + len(REF))`, or to `INFO/END` (structural variants).
  The anchor base of an indel (VCF's leading shared base) is kept, not trimmed.
- **Names.** The index's names are the VCF's: ClinVar and Ensembl say `1`, UCSC says `chr1`. A source's
  coordinate system key is the index's; the host maps one to the other (ST4b).
- **Lengths** are not in a tabix index. The manifest takes `lengths` (the host's, from a 2bit or chrom.sizes), else
  the header's `##contig=<ID=1,length=...>` lines (ClinVar has none), else where the linear index ends (a multiple of
  16 kb below the end).
- **Columns:** `ref`, `alt` (one row per ALT allele), `name` (the VCF `ID`, empty for `.`), `category`
  (`SNV`, `MNV`, `insertion`, `deletion`, `structural`, by lengths and a `<` alt) and the INFO keys named in `info`.

## Recorded slices (`_test/fixtures`)

Tests read small byte ranges of real hg38 files, recorded once by `scripts/record.ts` (by hand; the network):

| Key | File | Recorded |
|---|---|---|
| `twobit` | `goldenPath/hg38/bigZips/hg38.2bit` (835 MB) | 27 ranges |
| `bigwig` | `gbdb/hg38/bbi/gc5BaseBw/gc5Base.bw` (1.76 GB), GC percent in 5-base windows | 28 ranges |
| `bigbed` | `gbdb/hg38/gencode/gencodeV50.bb` (57.5 MB), the `knownGene` track | 10 ranges |
| `vcf`, `tbi` | NCBI `clinvar_20240107.vcf.gz` (91 MB, an archive, so it does not change) and its `.tbi` (520 kb) | 5 ranges, the whole index |

`slices.json` lists each file's ranges and `<key>.bin` holds the bytes; a read the recording did not make fails
naming its range. The ranges are exactly what the readers read, so a change to a reader's reads is a re-recording
(`node scripts/record.ts [twobit|bigwig|bigbed|vcf]`). `expected.json` holds the answers of sources that do not use the
readers: UCSC's REST API (`api.genome.ucsc.edu/getData`) for letters, GC windows and genes, and for the VCF the first
and last megabytes of the file inflated with Node's gunzip and filtered by line.

A small synthetic genome in all four formats, whole (`_test/stream/gmod/demo-files.ts`: two sequences, a 2bit with an
N run and soft-masked repeats, a BigWig of GC percent, a BigBed of BED12 genes, a bgzipped VCF and its `.tbi`), is made
in memory, so any region of it reads: the `Seq*/File formats` story and `demo-files.test.ts` read it.
