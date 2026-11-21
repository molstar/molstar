# BGZF and tabix

`@molstar/io/common/bgzf` and `@molstar/io/common/tabix`: blocked gzip read whole or by virtual offset, and TBI and CSI
indexes of a BGZF text file sorted by position, a region's lines as bytes. Lines and offsets, never tables: a Mol\*
reader (a bgzipped trajectory index, a tabix-indexed annotation) uses them as a genome runtime does. The implementation
is in `src/common/internal/` (not published): `bgzf-blocks.ts` (the block walk, the reader) and `tabix-bins.ts`
(`reg2bin`, `reg2bins`); tests reach them there.

## 1. BGZF

- **A block is a gzip member** with FEXTRA and a `BC` subfield of 2 bytes holding BSIZE (the block's size minus 1; SAM
  spec §4.1). A block that is not one (no gzip magic, no `BC`, BSIZE smaller than its header and trailer, ISIZE over 64
  KiB) is `malformed`; one that runs past the bytes is cut short, `malformed`. `isBgzf` tells it by the first block's
  header, never by a name.
- **Each block is inflated by core** (`inflateStream(…, 'gzip', buf)`, the platform checking CRC-32 and ISIZE; core's JS
  inflater and a CRC-32 where the platform has none), into a buffer of its ISIZE: never a second inflater, never a
  buffer grown. A block whose bytes are not its ISIZE or CRC-32 is `malformed`.
- **A virtual offset** is a block's file offset × 2^16 + an offset into its inflated bytes, held as a number: files
  under 2^37 bytes (128 GiB); a uint64 past 2^53 is refused. An end at its block's ISIZE reads as the next block's start
  does (htslib writes either).
- **Read whole** (`inflateBgzf`): the blocks walked by BSIZE first, their ISIZEs summed against `maxInflatedBytes`
  before the output is allocated.
- **Read by range** (`bgzfReader(ranges)`): over core's `AssetRangeReader` (`openRanges`, 64 KiB pages and their cache:
  no page cache here). A read's compressed span (from its first block to its last) is bounded by `maxInflatedBytes`
  (plus a block) before it is read, so an index's offsets never size a read: a block that does not shrink is stored, so
  compressed bytes are at most about the inflated ones. Then the blocks' headers are walked, the exact length summed
  against the bound, and the blocks inflated; a virtual offset off the blocks or past its block's ISIZE is `malformed`.
  It keeps `cacheBlocks` inflated blocks (16) for a pan back; `dispose` lets them go, and the range reader stays its
  owner's.
- **The EOF block** (28 bytes) is checked by `hasEof` with one read of the file's tail; a missing one is the caller's
  warning (the file may be cut short), the file still read (htslib's rule).

## 2. Tabix: the index

- **One reader for TBI and CSI**, told by the magic (`TBI\1`, `CSI\1`) of the index inflated (BGZF) or as given. A TBI
  is `minShift` 14 and `depth` 5 (positions under 2^29) with a linear index; a CSI has any scheme (`minShift` 1–32,
  `depth` 1–10, `minShift + 3 × depth ≤ 52`), a bin's `loffset` in place of the linear index, and the tabix header in
  its aux data (a CSI without one is not a text file's index: `malformed`). Why both: a plant chromosome passes 512 Mb.
- **The header**: format (`generic`, `sam`, `vcf`; the `0x10000` flag: positions 0-based, ends exclusive), the 1-based
  sequence, begin and end columns (end 0: none), the meta character, lines skipped, the names (one a reference).
- **Every count and offset is checked before it allocates**: the index inflated within `maxIndexBytes`; each count
  (references, bins, chunks, linear windows) against the bytes left at its item's least size and against the scheme
  (bins: `binCount(depth)` and the metadata pseudo-bin, `binCount + 1`, 37,450 in a TBI); a bin listed twice; a chunk
  `begin ≤ end` (but the pseudo-bin's counts); names one a reference. A failure is `TabixError` `malformed`, a bound
  `too-large`; fuzzed and truncated indexes throw nothing else.
- The bins are kept off the public `TabixIndex` (a module's weak map): `regionChunks` takes an index `readTabixIndex`
  read.

## 3. Tabix: a region

- **Chunks** (`regionChunks`, 0-based half-open): the bins `reg2bins` gives (htslib's, by division so a CSI may pass
  2^31), each chunk ending past the region's least offset (TBI: the linear index at the region's window, a zero entry
  taking the last nonzero before it; CSI: the `loffset` of the smallest bin holding the start, up the tree), sorted,
  overlapping chunks and chunks meeting in one block merged (each line read once).
- **Cost** (`regionBytes`): each merged chunk from its first block's start to its last block's start, as the index's
  offsets say (the last block's own bytes, at most 64 KiB a chunk, are not in the index). A read over `maxRegionBytes`
  is `too-large` before a byte is read; the chunks inflated together stay within `maxInflatedBytes`.
- **Lines as bytes** (io's rule: text is tokenized as bytes): each line's span (a `\r` before its `\n` dropped) and its
  extent, 0-based half-open, by htslib's `tbx_parse1`: generic `begin` (− 1 unless 0-based) to the end column, or
  `begin + 1` without one; SAM `POS − 1` plus the CIGAR's reference length (M, D, N, =, X); VCF `POS − 1` to `INFO/END`
  when INFO starts with `END=` or holds `;END=` and the value is past `POS − 1` (whatever the alleles), else
  `POS − 1 + len(REF)`. So `readRegion` answers exactly what the index would; a line outside the region from a chunk's
  edge is dropped there, and a region read equals the whole read clipped.
- **A line naming another sequence** (an index for another file), without a position, or past 2^31 − 1 is `malformed`
  for that region, never a row in the wrong place; a line over `maxLineBytes` is `malformed`.
- **Header lines** (`readHeaderLines`): the bytes before the first indexed line (or, with none, blocks from the start
  until a line that is not a header), within `maxHeaderBytes`; the lines starting with `meta`, and the first `skip`.
  Their extents are −1.
- **Positions stay int32** (the protocol's): a sequence longer than 2^31 − 1 is the caller's to refuse at open.

## 4. Bounds

| Bound                                                             | Default | Over it                   |
| ----------------------------------------------------------------- | ------- | ------------------------- |
| a block, compressed and inflated                                  | 65,536  | `malformed`               |
| a read's inflated bytes, its blocks together (`maxInflatedBytes`) | 64 MiB  | `too-large`               |
| inflated blocks kept (`cacheBlocks`)                              | 16      | the least recent let go   |
| an index inflated (`maxIndexBytes`)                               | 64 MiB  | `too-large`               |
| a region's chunks, compressed (`maxRegionBytes`)                  | 8 MiB   | `too-large`, nothing read |
| one line (`maxLineBytes`)                                         | 16 MiB  | `malformed`               |
| the header lines (`maxHeaderBytes`)                               | 16 MiB  | `too-large`               |

`BgzfError` and `TabixError` carry `kind` (`malformed`, `too-large`) and are told by `name` across copies of the
package. A range read's own failure (`AssetRangeError`) passes through: the region fails, never with wrong lines.
