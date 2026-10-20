# Vendored: biomsa

`src/biomsa/` is biomsa's source, vendored so the aligner is deterministic, erasable TypeScript and free of dependencies
(A1).

- **Origin**: <https://github.com/ppillot/biomsalign>, commit `6561931b597a04e7a5358bfbbbe8afd97ecfc2b1` (npm `biomsa`
  0.3.3).
- **License**: MIT, Copyright (c) 2021 Paul Pillot; the text in `src/biomsa/LICENSE`.
- **Copied**: `src/` without `src/**/__tests__/` and `src/dotplot/` (an unused, commented-out dotplot); `TODO.md` kept
  as it was.
- **Used through**: `src/biomsa-aligner.ts` (the `biomsa` aligner), which calls `BioMSA.align` (the named export of
  `src/biomsa/index.ts`; upstream's was the default) on two sequences and reads the gapped rows as a CIGAR, and
  `src/aligner.ts`, which uses `getSequenceType` to detect the sequence type.

## Changes

Each change in the code carries a `STORIES-UNTOLD(Bn)` comment, except B5, B6 and B18, which are mechanical and
everywhere. Nothing else was changed; the code keeps biomsa's names, comments and logic. B7 to B11 (A1 of the Seq*
rebuild; B11 from A1's review) and B12 to B14 (A1b) fix upstream bugs and move the results away from upstream's: the
tests in `_test/aligners.test.ts` pin them.

- **B1, `const enum`s as `as const` objects** (erasable TypeScript, Node's type stripping): `SEQUENCE_TYPE`
  (`align/params.ts`), `TRACE_BACK` and `ALIGNOPT` (`align/align.ts`), `EPATH_2_STRING` (`utils/estring.ts`),
  `NODE_TYPE` (`sequence/tree.ts`), each with a type alias of the same name for its values, the members keeping their
  numbers. `BioMSAClass.typeSeq` (`index.ts`) is typed `SEQUENCE_TYPE`, since a member of an `as const` object is a
  literal type where the enum's member widened to the enum.
- **B2, typed-array declarations** (`align/align.ts`, `MSAMSAAlignment`): `lProfAAAScores` and `lResList` typed
  `Float32Array` and `Uint8Array`, since TypeScript 5.7+ infers `Float32Array<ArrayBuffer>` from the constructor, which
  the later `subarray` assignments (`<ArrayBufferLike>`) do not satisfy. Types only.
- **B3, ambiguous letters encoded deterministically** (`sequence/sequence.ts`, `setRandomCode`, name kept): biomsa
  replaced a letter outside the standard alphabet by a random choice among the residues it may stand for
  (`Math.random()`), so the same input could align differently. Now it is always the first of those choices in the
  encoding's order (A C G T; A C D E F G H I K L M N P Q R S T V W Y):
  - nucleotides: M, R, W, V, H, D → A; S, Y, B → C; K → G; N and any other letter → A;
  - amino acids: B (Asx) → D; Z (Glx) → E; J (Xle) → I; U, X and any other letter → A. An ambiguous letter therefore
    scores as that residue (two X's score as two A's).
- **B4, the logger is a no-op** (`utils/logger.ts`): `start`, `add` and `summary` keep their signatures and do nothing;
  the original timed events with `performance.now()` and printed them with `console.table`. Its calls are behind `DEBUG`
  (false) anyway. biomsa's other console calls are left as they were: those behind `DEBUG`, and two on error paths of
  multiple alignment (`console.error`/`console.dir` before a throw in `align/noalign.ts`'s `centerStarNoAlign`,
  `console.trace` in `sequence/sequence.ts`'s distance matrix), which a pairwise alignment does not reach.
- **B5, imports for NodeNext**: relative imports end in `.js` (`".."` became `"../index.js"`), as the rest of Mol\*'s v6
  sources, resolved to the `.ts` beside them by the bundlers, Vitest and tsc; type-only imports are `import type`
  (TypeScript's `verbatimModuleSyntax`; applied by Biome's `useImportType` fix).
- **B6, formatting**: Biome's formatter with the repository's settings (`biome.json`: 2 spaces, single quotes, 120
  columns), commit "Mol*'s formatting"; the imports are not reorganized. No lint rule is turned off for `src/biomsa/**`
  (B18).
- **B18, Mol\*'s lint** (the port into the Mol\* repository; mechanical, no behaviour change): `let` that is never
  reassigned became `const`; `==` and `!=` on operands of the same type became `===` and `!==`; the default exports
  became named ones (`BioMSA` in `index.ts`, `Log` in `utils/logger.ts`, `BitArray` in `utils/bitarray.ts`), with their
  importers; the unused `score1` (the result of `treeAlign` in `align/progressive.alignment.ts`) was dropped, the call
  kept; the two `throw`s of a string in `index.ts` stay, each with a `biome-ignore lint/style/useThrowOnlyError` and the
  reason (callers see the string).
- **B7, a leading gap in `b` opens like one in `a`** (`align/align.ts`, `pairwiseAlignment`, the start of each row of
  the DP): `lPrevMatch = GAP_OPEN + GAP_START_CORRECTION` where upstream had `GAP_START_CORRECTION` alone. Row 0 of the
  DP starts at `GAP_OPEN + GAP_START_CORRECTION` (-150 by default), column 0 started at `GAP_START_CORRECTION` (+150):
  leaving the first letters of `a` unpaired was rewarded, so `align("KV", "MKV")` gave `1I1M2D` and `align("M", "MKV")`
  gave `1I3D`, pairing nothing. Now `1D2M` and `1M2D`; `align("MKV", "KV")` is `1I2M`. The default result for the DNA
  pair in the tests changes from `4I11M` to `6M4I5M` (`gap_open: 800` gives `4I11M`).
  - **The same asymmetry in the sequence-to-profile DP** (`MSASeqAlignment`) is B15. The profile-to-profile DP
    (`MSAMSAAlignment`) starts both from the profiles' own gap scores and is symmetric.
  - **Not fixed, protein k-mers alias**: the diagonal heuristic (`extractMinimizers`) packs a k-mer at 2 bits a letter
    into 16 bits, right for nucleotides (4 codes) and too small for amino acids (20 codes), whose k-mers collide. Seeds
    for long proteins are weaker than for nucleotides; the alignment is still a valid global one.
- **B8, positions in `Int32Array`** (`align/noalign.ts`, `TMinzStore`, `extractMinimizers`, `extractLongestPath`):
  `kmerPos`, `winPos` and `winPosEnd` were `Uint16Array` and wrapped past 65,535 letters (a 2 kb fragment at offset
  66,000 of a 70 kb sequence was placed near 464, without an error); `lTraceBack` of the diagonal chain was `Int16Array`
  and wrapped past 32,767 diagonals. All `Int32Array`. The k-mer values (`kmer`, the k-mer array) stay `Uint16Array`: 16
  bits is what an 8-mer of nucleotides takes.
- **B9, an empty minimizer store for a sequence shorter than k** (`align/noalign.ts`, `extractMinimizers`): a sequence
  under 8 letters made `new Uint16Array(length - 8 + 1)` throw "Invalid typed array length" (`method: auto` with a
  sequence over 1,600 letters against one under 8, or `method: diag` on any short pair). It returns an empty store, so
  `noalignPair` finds no diagonal and takes the existing full Needleman–Wunsch fallback.
- **B10, a size guard in `pairwiseAlignment`** (`align/align.ts`): `checkCells("biomsa", n, m)` (our `src/limits.ts`)
  before the traceback is allocated throws `AlignerSizeError` when (n + 1) × (m + 1) is over `MAX_CELLS` (2^26). It
  covers `method: complete`, the fallback of the diagonal method, and the gaps between diagonals. The multiple-alignment
  DPs (`MSASeqAlignment`, `MSAMSAAlignment`) are B16's.
- **B11, the minimizer store sized to the number of windows** (A1's review; `align/noalign.ts`, `extractMinimizers`):
  upstream sized it to twice a random sequence's minimizer density, 4 / (w + 1) of the length. A tandem repeat has a new
  minimizer every few windows (an `AC` repeat about one per two letters), so the store filled; the writes past its end
  were dropped while `count` grew, and every seed after the repeat read undefined positions: NaN diagonals, or none.
  `head + "AC" × 3000 + tail` against `head + tail` (10 kb and 5 kb) gave `15000M6000D`; now `10000M6000D5000M`. A
  window adds at most one minimizer, so the store now holds one per window (14 bytes a letter).
- **B12, diagonals chained in their order in `a`** (A1b; `align/noalign.ts`, `noalignPair` and `extractLongestPath`).
  Upstream sorted the diagonals by `begin + diagId` (2 × the position in `a` − the position in `b`): a diagonal after an
  insertion in `b` longer than the shared stretch before it sorted first, and the chain could not reach it from the
  earlier ones (5 kb + 6 kb + 5 kb against 5 kb + 5 kb paired most of the first half with the insertion; 25 of 100
  random insertions of 1 to 10 kb in 20 kb came out split, some with an `AlignerSizeError`). They are sorted by `begin`,
  then `diagId`, so every possible predecessor of a diagonal comes before it; the garbage collection's bound
  (`pLenA - lDiag.begin`, what a path can still gain) holds only in this order. The chaining's other rules were checked
  against it; three changed, each needed for an insertion to come out as one run:
  - a diagonal extends the path it extends best (score after the penalty), not the first from the top that accepts it;
  - a jump to another diagonal costs `min(KSIZE, |Δ diagId|)` on top of upstream's tie-break: a chance k-mer match
    inside an insertion (a 7-letter range) adds less than the second jump costs, where it split the insertion in two
    (about one in eight); a one-letter indel still costs one, so dense short indels chain as before (the 33,000-diagonal
    test, 8-letter seeds a letter apart);
  - a diagonal may overlap the path's tip, in `a` or in `b`, by fewer letters than its length, the overlap not counted;
    `noalignPair` trims its start. Equal letters at an insertion's ends make the seeded ranges on either side overlap by
    a letter or two, and the diagonal after it could not be chained (`12426M` and then fragments). The rule that starts
    a new path from a diagonal below every path's bottom (`lDiag.diagId >= lMinDiag.diagId`) holds in the new order
    (such a diagonal cannot extend any path) and is kept; dropping its `diagId` condition changed none of about 1,800
    random insertion cases tried.
- **B13, the tie-break penalty uses the tip's diagonal** (`extractLongestPath`): `Math.abs(lDiag.diagId - lMaxDiagId)`
  subtracted the path's index in the list from a diagonal number; it is `lMaxDiag.diagId`. `extractLongestPath` is
  exported for the test.
- **B14, the ends aligned, not anchored at the corners** (`noalignPair`): upstream chained the first diagonal to an
  empty one at (0, 0) on the main diagonal and extended that over equal letters, so an insertion's first letters equal
  to the sequence's own gave `3M6000D4997M` for `6000D5000M`; the segment before the first diagonal was aligned with
  full-price end gaps, and the tails after the last were appended unaligned (an insertion 9 letters before the end gave
  `20000M6184D`, the last 9 letters paired with the insertion's first). Now the empty diagonals at (0, 0) and at the end
  of both sequences are not extended, an empty one at the end is added, and the leading and trailing segments are
  aligned with their outer gap favored, as the full DP does. The backward extension of a diagonal compares the letters
  before its start, so it reaches the previous diagonal's end when they are equal; upstream stopped a letter short
  (`begin = m + 1` at the bound), leaving a letter to the DP and, with the empty end diagonal, putting its start past
  the end. The trailing segment's DP is new: unseeded tails over `MAX_CELLS` (about 8 kb each) now throw
  `AlignerSizeError`, as an unseeded leading segment already did.

- **B15, a leading gap in the sequence opens like row 0's in the sequence-to-profile DP** (NP1b; `align/align.ts`,
  `MSASeqAlignment`, the start of each row of the DP):
  `lPrevMatch = profA.m_ScoreGapOpen[0] + profA.m_ScoreGapClose[i - 1] / GAP_START_FACTOR_A`, the profile's own open at
  its first position and close at this one, halved unless `DISABLE_FAVOR_START_GAP`, as `MSAMSAAlignment`'s column 0 has
  it. Upstream started the column at `GAP_START_CORRECTION_B` alone (+75 by default where row 0 started at -225),
  rewarding profile columns left unpaired at the start. Four sequences of one family, the last with a different start,
  came out with two columns of gaps before the last row's start and the others' first letters unpaired; now the gapless
  alignment (`_test/msa.test.ts`).
- **B16, a size guard in both profile DPs** (A1's review; `MSASeqAlignment`, `MSAMSAAlignment`):
  `checkCells("biomsa", profile length, sequence or profile length)` before the traceback (half a byte a cell) is
  allocated, as B10 does for `pairwiseAlignment`; over `MAX_CELLS` it throws `AlignerSizeError`, and the multiple
  alignment rejects with it.
- **B17, a multiple alignment that yields** (NP1b; `index.ts`, `align/progressive.alignment.ts`):
  `align(seqs, opt, step)` awaits `step` before each alignment of a node's two children, so the caller can yield to the
  event loop and throw to stop a long run (the `biomsa` aligner's `alignMany` does, about every 20 ms, on the stage's
  abort signal). `treeAlign` is async for it, and `align`'s promise executor became an async function whose throws
  reject. The singleton's state (`sequences`, the gap character) is read into locals before the first await, since the
  next call replaces it while this one waits. `centerStarNoAlign` (the diagonal method for sequences over 1,600 letters)
  does not yield.

**Our side of A1** (not in `src/biomsa/`): `src/limits.ts` (`MAX_CELLS`, `AlignerSizeError`, `checkCells`),
`AlignerOptionsError extends RangeError` in `src/aligner.ts` (thrown by `readOptions` and by the `biomsa` aligner's
matrix and type check; the only error a caller maps to a bad-options report), `Object.hasOwn` in `readOptions`, a
`matrix` that names the type when `type` is `auto`, `myersAlign` validating its input (letters only, `maxDistance` null
or a non-negative integer, equalities of two letters; `RangeError`), and `needleman_wunsch` checking `MAX_CELLS` before
it runs Mol\*'s shared DP (`@molstar/core/data/util/needleman-wunsch`, two rows of scores, a byte of traceback a cell).

# Ported: edlib (A2)

`src/myers.ts` is a TypeScript port of edlib's Myers bit-vector aligner, the `myers` aligner (A2). A port, not a copy:
the algorithm, its control flow and its names (`buildPeq`, `calculateBlock` as `advance`,
`myersCalcEditDistanceSemiGlobal` as `semiGlobal`, `myersCalcEditDistanceNW` as `nw`, `obtainAlignment`,
`obtainAlignmentTraceback` as `traceback`, `obtainAlignmentHirschberg` as `hirschberg`) follow edlib line by line.

- **Origin**: <https://github.com/Martinsos/edlib>, commit `0ddc23ea06abd95ab6227442acb4c27bd7607b68`,
  `edlib/src/edlib.cpp` (and the `edlib/include/edlib.h` it implements).
- **License**: MIT, Copyright (c) 2014 Martin Šošić; the notice is in the file's header.
- **Mapping**: edlib's query is `a`, its target `b`; its modes NW, SHW and HW are `global`, `prefix` and `infix`; its
  `k` is `max_distance` (`null`: edlib's `k = -1`, k doubling from a word until a solution is found); its additional
  equalities are `equalities`. Its INSERT (a letter only in the query) is our `I` and its DELETE (a letter only in the
  target) our `D` (seqspec:model#coordinate-maps); MATCH and MISMATCH are both `M`, as the other aligners give. In
  `prefix` and `infix` the target's letters outside the aligned part become leading and trailing `D`, so the CIGAR
  covers both sequences.

## Differences

- **32-bit words**: a block covers 32 query letters, not 64; words are int32 in `Int32Array`s, shifted with `<<` and
  read with `>>>`, so no `BigInt`. The high-bit test is `x < 0` (or `x >>> 31`). Int32 rather than `Uint32Array`, so the
  values stay small integers for the engine; the bits are the same. The unbounded search starts at k = 32, not 64.
- **One allocation per call**: `myersAlign` allocates every buffer once (the Peq tables, the band's blocks, the stored
  columns, Hirschberg's two score columns, the path) sized for the whole problem; the traceback's columns are capped at
  its 1 MB limit, and sub-problems use prefixes. edlib allocates per sub-problem. Sub-sequences are `subarray` views
  where edlib offsets pointers. The path is appended in place, where edlib allocates and concatenates halves.
- **The traceback's threshold** counts 32-bit words: (2 × 4 + 4) bytes per block and column plus 8 per column, under 1
  MB, so the split between traceback and Hirschberg falls at slightly different sizes than edlib's, and a path may
  differ from edlib's among equally optimal ones.
- **Positions**: `semiGlobal` keeps the first and last best end positions, the only ones edlib reads (the first for the
  path's end, the last of the reversed search for its start), not the whole list.
- **Smaller**: `getBlockCellValues`, `readBlock`, `readBlockReverse` and `allBlockCellsLarger` read cells in loops
  without the 64-element array; the Peq table drops edlib's wildcard row, which is never read; `transformSequences`
  upper-cases first (letters are case-insensitive) and an unknown sequence character is refused before (`alignable`).
- **Errors**: where edlib returns `EDLIB_STATUS_ERROR` from Hirschberg's split (no optimal crossing found, which a
  correct best score rules out), the port throws.
