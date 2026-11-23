# Lane BE: the BED family

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/be`. Model:
strong. Runtimes only, never the language (GC1 declared every param this lane reads).

**Task**: build [genome-formats](../../../design/genome-formats.md)'s `runtimes/bed` (its Decision: "What every runtime
shares" and "BED, bedGraph, BEDPE, chrom.sizes"; Bounds; Configuration rows marked BE), on BZ's io (`@molstar/io/common/bgzf`,
`common/tabix`) and `runtimes/internal/` (text, lines):

- BED (BED3 to BED12, `bed N+M` by `fields` with a column's 0-based number; BED12's blocks as pieces with the item's
  extent; identity `name@start-end`), bedGraph as `values`, BEDPE read whole and **indexed by both ends** (two
  `Int32Array` permutations; past `summaries.maxBytes` `too-large`), chrom.sizes as a coverage; whole, and by region at
  level 1 through tabix (one line parser either way); the fixed conventions (`base`, `ends` said once if written);
  issues `seq.source.bed.*` by the shared rule; every bound refused before it allocates;
- `runtimes/bed` an opt-in entry in no bundle that does not list it; its own size row under 30 kB in
  `packages/seq/_measure/sizes.mjs` (about 12 kB with io's BGZF and tabix);
- tests: the conformance suite (`data/_test/source-suite.ts`) on a generated fixture and on the recorded real slice;
  **a region read equals the whole read clipped** (a property over generated files and random regions); fuzzed and
  truncated files without a throw; BED12 pieces and `fields` by number; BEDPE pairs found from either end, deduped
  once; EC ids the design assigns to BE named in titles (read its "Edge cases" section: EC-190, EC-195, EC-196 at
  least); performance: 600,000 BEDPE rows indexed by both ends under 300 ms (heavy, `itWithin`; 12,000 quick);
- the **recording of real slices as pages** (design, "Real slices to record"): the recording script in `seqx-examples`
  (pages of 64 KiB keyed by file and page; a read not recorded fails naming its page) and a `sources.json` row with the
  licence. This lane builds the script and records **#1, UCSC hg38 `cytoBandIdeo.txt.gz`** (whole). **Downloads
  approved by the author, 2026-11-21 ("Fetch all four")**: fetch only #1 here, check its licence (UCSC's terms) when
  fetched and say it in the row; download into a fresh scratch folder, never run anything from it. GG and VC will reuse
  the script for #7, #8 and #2: keep its surface small and documented;
- stories: `Features/Genome formats/BED` (the file's head and the source node beside the canvas and the SVG; the
  cytobands as an ideogram's `bands`; a BEDPE's arcs if the pairs layer draws them; a misread file and its issues),
  light and dark, with the headless twin on the recorded and generated files;
- core's `_spec` Data section for the runtime; `files.md`'s and `sources.md`'s amendments listed as hot-file needs if
  they are under `.seq/`.

**Where**: a worktree on `seq/be`. Headers updated (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); hot files
never staged (`package.json`s, lockfile, `biome.json`, `.seq/`, `.v6/changelogs/`, `storybook/packages.ts`,
`packages/seq/TODO.md`): list them, with the exact edit (a new package export for `runtimes/bed` is a hot `package.json`
edit: give the line). Public API changes allowed in `packages/seq`; anything in `@molstar/core` or `@molstar/io`
stops the lane and comes back as a question. GG and VC may run beside you later in `runtimes/gff` and `runtimes/vcf`:
shared files (`runtimes/internal/`, `data/protocol.ts`) only grow, listed.

**Checks**: as `fx1-checks-and-carried.md` lists them, `pnpm -s -C packages/seq/core size` (the minimal and default
viewers unchanged by this lane), `pnpm test:heavy packages/seq/`, both perf tiers if a frame path changes, the four
`.seq/_checks` and `features-planned.ts --plan .seq/execution/fifth-slice/plan.md`. A fix's test fails on its parent.
Do not claim a failure is the base's without running it on the base.

**Report**: under 300 words: commits; checks with numbers; the story ids for the orchestrator's look; the public
surface added; the recording script's surface; the licence found; sizes; hot-file needs; the `.v6/changelogs/` text.
