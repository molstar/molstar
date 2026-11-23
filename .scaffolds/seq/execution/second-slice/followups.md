# Followups

- (M) `maxBytes` on `file://` reads in Node (`fetchData`'s file path ignores the four options).
- (M) The asset conformance suite sees a manager through `fetch`; a manager that downloads by other means is outside it.
- (L) The real p53 fixtures round-tripped JSON ↔ MessagePack in a seq package test (io cannot import them).
- (L) `repr.codec.msgpack.syntax` cataloged in core's `CodecIssues` beside the JSON and S-expression codes.
- (L) Duplicate map keys in MessagePack keep the last silently (JSON reports them): an option in `common/msgpack` later.
- (L) Payload alignment is by position: pad in `common/msgpack` if zero-copy decoding ever matters more.
- (L) `checkCodec`'s trees in core gain a leading-U+FEFF string, key and category entry; the codec's issue codes
  (`syntax`, `duplicate_key`, `key_type`) move into core's `CodecIssues`; the Codec's documentation says `encode` throws on
  a tree no reader would give back (a lone surrogate).
- (L) Encoder padding for dependable zero-copy; a canonical-category fast path; key-order determinism (with JSON).
- (K) Derive the next revision's ref index at the end of each apply: a structural edit on a 20,000-ref spec costs
  about 10 ms, mostly rebuilding the old root's index (0.7 ms with it cached); the 16 ms edit target.
- (K) `render/export.ts` still requests data for hidden layers.
- (K) `repr/tree/address.ts` (the child index) is reachable through core's `./*`; io imports `repr/tree/index` by
  path, so closing `./repr/tree/*` needs io moved to `./repr/tree` first.
- (R) Core's other `_test` folders have 73 type errors; `tsconfig.test.json` covers `data/immutable` only.
- The Complex dataset (examples.md, the author 2026-11-18): CDK2 with cyclin A, its interface, PAE and ligands with
  AlphaFill's transplants. Its fixtures and views go to the slice that builds pairs across chains, dot plots and ligands
  (the third slice's draft or later); the entries chosen then and checked against their sources.
- (K2 review) Many moves under one parent are quadratic in `applyPatch` (reversing 16k children: 105 ms alone, 229 ms
  in a merge); inserts and diffs are linear.
- (K2 review) `invalid` discards every change of theirs from `partial`, an unrelated one included (by design, N-8);
  lane N's UI should say so.
- (R+K) Law 3 fails with stored patches in a tree: the patch fuzz's derived mode, seed 3, case 572 (`fuzz(600, 3, true)`),
  the inverse refused as `repr.patch.dangling` with `inverseNeeds`. Skipped there with a comment; the orchestrator fixes it.
- (R+K) `RefSite.path` holds child indices, so an insert in a wide parent re-creates its later siblings' sites; nothing
  reads `path`: drop it or make it lazy.
- (O) An item's coordinate system as a prefix of location forms (`(item :id k)` gives elements): a `(system :item k)`
  form recommended. Refs in data nodes' params and in a section's `of`/`follow` are refused for now.
- (O) Declared numberings: the way back prints `ordinal` or `index` only until the language declares others.
- The frame lane's dropped-tick leak is fixed in `@molstar/ui/dom` (`frame.ts`, the scheduler remembered); lane O's
  twins can drop their per-test `DomHost.requestAnimationFrame` swap.
- `packages/sequences`' `stream/tiles/server.test.ts` fails under `--sequence.shuffle` (seed 9), on the parent too.
- Navigators (strip, miniature, ideogram, vertical overview; EC-110): no design yet (the author asked, 2026-11-18).
  A `navigators.md` after `stories.md`, its stage in the third slice.
- The stories library and the animation port (`stories.md`): placed at the next gate.
- `render/_test/frame-time.test.ts`'s per-frame max bounds flake under a parallel heavy run (lanes O and N: 18 and
  50 ms against 16); best of three, or a percentile instead of the max.
- (N) A track has no `hidden` param: in N's fix round (`track.hidden`, default false, the row collapsed; the
  orchestrator's call as for `layer.hidden`, in gaps).
- (N, for the stories lane) Its paused story steps against `stories.md`: a step that no longer applies shows the ops
  that still apply, flagged (today: the root without the step); `from`, a step's view and timing need language fields.
- (K3, for the stories lane) An op of a story step that no longer applies has no issue in `repr/patch`: the Seq side
  reports it on the step (stories.md). A host move can leave a step's later refs dangling in the step's own
  coordinates (its own move relocates the node): refused unless `allowDangling`; watch for it in the editing UI.
- (K3) A key typed as a ref with different types in different kinds is walked with each: a conformance case when a
  language has one.
- (Q review) `@molstar/core`'s `util/data-source.ts:482` opens a `file:` URL's path without percent-decoding: a root or
  file with a space or non-ASCII character is unreadable, and Windows paths open as `/C:/…` (fails safe today). When
  it moves to `fileURLToPath`, `checkFile` is reviewed again for Windows segments ending in a dot or space, device names
  and `:` streams.
- (Q) Bounds for regular expressions in search and for archives land with those features; nested-source breadth
  (fan-out W^depth) and joins allocating before `checkAnswer` (at most the file's size) recorded.
- (P) Drag snapping to whole codons over nucleotides (ux-foundations section 8); highlights holding their dynamic
  nodes; rects (dot plots) projected and drag-selected.
- (P) Projected selection bands are rounded at intron edges; the frame rule wants a piece's cut ends square (a
  projection-aware cut flag). The keyboard cursor reads as a full-height hairline below 1 px a cell (a minimum width
  or a caret).
- (S) A link followed before the viewer's spec is loaded is refused ("names nothing"), not kept; and a focus kept until
  the view is placed does not record its address (`focusAddress()` reads the shown range). The example story holds
  the link itself; the viewer should treat "no spec yet" as waiting and record the focus when the kept move lands.
- (S) The Features links story's `describeRun` prints a wait's message twice (the example uses `describeLinkRun`);
  `example-root.ts` keeps its own copy of the SVG-follow block (`_stories/follow-svg.ts`).
- (S) A `track.add` command in the viewer (the example inserts a track with a plain patch).
- (Stage 4 gate, the author 2026-11-18: third slice) A minus-strand gene following its protein moves opposite to it
  (maps.md: the hull in genomic order). A section's direction as a setting, `genomic`, `reversed` or `follow` (the
  leader's): a follower defaults to `follow`, a leader to `genomic`; a reversed axis counts down and says it is
  reversed; reversed DNA letters drawn as their complement (on by default). Built on reversed pieces of an axis.
- (Stage 4 gate, third slice) The protein and its gene on one axis, one under the other: the protein's letters drawn
  on the gene section through the codon map (a residue a codon wide, introns skipped), as a letters layer on another
  system's axis (today refused: `seq.render.layer.letters`). Not in the third slice's draft; its p53 views want it.
- (T) Export `projectPoint` from `@molstar/seqx-core/resolve`; drop the viewer's copy (`pathPoint` in `models/view.ts`).
- (T, third slice with the reversed follower) `Follow.point` runs a codon split by an intron straight from its first
  nucleotide, so some gene positions are never a follower's edge (a pan from the gene jumps them): follow each
  nucleotide of the codon.
- (T) The `genome-browser` preset without the add, subtract and extend drags (a drag pans there): decide with the
  presets; context-menu entries for `selection.add` and `selection.subtract`; a cursor moved off-screen on a follower
  revealing itself by moving the leader.
- (T) Alt+drag is bound (subtract): a host using Alt+drag for a copy-drag overrides the binding.
- (Stage 4 gate, the author 2026-11-18, third slice with collapsed introns) A gene following its protein drops its
  letters while an intron is in the window: the hull spans the intron at full width (a 20-residue window: 60 coding
  nucleotides at about 10.7 px each, 140 with the intron at about 4.6, below `density.letterMinPx` 7). Collapsed introns
  keep three nucleotides a residue; the default for a section following through a codon map.
- (U) The view history has no chrome yet (toolbar Back and Forward, their disabled state); `view.history` is reactive.
- (U) `Hit.hull` is new in core's public `Hit` (the item's pieces on the section): kept, recorded in the core `_spec`.
