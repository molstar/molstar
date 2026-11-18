# first-slice: the end-to-end slice of Seq\* (next) (2026-11-17 to 2026-11-18)

## Asked

- Build `packages/seq`: a reimplementation of Seq\* (SeqViewSpec and the viewer), from the oz principles' scaffold.
  Headers carry the Mol\* header and authors only. No mention of the previous implementation, no `.sexpr` file type, no
  backwards compatibility. Its Storybook group is "Seq\* (next)".
- Plan in the work folder, then run the loop.
- "continue all the way to the end of stage 7 as an overnight run."
- "also for gate decisions: if there is a way to make something configurable, then do that as choose your
  recommendation as the default."
- The stage 7 gate: "As recommended."
- "continue to finalize without waiting for me, including figuring out "foldback.md" for rules we should bring back to
  principles/design from what we learned. i will approve the foldback before committing it." Then: "foldback approved,
  commit it after lane J and fold it into the principles, etc."
- Plan the second slice, and draft the third, for direction.
- History is exposed by Seq\* as a port, because the viewer lives inside a bigger app: "yes, go ahead."

## Decided

- Gate 0 (2026-11-17, "approve both"):
  - meso\*'s tree core forked into `@molstar/core/repr`;
  - the slice's kinds and params;
  - three result types, with labels 1-based by default;
  - tiles of 1,024 bins;
  - FASTA read whole, so range reads wait for bigWig.

  These now live in tree-core, language-description, locations-and-queries and protocol.
- Configurable over fixed: a choice that could go more than one way becomes a config item with the recommendation as
  default (the author, 2026-11-17).
- The stage 7 targets stand ("as recommended", 2026-11-18). 3M bars gets a values direct path and `values.coarser`
  (default 2). Draw and hover allocate nothing warm, and a layout's allocation is bounded per frame. This lives in
  performance.md and the core `_spec` (Render).
- The fold-back (the 14 proposals of the first slice's fold-back, and evidence for the 15th) was approved, then applied
  to the manifest, the designs, the `_spec`s and the EC-id check.
- History is a Seq\* port, and the host owns the undo stack. `@molstar/ui/history` is one adapter. Lane N of the
  second slice applies this to changes.md.

## Delivered

| Lane | What | Merge | Verified by |
|---|---|---|---|
| A | the tree core in `@molstar/core/repr` (tree, schema, builder, codecs) | `fbc27cde5` | review (strong), 71 tests |
| B | the language, builder and issue catalog | `0eaf92fb2` | review, fix round |
| C | locations, element keys, queries, numberings | `5fa7c43d4` | review, fix round |
| D | the source protocol and the data store | `6e885495a` | review, fix round |
| C2 | half-open range sets in `@molstar/core/data/int` | `7bf4aeaf9` | review |
| E | the static model and the resolved view (paths, maps, follow) | `5592a743a` | review, fix round |
| F | render: draw lists, layer types, tracks layout, SVG | `e01f9f312` | review, fix round |
| G | the viewer's models | `02235fc67` | review, fix round |
| H | the viewer, Canvas2D, the p53 examples and stories | `ff7f8addb` | review, fix round; the look (light and dark) |
| fixes | buffers kept when they grow; the theme from the UI root | `855bddb8d`, `5f2168509` | tests failing before, screenshots |
| I | by region: async sources, summaries, FASTA, dynamic nodes, templates, nested sources | `efb8f9a94` | review, fix round |
| fixes | shared overlay sources; no stuck tiles; issues on change; `rested` after data | `306e70daa`, `27aaf5dfb`, `cdf667452` | the Gene by region test |
| M1 | measurements (bounds, allocation cases, bundle sizes) | `8af4eabdc` | the orchestrator's runs |
| J | the stage 7 fixes | `47a27cca7` | review (520 seeded frames identical), fix round, heavy and perf tiers |
| fix | `rested` once a rest, after the settled frame | `e16fe448c` | a test failing on its parent |

The final tree:
- 743 quick tests in under a second (`core/repr`, `core/data/int`, `packages/seq`);
- 32 heavy tests and 45 allocation cases;
- tsc, the test typechecks, lint, format and the scaffold checks (with `ec-ids`) all green;
- `pnpm build` green; `pnpm check:workspace` red only on `@molstar/ui`'s pre-existing `glass-wordmark` export, as on
  master.

## Measurements

Headless, Node 22, 1440 CSS px wide, dpr 2, after lane J:

| What | Result |
|---|---|
| Pan and zoom frames | every scenario's maximum under 2.3 ms (p53, 3M letters, cells and bars) |
| Hover | about 0.001 ms |
| SVG export (p53 Tracks) | about 1 ms |
| Validate and resolve 1,000 nodes | 7 ms |
| Allocation | 0 bytes warm for draw, hover, selection and cursor redraws; 25–40 KB a layout (p53), 0.7 MB (3M bars, bounded by the items in view) |
| Bundles, minified | 133 kB reading a spec, 229 kB headless SVG, 318 kB the viewer |

Browser frame times and first paint were not measured; they come in the second slice, in a visible pane.

## Learned

The approved fold-back, now in the manifest ("Folded back from the first slice"), the designs and the `_spec`s:

- Look at what is drawn before the merge, not only at the gate: five drawing bugs passed tests and review.
- A fix's test is shown failing on the fix's parent: one first test did not.
- A visual or timing observation states its conditions: a hidden pane runs no frames.
- Configuration tables in the `_spec`s, and the gaps queue lists only the unsure defaults.
- Merges name the refuted findings.
- A performance target names its path and its test, with a quick case when the path is built.
- Reused buffers are tested against fresh ones.
- What consumers share is held per holder and tested with two consumers live.
- Signals publish on change. A signal named for the picture fires on every cause; writing that rule found the
  double `rested`.
- A headless twin for each composed story, and a synthetic twin for each async source.
- An EC-id check.
- Verdict: hot files stayed with the orchestrator, and eleven lanes had no conflicts.

## Open decisions

For the author, with the recommendations as built:
1. Undeclared typed columns are tagged.
2. A cyclic range with start > end wraps.
3. The query decisions: `chain` as Nodes, quoted `.n` labels, the boundary rule, structural Location equality.
4. The protocol and store details: malformed answers not retried, `list`, row identity, per-consumer want sets,
   cancellation, `AssetReader`.
5. The `seq.lang.*` codes and alphabets.
6. The configurable defaults of every lane: render, viewer, canvas, the asset policy, summaries, nesting, templates.
7. The theme source (the UI root).
8. The selection tint drawn over the data.
9. Letters' density in CSS px.
10. `values.coarser`'s default: 2 as decided, but the review recommends 1, now that the direct path makes level 1
    cheap.
11. The asset manager change (credentials omitted, redirects refused, the signal, a streamed limit): a public API change
    to `@molstar/core`.
12. The TP53 exon coordinates in the test geometry, written from memory: check them against RefSeq.

These are re-presented at the second slice's Gate 0.

## Left open

- `packages/seq/TODO.md`:
  - projected hover is not drawn, and one window serves all leading sections;
  - a source-backed leader's length;
  - template re-expansion per arrival;
  - the ranges layout's allocation per item;
  - an allocation tier for the viewer;
  - what is deferred to the Genomes milestone.
- The Gene by region story was not rechecked in a visible browser after the last fixes: the pane was hidden, and its
  headless test covers them. The second slice's stage 0 does it.
- The drift notes are resolved into the designs or the TODO, and the drift folder is emptied.
- The second slice is planned and waits for its Gate 0; the third is drafted.

## Self-check

The claim was that the fold-back makes the principles truer to what was built.

- **Integrated**:
  - the rules on sharing, composition and signals, and the EC-id check, which stand whether or not the first slice's
    bugs recur;
  - the fix-fails-on-parent rule, which caught one miss in this slice;
  - the hot-file verdict, recorded with its count.
- **Carried**, as standing tensions in the manifest's "For the package's home":
  - configurable over fixed, against a config surface the author can confirm (about sixty defaults in one slice);
  - a look before every merge, against its cost per lane.
