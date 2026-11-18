# Plan: the third slice

**Rewritten 2026-11-18 as the second slice closed**, from what that slice taught (its record:
[records/second-slice](../../records/second-slice/2026-11-18-second-slice.md)). It is a plan for the author to agree at
Gate 0; no loop runs on it before then, and its loop file is written after.

Work id `third-slice`. Version 1's second milestone, **figures** ([version-1](../../design/version-1.md)): the axis
finished (direction, letters through a map, collapsed introns), the layouts beyond tracks, the remaining layers and
themes, the stories library and its animation port, navigators, PNG and the render CLI, and the other datasets with
their views ([examples](../../design/examples.md)).

## Where this is heading

| Slice | Milestone | In one line |
|---|---|---|
| first (closed) | the end-to-end slice | the tree core, the language, resolve, draw lists, Canvas2D and SVG, the viewer's models, data by region, measured |
| second (closed) | 1, core editing and addressing | patches, the history port, previews, addresses, refs and links, selection, security bounds, MessagePack, navigation |
| **third (this plan)** | 2, figures | the axis, the text, ring and block layouts, alignments, dot plots, arcs and logos, themes and legends, stories, navigators, PNG and the CLI, the datasets |
| fourth | 3, integration | structures (the port, structure maps, the Mol\* showcase), services, bigWig and range reads, the Nightingale reference, archives |
| later | 4, Genomes | layers across sections, synteny, contact maps, MAF, two assemblies on one axis, the genome formats |

## What the second slice changes in how this one runs

From the principles ("Folded back from the second slice"), written into this plan's gates and its loop file:
- **Each new design gets the self-check when it is agreed** (the axis, the layouts, the theme sublanguage, navigators,
  the stories library): its premise negated, what integrates planned, the rest kept as a tension in the design.
- **The interaction matrix** ([ux-foundations](../../design/ux-foundations.md), section 8): this slice fills the cells
  of a track or layer (select a row, focus a track), the selection's bookmark command, a bookmark's removal and label,
  rows and clades (collapse, pin), the dot plot's rectangle, the ring's window; Gate 0 names them, each gate checks
  them.
- **Precedents cover interactions** at Gate 0: what Jalview, ESPript, IGV, SnapGene-style plasmid viewers, R2DT/forna
  and the current Seq\* let a user do with these figures, each listed as built, planned or refused.
- **Use gates**: a short one after stage 1 (the axis and the first themes, on the p53 views), and a novice pass by a
  light agent before the author's use gate at stage 6. Each gate's findings are counted by kind (bug, missing design,
  look); the fold-back reads the count (the carried question: does an early use gate pay for figures?).
- **At every gate**: the size check (`pnpm -C packages/seq/core size`) beside the tests; the gaps raised since the last
  gate asked as yes or no; drawing checked against the data in the example views; heavy runs never two at once.

## Done when

- Every dataset has its views from examples.md running three ways (Storybook with canvas and SVG side by side, a
  quick-lane test that also checks drawing against the data, a spec file a reader without Seq\* could take):
  - **p53**: the protein and its gene on one axis (the protein's letters through the codon map), the gene following
    its protein reversed, introns collapsed; the text figure, secondary structure and variants on the tracks, the
    themes view;
  - **Globins**: the family alignment with logo and rows by tree, collapsed clades, the ESPript figure, hemoglobin with
    copies, the 10,000-row alignment with lazy rows and the keyboard cursor;
  - **Antibody**: Kabat against IMGT, the design against the wild type, CDR bookmarks with undo;
  - **Plasmid**: the ring with a linear detail following its window, a feature across the origin;
  - **tRNA**: arcs nested and crossing, the base-pair dot plot with crosshair hover and rectangle selection;
  - **Complex** (CDK2–cyclin A with an AlphaFill ligand): the two chains with their interface and the ligand's
    contacts, as examples.md describes.

  Each dataset has an editing view.
- The section layout interface carries four layouts (tracks, text, blocks, rings) and the dot plot as its own section,
  with no layout-specific code in the viewer.
- A section's direction (`genomic`, `reversed`, `follow`), letters through a one-to-one or codon path, collapsed
  introns: built on the axis's pieces, hit-testing, selection and export following.
- **Stories**: `stateOf`, `plan` and `sample` in `seqx-core`, the `AnimationPort` with its default and headless
  adapters, the p53 story played, scrubbed and exported ([stories](../../design/stories.md)).
- **Navigators**: designed (`navigators.md`, with its self-check) as one part over an axis or two, and built in all four
  forms of EC-110: the strip (a navigation track with window handles), the miniature, the ideogram, and the vertical
  overview of an alignment's rows (the author, 2026-11-18: "we forgot to even mention the "vertical maps"? ... it
  should prolly be built with the "navigation tracks""); each moves the view, idle when there is nowhere to go.
- Themes: field palettes (binned, continuous, log), rules, the residue schemes, legends, shared domains.
- PNG renders headless with the same picture as SVG; the render CLI turns a spec file into SVG or PNG under the headless
  asset policy (EC-126).
- Accessibility beyond the chrome: titles in exported SVG, the cursor and hover announced, forced colors.
- Measured: the 10,000 × 300 alignment, a ring's frame, the text layout's reflow, PNG export, a story's frame, the
  bundles at every gate. The author has used the views at both use gates.

## Gate 0 (the author, before any lane)

- **The datasets' exact data**: which globins, antibody, plasmid, tRNA and p53 structure; the Complex's entries.
  Recommendation: small public entries with permissive licenses, fixtures under 1 MB each, the large variants generated
  in code.
- **PNG in Node**: rasterize the SVG (resvg as a dev or optional dependency) or a Node canvas. Recommendation: rasterize
  the SVG (one picture, no native canvas in the library).
- **The render CLI's home**: `cli/` as `seqx-render`. Recommendation: yes.
- **Designs drafted first, each with its self-check**: the axis (direction, letters through a map, collapsed introns),
  the theme sublanguage (the `layer.color` slot's grammar, legends' place), `navigators.md`, the layout interface's
  extension for text, blocks and rings. Each agreed before its lane.
- **The interaction matrix's cells this slice fills**, and the precedent capability list.
- **Glyph atlas and band copy**: not needed by the second slice's browser numbers (frames about 1 ms). Recommendation:
  out of this slice unless the alignment's measurement asks.
- **Size targets with room**: headless rendering has 32 kB left under 300 kB, and this slice adds layouts, themes and
  layers to it. Recommendation: measure at stage 1's gate and bring a target change, if any, then.

## Testing strategy

- **Every lane**: its quick-lane tests in `_test/` beside the code; a story that composes more than one consumer has a
  headless twin on a fake clock; the example views' tests check drawing against the data (every item in view has its
  marks, every drawn mark is hit) at the fit view and across a tile edge; a fix lands with a test shown failing on its
  parent.
- **A path's target, in the stage that builds it**: a quick representative under `itWithin` (the heavy lane, the best
  of three): a text reflow, a ring's frame, a dot plot's frame, a story's `sample`, a PNG export; draw and hover warm at
  0 bytes (the opt-in tier).
- **New layouts and the axis**: property tests over random sections (every position placed once, hit-testing the
  inverse of placement, selection and export agreeing) and a record per layout (`recordDrawLists`) committed for the
  p53 and one other dataset.
- **Lanes**: quick in the edit loop; medium for the GL-free story twins and whole sessions; heavy for the bounds and
  the 10,000-row scaling. **At every gate**: quick, medium, heavy for the touched packages (never two heavy runs at
  once), the allocation tier, `pnpm -C packages/seq/core size`, the Storybook typecheck, lint, format, the `.seq`
  checks; `pnpm build` then `pnpm check:workspace` when `@molstar/core` or `@molstar/io` changed.

## Stages and lanes

Lane names are new to this slice. Each lane: a brief in `briefs/`, its `_spec` and Configuration table kept, a v6
changelog entry (written in its report), an independent review and one fix round, the look in a visible pane in light
and dark before the merge.

### Stage 0: carried from the second slice (orchestrator and one fast lane, beside Gate 0)

- [ ] C0 (fast) the carry-overs in `packages/seq/TODO.md`: export `projectPoint` from `resolve` and drop the viewer's
      copy; a link followed before the spec is loaded kept waiting (the address model treats "no spec yet" as waiting)
      and the focus recorded when a kept move lands; the Features links story's single message; `track.add` as a
      command; `view.back` and `view.forward` in the default toolbar with their disabled state; drawing checked against
      the data in the existing p53 example views.
- [ ] GATE: merged, reviewed, the checks green, sizes measured.

### Stage 1: the axis, then themes and the small layers (after Gate 0)

- [ ] AX (strong) the axis, as the axis design agreed at Gate 0 (with its self-check), first in the stage:
  - a section's direction `genomic`, `reversed`, `follow` on reversed pieces; a follower defaults to `follow`; a
    reversed axis counts down and says so; reversed DNA letters as their complement (on by default);
  - letters through a one-to-one or codon path (a residue a codon wide, a codon split by an intron in each piece),
    any other path still `seq.render.layer.letters` with its reason;
  - collapsed introns (the default for a section following through a codon map), and `Follow.point` following each
    nucleotide of a split codon;
  - hit-testing, selection, follow (both ways), focus and the SVG export through the new axis; the p53 views: the
    protein and its gene on one axis, the gene reading along its protein.
- [ ] TH (strong, after the theme sublanguage's design) themes: the `layer.color` slot's grammar (fields, palettes
      binned, continuous and log, rules), the residue schemes, legends as parts a host places, shared domains, size
      channels; the themes view.
- [ ] SL (fast, beside TH) the small layers: the ruler and numbering, secondary structure, variants, values' line and
      area, ranges' line styles and arrowheads; each with its fallback.
- [ ] GATE: merged, reviewed, the checks green, sizes measured (headless rendering's room, Gate 0's question), gaps
      asked, the matrix cells this stage fills checked.
- [ ] USE GATE (early): the author tries the p53 views for about ten minutes; findings fixed or placed, counted by kind.

### Stage 2: layouts (two strong lanes, after stage 1)

- [ ] TX (strong) the text layout and alignment figures, on the layout interface's extension agreed at Gate 0: groups
      of ten, lanes, inline labels, numbers per row, frames; reflow on resize within its target.
- [ ] RB (strong) blocks and rings: labels with leaders, a window a linear detail follows (EC-108, EC-114), a feature
      across the origin by its shortest arc.
- [ ] GATE: merged, reviewed, the checks green, sizes, gaps, matrix cells (the ring's window).

### Stage 3: pairs, grids and rows (strong, after stage 2)

- [ ] PG (strong) arcs and chords; the grid and dot plot as a section of its own, crosshair hover and rectangle
      selection.
- [ ] RW (strong) the rows layer: tree or data order, collapsed clades with a summary, pinning, lazy rows for 10,000
      (the store's want sets, bounded layout); the logo.
- [ ] GATE: merged, reviewed, the checks green, sizes, gaps, matrix cells (rows and clades, the rectangle).

### Stage 4: stories and navigators (SY after stage 2, beside stage 3; NV after RW)

- [ ] SY (strong) the stories library, as `stories.md`: `stateOf`, `plan`, `sample` and the channel kinds in
      `seqx-core` (property tests: `sample(plan, 0)` is `from`, `sample(plan, 1)` is `to`, monotone, deterministic,
      allocation-free warm); the frame slot over the shown state; the `AnimationPort` with its default and headless
      adapters; the view history's jumps and focus through it; the p53 story played, scrubbed and exported.
- [ ] NV (strong, after `navigators.md` and RW) navigators, one part by the axes it covers (positions, rows, or
      both), drawn from the coarse summaries the views already read, the window dragged and resized to move the view,
      idle when there is nowhere to go, a floored window drawn wider (EC-110), a part a host places:
  - the strip: a navigation track with window handles over a section's positions;
  - the vertical overview: an alignment's rows (columns compressed across), the Globins 10,000-row view's;
  - the miniature: both axes; the ideogram: a strip with a chromosome's bands.
- [ ] GATE: merged, reviewed, the checks green, sizes, gaps.

### Stage 5: datasets and views (fast lanes, one per dataset, after stages 2–4)

- [ ] DG Globins, DA Antibody, DP Plasmid, DT tRNA, DC Complex: fixtures (Gate 0's entries), builders, the views of
      examples.md, their stories and headless twins, drawing against the data, an editing view each.
- [ ] GATE: merged, reviewed, the checks green, sizes, gaps.

### Stage 6: output, accessibility and use

- [ ] OP (fast) PNG headless (the picture SVG draws, rasterized as Gate 0 decides) and the render CLI (`cli/`) under
      the headless asset policy (EC-126), with its cold-start target.
- [ ] AY (fast) titles in exported SVG, the cursor and hover announced, forced colors.
- [ ] NOVICE PASS: a light agent plays a user in a visible browser with a task script per view; friction reported.
- [ ] USE GATE: the author uses the views; findings fixed or placed, counted by kind.

### Stage 7: measure (orchestrator, light)

- [ ] M3 the heavy lane and the opt-in tiers; the browser in a visible pane at dpr 2 (the alignment of 10,000 × 300,
      a ring, the text reflow, a story's frame); PNG export; the bundles.
- [ ] GATE: the author confirms or moves the targets the slice touched.

### Finalize

- [ ] FINALIZE: reconcile (designs, `_spec`s, status, features, edge cases, `TODO.md`, drift emptied); the fold-back
      as its own step (this slice's lessons, the earlier fold-backs' rules against this slice's evidence, the count of
      what each gate found by kind, the self-check), decided by the author; the record; the fourth slice's plan.

## Order and parallelism

Gate 0 (with C0 beside it) → AX → TH and SL in parallel → early use gate → TX and RB in parallel → PG, RW and SY in
parallel (two heavy, one light at once), NV after `navigators.md` and RW (the vertical overview reads RW's
rows) → the datasets' lanes (two or three at once, fast)
→ OP and AY → novice pass → the author's use gate → measure → finalize. Never two heavy test runs at once.

## Hot files (the orchestrator's only)

- `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`, `biome.json`, `vitest.config.ts`;
- `storybook/packages.ts`, `.v6/changelogs/`;
- `packages/seq/*/package.json` and `tsconfig*.json`, `packages/core/package.json`, `packages/io/package.json`;
- a new `cli/` package's inventory entry (OP asks in its report);
- `.seq/`.

## Known risks

- **The layout interface** has carried one layout (tracks); the text layout and rings are its real test, and a change
  to it touches the draw lists' key and the hover index.
- **The axis lane first**: direction and collapsed introns change what an offset means on a section; hit-testing,
  selection, follow and export all read it. A twin per consumer, and drawing checked against the data.
- **Lazy rows for 10,000 × 300**: the store's want sets and the layout's bounded allocation meet here.
- **Headless rendering's size**: 32 kB of room for layouts, themes and layers (Gate 0).
- **Seq\* inside a bigger app**: legends, navigators and the animation loop stay host-neutral (parts a host places, a
  port the host drives).

## Out of scope

Structures and services (the fourth slice), bigWig, Genomes, writers beyond SVG and PNG, settings layers.
