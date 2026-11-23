# The genome ring as an extension: `rondo`, bundled and count-weighted

Lane RN's design (the fifth slice, stage 0; agreed by the author 2026-11-21 and promoted from its draft in `execution/fifth-slice/drafts/`). It amends
[extension-openings](extension-openings.md) (Part 1's "Enough for the genome ring?", the entry's
exports, extension bindings), [layouts](layouts.md) (Rings: the hole's title as a spec param) and
[themes](themes.md) (an extension's scheme); what it needs of the core is in GA's and RR's designs
(`section.include`, `navigator.switch`, the `system` domain before `include`, `pairsLevels`' arguments, `coarsenBins`,
a filter's rebuild). Built on Gate 0's answers as recommended (the fifth slice's plan): its own package `packages/seq/ext-rondo` (3);
bundling over **core's pairs level** (4); ring zoom as **a cut** (5: `section.include` in place, a second section as the builder's option); cytobands synthetic in tests,
UCSC's in the showcase if its licence allows, asked before the fetch (6); the navigator over several systems core's (8).
Built on the four designs before it and taking what they decided: [genome-coordinates](genome-coordinates.md) (GA:
assemblies, blocks, `placesOf`, runs, a ribbon's twist, the `block` and `system` fields' block-ordered domain, the
square's two-axis window over blocks × blocks, stubs on a ring RN's), [region-reads](region-reads.md) (RR: pairs tiled
in 1D by either end, the pairs level of bin-pair counts by strand, `PairsReader` over held and region pairs, items gated
by counts with hysteresis, "a lattice, not a clustering: clustering over bins is RN's picture"),
[genome-formats](genome-formats.md) (GF: BED's `fields` by column number for cytobands, BEDPE read whole and indexed
by both ends, GFF by region for a host's gene detail, the real-slice table), [across-sections](across-sections.md) (LX: ribbons and stubs are core's in a bridge only, on a ring RN's; the
twist drawn by corresponding positions; one hit rule, the drawn run as a `rect`; the `ribbons.*` look). On the extension
contract as the fourth slice's EX builds it ([extension-openings](extension-openings.md), Part 1: the
descriptor, namespaced types with plain params and a `fallback`, one `extension` entry) and CN (Part 2: navigator
pictures as registry entries).

## Intent

The author, 2026-11-19, on the genome ring: "Rondo-like, extension" (research-genome-ring, third slice, Decided): a
genome ring bundled and count-weighted, built as an extension package on the public surface, so the showcase is also
the test of the extension contract at genome scale ([version-1](version-1.md): "`ext:` layer kinds …
land before it"). The plan's Done when: "`rondo:ribbons` (wide by count, trumpets as cluster stems) and the bin-pair
level, cytobands, a hue per block, outward links to other chromosomes as stubs, the square navigator; 23 chromosomes and
12,000+ links at frame time". What a user does with it: see which chromosomes talk to which and how much, find the
strong clusters, open one chromosome (or a pair) to see its own links and where the rest go, follow a cluster to its
partners, read it with the keyboard, export it as a figure; and a reader without the extension still sees the links,
honestly simpler.

**Library and application** (the author, 2026-11-23, [the engine, not the game](README.md)): the extension is
building blocks (`rondo:ribbons`, `rondo:cells`, the `rondo:square` picture, `rondo:hues`, commands, a scoped binding)
on the library's public surface; how they are composed (a square beside a ring, what `opens` names, which menus a host
shows, which kinds of links a story layers, the named translocations) is application-level. `rondoRing` is one such
composition, a convenience, not the definition of a genome ring; the showcase is an application built from them.

## Precedents

| Precedent | Taken | Not taken |
|---|---|---|
| Rondo (CSIRO Data61, Garvan; research-genome-ring): chromosomes on a ring, connections clustered hierarchically at several resolutions, ribbons cluster to cluster thick by count, trumpets as cluster stems, a cluster's own links only zoomed in, one chromosome as its own ring, gene detail from about 1 Mb, filters by region, Q-value and minimum connections, the square navigator | stems as a coarser level of the bins (Rondo's resolutions), ribbons stem to stem by count, links inside a stem shown only when the level parts it, one chromosome (or a pair) as its own ring, the square over chromosome pairs, a minimum count by the data filter | gene detail by default (a host's track: "a library, not a viewer"); a Q-value filter at the bins level (RR's bins do not keep it: a Tension); BEDPE, BED, matrix export (a host's) |
| Circos: ideograms as arcs with gaps, cytobands, ribbons with twist, links coloured by chromosome, ribbons sorted at their ends, a large inner radius | the ideogram ring of bands, ribbons sorted by their partner's angle at each end, twist (GA's rule, LX's geometry), the hue a chromosome | static figures; rule blocks by regex |
| d3-chord and d3.ribbon (Bostock's chord diagram): groups as arcs, ribbons through the centre, hover on a group fades the rest | the ribbon through the hole as cubics, the fade on hover of a group (an emphasis wash) | the matrix as the only input (here pairs and bins) |
| Hierarchical edge bundling (Holten 2006; d3's `curveBundle` β) | `bundle` (β, 0 to 1): how far a ribbon is pulled toward the centre, so ribbons of one stem stay together longer | a hierarchy given by the data (here the stems are the hierarchy: block, stem, bin) |
| UCSC's cytoBandIdeo and its ideogram | the `bands` layer (built: `giemsa`, `acen` triangles), the file as BED with `fields: { stain: '4' }` (GF) | — |
| The current Seq\*'s ideogram (a recording of UCSC's hg38 cytobands, its provenance kept, its licence not) | the data's shape | its bytes without a licence row (Data, below) |
| This codebase: the ring of blocks (`rings-fixtures`), `scope: 'section'` chords through the hole placed by `ctx.place` and the polar warp (`Placement.warpOf`), `pairs.maxArcs` and its cull, the `bands` layer, GA's square, RR's `PairsReader` | every one, read through the `extension` entry | a private reduction, layout or navigator geometry in the extension |

### Interactions in the precedents

| What a user does there | Here |
|---|---|
| hover a chromosome: its links lit, the rest faded (Rondo, d3-chord) | built: hover a block's arc or its stems; an emphasis wash over the rest of the hole |
| hover a ribbon: its two clusters and count (Rondo, Circos tooltips) | built: the card (Decision, Hover) |
| click a chromosome to open it as its own ring (Rondo) | built: a block label's menu "Show only chr4" (core's `section.include`, GA's), the square's click, or a double-click on a stem with `opens` (Decided, 7) |
| open a pair of chromosomes (Rondo's square) | built: the square's cell click (`navigator.switch: include`) |
| back to the whole genome (Rondo) | built: undo (one patch), the square's title click, the label menu's "Show all" |
| see a cluster's own links by zooming in (Rondo) | built: the chromosome's ring reads finer bins, then items, and the stems part |
| filter by a minimum number of connections (Rondo) | built: the layer's data filter on `count` (EC-225) |
| filter by Q-value (Rondo) | partly: on items, and on bins only where the pyramid is rebuilt for the filter (Decided, 8) |
| filter by region (Rondo) | built: `section.include`, or a selection's ribbons in the accent |
| gene detail from about 1 Mb (Rondo) | refused as a default: a host adds a genes track (a GFF by region, GF's `gff` runtime) to the chromosome's ring |
| export SVG, PNG (Rondo, Circos) | built: the export draws the ribbons and stems as paths |
| export PDF, BEDPE, BED, a matrix (Rondo) | refused here: a host's |
| rotate the ring (Circos) | built: the section's `rotation` (layouts) |
| reorder chromosomes (Circos's order) | refused as a gesture: the assembly's `order` by patch (GA) |

## Decision

### Whose each choice is

By [axis](axis.md)'s rule. What the figure says (which data, which width, which section a double-click
opens) is a spec param, plain on the extension's own type (extension-openings, Decided 2); how stems and ribbons are
built and bounded is the extension's config (`rondo.*`, the descriptor's `config`), the look shared with bridges is
core's `ribbons.*` (LX); what a click on the square does is the navigator's (core). No port: the extension reads the
store, the placement and the theme through the `extension` entry.

### The package

- **`packages/seq/ext-rondo`, `@molstar/seqx-ext-rondo`** (Gate 0, 3), in the inventory and the root references (hot
  files, the orchestrator's). Two entries: `.` (`rondo: SeqExtension`, everything headless: layer types, the picture,
  the scheme, params, config; the builder helper `rondoRing`) and `./viewer` (`rondoViewer: SeqViewerExtension`, the
  above plus commands, menu entries and the binding). Why two: the CLI and `createSeqContext` compose the first and pay
  nothing for the viewer's parts.
- **Imports only `@molstar/seqx-core/extension` and `@molstar/seqx-viewer`'s extension types**, checked by its
  `noRestrictedImports` (extension-openings, Decided 5). Why: the package exists to prove the contract; an import past
  it is a contract gap to bring to the core, never a reach into `src`.
- **`requires`**: the language version GC1 lands (0.6 expected); `version` 0.1.0. A host lists it
  (`spec.extensions: [rondoViewer]`) or registers it late (`viewer.extend`); never in the minimal or default viewer.
- **Its own size row**: under 30 kB minified (GR's estimate: 650 to 950 lines), the minimal and default viewers not
  grown. No lazy chunk inside: the package is itself the opt-in.
- **What it costs the core and the viewer** (the default viewer is a few hundred bytes under its bound): `section.include`
  in the language and resolve (about 0.5 kB, in reading a spec: every reader needs it), `section.center` (about
  0.1 kB), the hooks below (about 0.3 kB: three optional members and their callers); `navigator.switch`, the
  `section.include` command and the label menu's two entries (about 1 kB) **lazy**, in the chunk that loads with a
  two-axis navigator over blocks (GC3's) and the label menu's (HP2's menu chunk). Each a row measured at GC2's and
  GC3's gates; past a bound, to the author with the feature (plan, Sizes).

### What it ships

| Entry | Kind | Fallback | Reads |
|---|---|---|---|
| `rondo:ribbons` | layer type, `scope: 'section'`, on `rings` (and a `tracks` row, below) | `pairs` (written by `rondoRing`) | `pairs` (RR's `PairsReader`: bins at a level, items in a region) |
| `rondo:cells` | layer type for the square's picture: block × block cells | `none` (written for the rule; unreachable, since a reader without `rondo` draws the whole picture as `axis`) | `pairs` bins at the coarsest level, by block (GA's square) |
| `rondo:square` | navigator picture, `axes: both` | `axis` (the lattice and labels, fixed) | its tracks: one `rondo:cells` layer over the bound layer's data |
| `rondo:hues` | scheme over categories: n hues evenly round the oklch wheel, n the domain's size | the layer type's default colour | — |
| `rondo:partner` | command (`navigate`) | — | the stem under the cursor |
| `rondo:select-partners` | command (`select`) | — | the selection's ribbons |
| the stem menu, `opens`' double-click | viewer menu entries, a binding | — | — |

No node kind, no data shape, no format runtime (extension-openings, Part 1): RR's bin-pair counts are `pairs`, so
`rondo:binpairs` is not built (Gate 0, 4).

### `rondo:ribbons`: what it reads

- **The extension surface is EX's design, not yet built**: no `./extension` export or `SeqExtension` exists today
  (the fourth slice's EX is pending); this design assumes extension-openings Part 1 as agreed and lists below, as core
  changes, the hooks it needs beyond it (Consequences, `seqx-core`).
- **Two reads, RR's rule with this layer's numbers**: bins at the **stem level** (below), always; items at a tile size
  when the stem pairs with an end in view sum to at most `rondo.maxRibbons` (1,024), kept to `pairs.itemsHysteresis`
  × it (1,280), bins drawing until they land (EC-169). So RR's level function takes the gate and the bin width as
  arguments (`pairsLevels(view, { binPx, gate })`, exported; RR's numbers its defaults). Why its own gate: every item
  is one count, so past the drawn bound the cull could only rank items by span and would fan out links for no visible
  reason; at the gate every item read is drawn. The fallback `pairs` keeps RR's 4,096: two readers, each gating at
  what it draws. Held and region pairs alike (`PairsReader`).
- **On blocks**: ends placed by GA's `placesOf` (every block that shows the position, each place counted against the
  bound), runs and unread blocks by the joint rule (an end on an unread block at its place, read from its other end's
  tile). Two assemblies on one ring (EC-210) read a sided source by each side's assembly (RR's sides).

### Stems: a coarser level of the pyramid (the trumpets)

- **A stem is one bin of the stem level**, clipped to its block: the stem level is `2^j` × a level the source lists
  (the finest listed level `L` at most the target, then the largest `2^j · L` at most the target), the target being
  `rondo.stemPx` (16 CSS px) of positions at the ring's outer circumference in **CSS px of the section's layout**, never
  device px. Stem bins are aligned on position 0 of each system, so a stem is a union of whole `L` bins and its counts
  are exact sums (RR's levels are nested). Why a level and not a clustering by gaps: dense data leaves no gaps, so
  gap clustering cuts fixed chunks arbitrarily; a level is exact, cached by level (one `coarsenBins` reduction of the
  `L` bins by `2^j`, RR: bins to bins, kept by the entry's version), the same on screen, in the SVG and at any DPR, and changes on a resize only
  when the target crosses a power of two. RR's Tension ("a lattice, not a clustering") stays: a cluster across a stem
  edge is two stems, accepted for exactness (Tensions).
- **Its counts are sums**: a stem's count is the pairs with an end in it (each pair once: RR's dedupe across tiles and
  blocks); a stem pair's count the sum over the `L` bin pairs between them, by strand (straight and reversed never
  merge). Empty stems are not drawn.
- **Recomputed** when the stem level changes (a resize past a power of two, `section.include`) and when a tile lands,
  never on hover or on its window's move (layouts, "The ring under a moving window").
- **The trumpet**: from its mouth (the stem's arc on the band's inner edge, clipped to its block) to its neck at
  `rondo.stemDepth` (0.12 of the hole's radius) inward, two cubics flaring from the mouth's width to the neck's (below),
  filled in its block's hue at `ribbons.fillOpacity`. A dense narrow stem flares inward, the same geometry.
- **Links inside one stem** (both ends in it) are not ribbons at that level: counted in its card ("+12 within"), drawn
  once a finer stem level parts them (Rondo: "internal connections shown only zoomed in").
- **`stems`** (param): `level` (the default, above), `none` (no trumpets: ribbons from the band's edge at the stem
  level's bins, Circos's look).

### Ribbons: stem to stem, wide by count

- **A ribbon a stem pair a strand**, from neck to neck through the hole: each side a cubic whose control points sit at
  `(1 − bundle)` of the neck's radius along each neck's angle (`bundle`, β, 0.7: Holten's strength; 1 passes through
  the centre as d3.ribbon, 0 a straight chord).
- **Width: one px a count in the view, no floor inside a neck.** `k = rondo.maxWidthPx (24 CSS px) / the largest stem
  count in view`; a neck is `k ×` its stem's count, at least `ribbons.minPx` (1 device px) as a whole; inside it, each
  ribbon whose `k × count` is at least `ribbons.minPx` is drawn at exactly that width, and the stem's ribbons under it
  are one **fan** of their summed width (exact, however thin), hit as the stem ("+212 weaker links to 31 stems"). So
  a neck is always the sum of its ribbons and its fan, and a hub stem with hundreds of partners keeps its width and its
  count. Why this gives and not the cap or additivity: the cap keeps the hole readable, additivity is the check a
  reader makes; the floor is only for being seen, which the fan serves. Why linear: RR's arcs use square root for
  strokes that never add. (`width`: `count`, or a field summed, RR's `sum`.)
- **The twist at a neck** (LX's rule restated for a ring): a ribbon's slice of a neck has a **low edge** and a **high
  edge**, toward lower and higher positions in its block's display direction (a slice is a share of the count, not
  positions, so the edges stand for the slice's ends of the stem). On a ring two arcs face each other with opposite
  tangents, so an uncrossed ribbon joins A's low edge to B's high edge and A's high to B's low (d3.ribbon's and
  Circos's flat ribbon). With `s = strand × display direction of A × display direction of B` (GA): `s = +1` uncrossed,
  `s = −1` crossed (low to low, high to high), strand 0 uncrossed. Items, whose ends are positions, use their extents'
  ends the same way: on a ring an uncrossed item joins A's start to B's end, so LX's "corresponding positions" (start
  to start) is **not** what a ring draws; what holds is GA's meaning (an inversion crossed, a reversed block untwisting
  it), for items and bins alike, and a property test holds the rule equal to GA's formula over every direction
  combination (corrected 2026-11-23: RX's review found the earlier "holds for items exactly" false; the code follows
  GA).
- **Order in a neck**: its ribbons side by side in the order of their partner necks' angles clockwise from the stem
  (Circos's sorting), then its stubs by the partner's `order`, then its fan last. Fixed: it keeps ribbons of one stem
  from crossing near it, and the fan where the weakest would sit.
- **Items** (past the gate, a chromosome's ring): each pair a ribbon of one count from its ends' own extents, still
  grouped into its stems' necks, so the bundle holds while its members show.
- **Bounded**: at most `rondo.maxRibbons` (1,024) ribbons, the strongest by count across the view (ties by span, then
  key), the rest in their stems' fans; said once (`rondo.ribbons.limit`, info: "612 links drawn in fans"), not core's
  `seq.render.pairs.limit`, whose words are arcs and a density strip. Below the items gate nothing is culled.
- **Drawn widest first**, so narrow ribbons lie on top and stay hittable (LX's rule); fill at `ribbons.fillOpacity`
  (0.35), no outline.

### Stubs: links to blocks the ring does not show

- **An end with no place** (its system not in the section: another chromosome on a chromosome's ring; a contig
  `include` left out) makes its pairs **stubs**, aggregated a stem a partner sequence: one tapered ribbon from the neck
  toward the centre for `rondo.stubShare` (0.2 of the hole's radius; LX's `ribbons.stubShare` is a share of a bridge's
  band, so not reused), in the partner's hue, its width by count in the neck as a ribbon's (a fan member when under the
  floor), labelled at its tip with the partner's name and count when the tip has `rondo.stubLabelPx` (24 CSS px) of
  room ("chr7 · 37"). Why inward: the labels own the outside of the ring. Why aggregated: a chromosome's ring may have
  thousands of outward pairs; one stub a partner a stem is the question asked.
- **Hit as the drawn run** (LX): `(rect :a <the stem's range> :b <the hull of the partner ends>)`; hovered as "→ chr7 ·
  37 links, not shown". **`stubs`** (param): `partner` (the default), `none` (core's caret instead).

### Colour: a hue per block

- **The default channel is `(field :name system)`** of the data's `a` end (LX's rule for one section, taken), stems by
  their own block's. GA gives `system` a domain of the assembly's names in the assembly's order, taken
  **before `include`** (the assembly's and the section's), so chr7 has one hue on the genome's ring, on chr4's ring and
  on a stub pointing at it. Why `system`, not `block`: a block's domain is the section's blocks, and
  a chromosome's ring has one. This amends EC-191's block-ordered domain (GA) for `system`: the assembly's order,
  never the section's blocks.
- **The palette**: the layer's default is the theme's default categorical palette (`categories`, 8 tokens, safe in both
  themes, `repeat` true). On the whole ring neighbours are always distinct (chr1, chr9 and chr17 share a hue, eight
  apart); **on a subset they may not be**: a pair eight apart (chr1 and chr9, chr4 and chr12) draws in one hue. Stated,
  not fixed: re-deriving the domain from the blocks shown would give chr7 two hues across two sections, which is worse;
  the pair's arcs are labelled, its stubs named at their tips, and `rondo:hues` has no repeat. `rondo:hues` (n hues evenly round
  the oklch wheel at the categorical tokens' lightness and chroma, its dark variant derived) is offered for 24
  chromosomes, **not colour-blind safe** (adjacent hues 15° apart), so the picker says so and its substitute is
  `categories` (themes, Schemes); the showcase spec names it, the layer never defaults to it (Decided, 5).

### Cytobands and the hue ring

- **Cytobands are core's `bands` layer** (giemsa, `acen` triangles) on the outer ring, over a BED annotation with
  `fields: { stain: '4' }` (GF) whose `target` is the assembly (GA). Nothing in the extension.
- **The hue ring**: a thin second `bands` layer over the same data coloured by `(field :name system)` (the bands cover
  each chromosome, so the ring is solid hue with the centromere's notch); without bands, `rondoRing` writes a ranges
  annotation of whole sequences from the assembly's inline `sequences`. Why no layer type of its own: the bands already
  cover the block, and a reader without the extension draws it.
- **Synthetic bands for tests** (Gate 0, 6): `syntheticBands(sequences, seed)` in `seqx-examples` (alternating `gneg`
  and `gpos25`…`gpos100` of lengths drawn from the seed, an `acen` pair at 40 % of each length, a `gvar` on the last
  third of a chromosome named `Y` when present), so every test and story without the recording is deterministic.

### The ring zoom: a cut, by `section.include` (Gate 0, 5)

- **Opening a chromosome changes which blocks a section shows**, not a viewport: a ring shows its whole axis and its
  key never holds a viewport (layouts, GA). So the ring reads GA's **`section.include`** on a section of an assembly (the
  assembly's `include` grammar: names, patterns, `{ min_length }`; intersected with the assembly's own) and its core
  command (change kind `patch`, one undo step): "show only chr4" patches `include: ['chr4']`, a pair
  `['chr4', 'chr7']`, "show all" removes it. Why a section param and not new `block` nodes: one param to patch and to
  undo, the same grammar as the assembly's, the hue domain and `order` untouched; and it gives GA's standing
  "switching on a section of one chromosome" its answer for rows too (UCSC's one chromosome at a time).
- **In place by default** (Rondo's): one ring section; `rondoRing({ detail: 'below' })` writes a second ring section
  of the same assembly under it, its `include` the one patched, so the genome stays in view. Both are the writer's.
- **A cut**: the ring is laid out again with its labels and stems (GA: a change of the block list is a cut; stories);
  no tween (Gate 0, 5: a channel only if the use gate asks). A story step "zoom into chr4" is a step whose patch sets
  `include`; "back" the reverse step.
- **At the chromosome's scale** the same reads pick finer bins (chr4's 190 Mb on the circumference is about 42 kb a
  device px), its stems about 1 Mb, then items once its stem pairs in view fit `rondo.maxRibbons`; links to other chromosomes become
  stubs. History: viewports by `{ system, position }` (GA) survive the change; the selection's addresses unchanged.

### The square navigator

- **Core's square** (GA, Gate 0, 8): a navigator with `axes: both` over an assembly's blocks × blocks, its lattice laid
  out by the run rule. **Its lattice is the assembly's** (after the assembly's `include`, never the bound section's),
  so after "Show only chr4" the square stays 23 × 23 with every cell but chr4 × chr4 washed, and **washed cells stay
  clickable under `switch: include`** (the wash says "not shown", not "unreachable"). GA's square gives the two things the ring needs: a navigator
  may bind **one section on both axes** (`of: [ring, ring]`: its reachable cells are the section's block pairs, the rest
  washed), and **`navigator.switch`** (`none`, the default; `include`): a click on a block or a cell runs
  `section.include` on the bound section with that block or pair (a click on the cell of what it already shows: no
  change), a click on the navigator's title (the assembly's name) or Home on its focus removes the `include`. Why in
  core: a ring's window is idle unless followed, so a centring click does nothing there, and switching a section of
  one chromosome is wanted on rows too; the extension only draws.
- **`rondo:square` draws the cells**: a picture (CN's contract) returning one track with a `rondo:cells` layer over the
  bound `rondo:ribbons` layer's data (named by the picture's `of_layer` param, else the first `pairs`-shaped layer of
  the bound section): a cell a block pair coloured by its count on a square-root ramp (`seq-ramp-*`), the diagonal
  (links within a chromosome) included, upper and lower mirrored for unsided data; sided data (two assemblies) fills
  both triangles from its own sides. Read as GA's square reads pairs: no 2D summary, so one want a block at the coarsest pairs level (a tile holds every pair with an end in it), summed into one
  image of block × block cells, reduced to its pixels by the shared reduction; its cost is the blocks, not their
  square, and the image is always drawn. `navigator.square.maxCells` (4,096) is GA's: past it the lattice's lines are
  not drawn, never a cap on what is read. Hover: "chr4 × chr7 · 37 links"; the window outline over the cells
  the section shows.
- **Without CN** (Part 2 cut in the fourth slice): the square draws core's `axis` picture (the lattice and labels) with
  the same switch; `rondo:square` is dropped from the package, not rebuilt privately.

### Hover, selection, emphasis, the keyboard

- **Hooks the extension needs in core** (beyond EX's contract; each an optional `LayerType` member, so built-ins need
  nothing): `describe?(hit, lay): readonly string[]`, the hover card's lines for a hit of the layer (the card's
  `viewer/src/describe.ts` is private; it calls the hook first, its own text after); `summary?(lay): string`, the
  layer group's ARIA description, set when the frame settles; `announce?(lay, system, position): string | undefined`,
  what the keyboard cursor adds at a position (today a layer offers only `items()`). Why hooks and not exports: the
  card, ARIA and the cursor stay the viewer's; a layer only supplies its words.
- **Hover a ribbon**: its two stems' ends lit on the ring (a pair closes to its ends), the ribbon at full opacity with
  `seq-hover`'s outline; card "chr1:12–28 Mb ↔ chr17:40–52 Mb · 140 links (contact) · 1 Mb bins · zoom in to see
  them" (RR's bin text; items: the pair's own). **Hover a stem**: its ribbons lit and the rest of the hole washed
  (`seq-wash` at `rondo.fade`, 0.6: d3-chord's fade), card "chr4:50–58 Mb · 312 links to 9 chromosomes, 12 within".
  **Hover a block's arc or label**: every ribbon touching it, the same wash. All overlay passes: the base is not drawn
  again.
- **Hits** (closed, extension-openings): a ribbon is its drawn run, `(rect :a … :b …)` (LX's one rule); a stem its
  range on its system; a stub as above. Polygons flattened at layout (each cubic `rondo.flatten`, 12 segments), a
  uniform grid over the hole (`rondo.hitCellPx`, 16 CSS px cells) listing candidates, tested topmost first: under 2 ms, 0 bytes warm.
- **Select**: click selects the hit's element (a ribbon both ranges; a stem its range); the selection's ribbons (with an
  end inside it) drawn in the accent; `rondo:select-partners` adds their far ends (one selection change). Why a
  command: partners are added when asked (LX's `selectAcross`).
- **Emphasis by Location**: a story step's emphasis "pairs touching chr4" is a Location (chr4 whole); the layer lights
  every ribbon with an end in it and washes the rest. No new query form (research-genome-ring's ask answered).
- **The keyboard**: the ring's cursor (built: Left and Right round the ring) announces the stem under it ("chr4,
  50,120,000 · stem of 312 links"); `rondo:partner` (`.` and `,`, checked against the bindings at RX) moves the cursor
  to the strongest partner neck's end and on to the next by count, announced ("chr17, 40,000,000 · 140 links,
  reversed"); Enter selects the stem; Shift+F10 opens the stem's menu where the viewer has element menus
  (ux-foundations: "context menus stay later"), its entries the partners by count and "Show only this chromosome".
  ARIA: the layer's group described by its counts in view ("612 ribbons, 23 chromosomes; strongest chr1 ↔ chr17,
  140"), set when the frame settles; ribbons not focusable one by one.

### On other layouts

- **A `tracks` row of blocks**: the same stems below the band and ribbons as filled arcs between them (the polar warp's
  absence is the only difference: `ctx.place` and the cubics), LX's "ribbons inside one row: an extension may draw
  them". **Text or a grid strip**: not drawn, `rondo.ribbons.layout` (warning, once); its fallback is for readers
  without the type, not for layouts. Why rows: a genome in a row with its links above is Rondo's ring unrolled, and
  costs nothing more.

### What a reader without the extension draws (EC-235)

| Part | Drawn as |
|---|---|
| `rondo:ribbons` (`fallback: pairs`) | core's chords: a bin pair one chord stroked by count (RR's picture 3), items as chords when they fit, culled at `pairs.maxArcs`; an end with no place a caret |
| `rondo:square` | core's `axis` picture (lattice and labels), the window and `navigator.switch` as with it |
| `rondo:hues` | the layer type's default colour (the categorical palette over the same field), said |
| cytobands, the hue ring, `section.include`, the zoom | core's, unchanged |
| the stem menu, `rondo:partner`, `opens` | absent; the label menu's "Show only" (core) remains |

Said once a view and namespace, one issue for every part (EX's rule): `seq.view.layer.extension` ("`rondo` 0.1 missing:
1 layer drawn as `pairs`, 1 navigator picture as `axis`, 1 scheme as the default"); the spec's `extensions:
{ rondo: '0.1' }` written by `rondoRing`, which also writes the colour channel explicitly (`(field :name system)`
with its palette), so the fallback's colours are the same by the spec, not by two defaults agreeing. Why this is honest: the fallback draws the same bins in the same colours;
only the bundling and the stems are missing.

### The builder helper

`rondoRing(b, { assembly, links, bands?, detail?: 'inplace' | 'below', square?: true, hues?: true, opens?: true })`
writes the ring section(s), the bands and hue tracks, the `rondo:ribbons` layer with `fallback: 'pairs'`, the square
navigator with `content: 'rondo:square'`, `switch: 'include'`, `of: [ring, ring]`, the section's `center: none`
(below) and the spec's `extensions`. Spelled, never inferred; nothing added the call does not name (no genes track).

### The hole

The ribbons fill the hole, where the ring draws its title (`rings.center: 'title' | 'none'`, a `RenderConfig` preset,
built). Layouts' preset is **promoted to a section spec param `center`** (the same values; absent: the preset), written `none` by `rondoRing`, the
assembly's name and length then in the section's label line. Why a spec param: the title over ribbons is the figure's
choice, and a reader without the extension draws chords there too.

### Data for the showcase (no fetch in this lane)

- **The genome**: GRCh38's 23 chromosomes (1 to 22 and X; Y joins as the 24th by `include`) as an `assembly` with
  inline `sequences` (names and lengths: facts from the GRC's assembly report, typed into the generator, no file).
- **12,000+ links, generated** (GD's; RR: "12,000 generated links, held and as a source: one picture"): a seeded
  generator over the lengths, 60 % within a chromosome with a distance decay, 40 % between, concentrated on a dozen hub
  pairs and hot spots so stems and strong ribbons exist (`type`: `contact`, `other`; a `q` field), written held and as a
  BEDPE source (one picture, RR's property).
- **Cytobands** (Gate 0, 6; GF's table): **the fetch to ask**: `https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBandIdeo.txt.gz`
  (a few kB), recorded whole with its `sources.json` row in `seqx-examples`' fixtures and a recording script in
  `packages/seq`; **the licence to check**: UCSC's conditions of use for its database tables (the band coordinates come
  from NCBI's ideogram data). The repository already holds these bytes as the current Seq\*'s recording (provenance URL
  kept, no licence row), which cannot be cited from `packages/seq`; the ask is the licence first, then one fetch so the
  recording and its script live here. **If UCSC's terms do not allow it**: NCBI's ideogram file for GRCh38 (a US
  Government work, public domain; its path to be confirmed at the ask), converted once to BED.

## Edge cases this design takes or assigns

- **Taken here**: EC-108 (the bins' columns on the outer circumference; labels outside, stubs inside), EC-169 (bins
  drawn until items land; a pending block's stems as the band's hatch), EC-171 (no "load anyway": bins and stems
  instead, RR's outcome), EC-206 (the twist on a ring by GA's direction rule), EC-210 (two assemblies' synteny
  ribbons inside one ring, sided), EC-218 (a hue per chromosome, stable across sections), EC-225 (a minimum count by
  the data filter on bins), EC-235 (the fallback table above).
- **Taken from the other designs, not redesigned**: EC-62 and EC-208 (RR: pairs read by region from either end), EC-102 and
  EC-184 (GA: gaps at scale, runs; switching chromosomes, here by `section.include`), EC-110 (navigators: the idle
  window of a ring nobody follows), EC-12 (the shortest arc on a focus).
- **Assigned**: a Q-value filter at the bins level (RR's rebuild of held and summarized pyramids, Decided, 8); gene
  detail (a host's genes track over GF's `gff` runtime by region).

## Defaults, with a picture

What a reader sees from `rondoRing(b, { assembly, links, bands })` and nothing else.

1. **The genome at the fit, with the extension**: 23 arcs clockwise from the top, cytobands outside a thin hue ring,
   stems as trumpets inside, ribbons stem to stem wide by count, the strongest on the bottom of the stack.
   ```
                  chr1
            chrX ╭▓▓▓▓▓▓╮ chr2
              ╭─╯╲▼▼  ╲╰─╮          ▼ a stem's mouth narrowing to its neck
        chr22 │   ╲━━━━╲  │ chr3    ━ a ribbon: chr1 ↔ chr17, 140 links
              │    ╲  ┃ ╲ │         ┃ a thin one: 3 links
              ╰─╮   ▲▲ ▲╭─╯ chr4
                 ╰▓▓▓▓▓╯  …
   ```
2. **The same spec without the extension**: the arcs and bands as above, bin-pair chords stroked by count (RR's
   picture 3), no stems; the issues panel says once "drawn as `pairs`".
3. **One chromosome** (Show only chr4): chr4's bands round the whole ring, its own links as ribbons (items when they
   fit), its outward links as stubs inward in their partners' hues, labelled at the tips.
   ```
             chr4 ▓▓▓▓▓▓▓▓
           ╭──▼───────▼──╮
           │  ╲   chr7·37 │     stubs: short tapers in chr7's and chr11's hues
           │   ╰━━━╮      │
           │ chr11·12 ▲   │
           ╰─────────────╯
   ```
4. **A pair** (the square's chr4 × chr7 cell): two arcs, their links between them as ribbons, the rest stubs.
5. **The square**, beside the ring (the host's place): 23 × 23 cells by count on the ramp, the diagonal darkest, the
   cells the ring shows outlined (the window), the rest washed only when the ring shows a subset.
   ```
        1 2 3 4 …
     1 ▓ ░   ▒
     2 ░ ▓ ░
     3   ░ ▓
     4 ▒     ▓   ┌┐ the ring's cells
   ```
6. **Hovering a stem**: its ribbons lit, the rest of the hole washed, the card with its counts.
7. **Two assemblies on one ring** (EC-210, peach and grape): GA's group gap, synteny ribbons inside the ring, an
   inversion crossed.
8. **Crowded** (3,000 contigs): GA's runs; stems only on blocks at least a device px wide; ribbons to unread blocks at
   their places in the run.
9. **Pending, failed, nothing**: the band's pending hatch where bins have not landed, no stems there; the placeholder
   with its issue; no pairs: an empty hole, its ARIA description "no links".
10. **Dark**: tokens only (the categorical palette's dark variant, `rondo:hues`' derived dark, `giemsa` from
    `msui-text`, `seq-wash`, `seq-hover`); nothing hard-coded.

## The interaction matrix's cells

Cells this design fills; the others stand as layouts, GA, RR and LX give them.

| | select | extend | add / remove | focus | back | zoom | pan | hide | edit |
|---|---|---|---|---|---|---|---|---|---|
| **a section: the genome ring** | no (its blocks or a range) | — | — | 0 fits; the label menu's "Show only" | undo returns the blocks shown | no: a ring shows its whole axis; opening a chromosome is `section.include` | no; `rotation` by patch | no: `include` by patch | `include` (patch, one undo step), `center`, `rotation` |
| **an element: a ribbon** | click: its drawn run, both ends lit | Shift+click: from the anchor end | Cmd/Ctrl+click | double-click: core's focus (the ring's window, idle unless followed); with `opens`, the pair shown there (Decided, 7) | Alt+Left; undo | — | — | no: a data filter (`count`, `type`) | no: data |
| **an element: a stem** | click: its range | Shift+click | Cmd/Ctrl+click | double-click: as a ribbon's, the stem's block with `opens` | as above | — | — | no | no |
| **an element: a stub** | click: its rect (the stem's range, the partner's hull) | Shift+click | Cmd/Ctrl+click | double-click with `opens`: the pair | as above | — | — | `stubs: none` by patch | no |
| **a block** (a chromosome's arc or label) | click on its label: the whole block (GA) | Shift+click | Cmd/Ctrl+click | label menu: "Show only" (`section.include`) | undo | — | — | no: `include` | no |
| **the selection** | — | — | `rondo:select-partners`: the far ends of its ribbons | Enter: as built | Alt+Left | — | — | — | as built |
| **the keyboard cursor** on the ring | Enter: the stem under it | Shift+arrows: as built | — | `.` and `,`: `rondo:partner`, to the next partner by count, announced | Alt+Left | — | — | — | — |
| **a track or layer** (`rondo:ribbons`) | — | — | — | — | — | — | — | `layer.hidden` | `width`, `stems`, `bundle`, `stubs`, `opens` by `layer.set` (one patch) |
| **a navigator: the square** | no | no | no | one click on a cell: show that pair at once (`navigator.switch: include`); its title or Home: show all | undo | over a ring: drag an edge or corner of the window to grow or shrink the run of blocks shown (Decided 16) | over a ring: drag the window along the diagonal to move the run shown (Decided 16) | as navigators | its params by patch |
| **the view** | — | — | — | — | — | — | — | — | undo returns the blocks shown |
| **a bookmark or link** | — | — | a ribbon's rect or a stem's range as a ref (positions, never the level or the stem) | runs it; a ring nobody follows highlights it | Alt+Left | — | — | — | — |

## Configuration

| Item | Whose | Default | Why the default; the alternatives |
|---|---|---|---|
| `type: 'rondo:ribbons'`, `fallback`, `data` | spec params | `fallback: pairs` (written by `rondoRing`) | the chords a reader without it draws |
| `width` | spec param (plain) | `count` | links counted; a field summed (RR's `sum`) |
| `stems` | spec param (plain) | `level` | Rondo's trumpets; `none` |
| `bundle` | spec param (plain) | 0.7 | ribbons of a stem stay together into the hole; 0 to 1 |
| `stubs` | spec param (plain) | `partner` | "where else does it go"; `none` (core's caret) |
| `opens` | spec param (plain) | absent | a section whose `include` a double-click sets (Decided, 7) |
| the colour channel | spec param (theme slot) | `(field :name system)` of the `a` end, `categories` | a hue per chromosome, safe; `rondo:hues` named by the writer |
| `section.include` (GA's) | spec param | absent (every block of `of`) | the chromosome or pair shown; names, patterns, `{ min_length }` |
| `section.center` (layouts amended) | spec param | the preset (`title`) | `none` under ribbons |
| `navigator.switch` (GA's) | spec param | `none` | `include` on the square over a ring |
| `content: 'rondo:square'`, `of_layer` | spec params | the first pairs-shaped layer of the bound section | the cells read the ring's data |
| `rondo.stemPx` | extension config | 16 CSS px of the layout, never device px | a stem a mouth a reader can point at; the stem level the largest `2^j` × a listed level under it |
| `rondo.stemDepth` | extension config | 0.12 of the hole's radius | room for the flare before ribbons part; any |
| `rondo.maxWidthPx` | extension config | 24 CSS px, the largest stem in view | one px a count for the view; ribbons under `ribbons.minPx` go to the stem's fan |
| `rondo.maxRibbons` | extension config | 1,024 | the items gate and the drawn bound in one: past it the hole is one colour; fans keep totals |
| `rondo.stubShare`, `rondo.stubLabelPx` | extension config | 0.2 of the hole's radius; 24 CSS px of room at the tip | a stub says "elsewhere" short of the centre; a label only where it reads |
| `rondo.hitCellPx` | extension config | 16 CSS px | the hit grid's cell over the hole; candidates few at 1,024 ribbons |
| `rondo.fade` | extension config | 0.6 | the hovered stem's ribbons read over the rest; 0 for none |
| `rondo.flatten` | extension config | 12 segments a cubic | hit within a device px at 720 px; any |
| `ribbons.fillOpacity`, `ribbons.minPx` | `RenderConfig` presets (LX's) | 0.35, 1 device px | one look for ribbons in bridges and rings |
| `pairs.itemsHysteresis` | `RenderConfig` preset (RR's) | 1.25 | no flip at this layer's gate either |
| keys of `rondo:partner` | `ViewerConfig` (`view.keys`) | `.` and `,` | a keyboard way along the links; checked at RX |
| neck order (partners by angle, stubs, the fan last); the twist by low and high edges; width linear, no floor inside a neck | fixed | — | the reasons in the Decision |

Not built: every row (RX; `section.include`, `navigator.switch` in GC3; `center` in GC2).

## Tensions

- **Stems are a lattice, not clusters**: a cluster across a stem edge is two stems, and a resize past a power of two
  changes the stem level. Accepted: exact counts, the same stems on screen, in the SVG and at any DPR; a clustering
  over the lattice (by partner) is a later param if the use gate asks.
- **Two items gates in one figure**: core's chords gate at 4,096, the ribbons at 1,024, so a chromosome's ring with
  2,000 links shows items without the extension and bins with it. Accepted: each draws every item it reads.
- **Repeated hues on a subset** (chr4 and chr12 as a pair): stated in the Decision; `rondo:hues` avoids it.
- **The colour of a ribbon is one end's**: a symmetric Hi-C's `a` end is the upper triangle's, so a ribbon takes the
  earlier chromosome's hue. Accepted (LX's rule); a gradient needs a draw-list form no other layer asks for.
- **24 chromosomes, 8 safe hues**: the safe default repeats; distinct hues are unsafe. Standing; the writer chooses
  (Decided, 5).
- **A Q-value filter at the bins level**: bins count every pair; a filter on an item's field cannot narrow a source's
  own levels. Open (Decided, 8).
- **A cut, not a zoom**: opening chr4 relays the ring with no tween; Rondo animates. Standing until the use gate asks
  (Gate 0, 5).
- **The square over a ring is a picker, not a window**: its window is idle there; the switch is the click. Accepted:
  the same navigator over a grid section is GA's two-axis window.
- **Fans hide partners**: ribbons under a device px or past `rondo.maxRibbons` are drawn as their stem's fan; their
  partners are named only in the stem's card and menu. Accepted: the card lists them by count, the neck stays exact.

## Self-check

The claim: **this design covers what a user does with a Rondo-like genome ring, and what every other reader draws.**

1. **Assume it holds.** A ring of chromosomes with bands and hues; stems as a coarser level of RR's bins, ribbons stem to stem
   by count, stubs for partners not shown; opening a chromosome or a pair by `section.include` from a label, the
   square or a double-click; the keyboard and a screen reader reach the links; a reader without the extension draws
   chords over the same bins.
2. **What it rests on, negated.** (a) One genome, two dozen blocks; plural: two assemblies, 3,000 contigs. (b) Links fit
   a lattice; negated: a cluster across a stem edge, dense data with no gaps, a partner not shown, both ends in one
   stem, a hub with hundreds of partners. (c) The reader has
   the extension; negated: none, an older one, no CN. (d) Zoom is continuous; negated: a cut; inverse: back. (e) The
   user points; negated: the keyboard, a screen reader. (f) Counts are the weight; negated: aligned bases, a score, a
   Q-value filter. (g) Data is held; negated: tiles arriving, pending, too many ribbons. (h) The ring leads nothing;
   symmetric: a linear detail following its window. (i) Colour tells chromosomes apart; negated: a colour-blind reader,
   more chromosomes than hues, a pair eight apart. (j) What is drawn is true; negated: a neck wider than its ribbons,
   a ribbon at the wrong stem, twisted wrongly, drawn twice, a stem that differs between screen and SVG.
3. **The negated world.** (a) sided ends and GA's group gap; runs and unread blocks; (b) stems as a nested level (exact, DPR-free),
   stubs aggregated by partner, links within a stem counted and shown finer, the fan holding a hub's weak partners; (c) the fallback table, the version
   said, the square on core's `axis` picture; (d) `section.include` as one patch, undo as back, a story step's patch;
   (e) the cursor's stem, `rondo:partner`, counts in ARIA; (f) `width` a summed field, the bins filter gap; (g) bins
   until items, gated at what is drawn (`rondo.maxRibbons`), the hatch, fans; (h) core's follow of a ring's window, unchanged; (i)
   the safe default, `rondo:hues` with its substitute, the subset's repeat stated; (j) a neck the exact sum of its
   ribbons and fan (no floor inside), stems as unions of bins in CSS px, the twist by low and high edges, one hit rule.
4. **Both held.** Integrated: every item of (a) to (j) in the Decision. Standing: the Q-value filter at bins (f, Decided
   8); the cut's missing tween (d); distinct safe hues past eight and the subset's repeat (i); a lattice, not clusters (b).

## Consequences

- **`seqx-ext-rondo`** (RX, stage 4, after RG and GC3; on EX and CN): the descriptor and the viewer extension; the
  layer types, picture, scheme, commands, menu entries, binding; `rondoRing`; `_spec/README.md` (the stems, ribbons,
  stubs, the fallback table, its Configuration table); its `noRestrictedImports`; its inventory row and size row; a
  changelog entry file in `.v6/changelogs/` (an addition). Tests (quick lane, `_test/`, exports `null`):
  - drawing against the data: every bin pair in view in exactly one ribbon, stub or fan; a neck's width `k ×` its
    stem's count and the sum of its ribbons, stubs and fan (a hub stem of 300 partners included); a stem's count the pairs with an end in its bins; a ribbon's ends inside its stems' mouths; stubs only
    for ends with no place, one a stem a partner; at the fit and across a tile edge, a reversed block, two assemblies;
  - property tests: stems a partition of the non-empty bins, the same at DPR 1, 2 and in the SVG; widths additive;
    items asked only at or under `rondo.maxRibbons` (no fan below the gate); the twist by
    geometry equal to GA's rule over every direction combination; the closed hit equal to the drawn polygon's topmost
    cover (sampled); held and region pairs give the same stems;
  - the fallback: the same spec drawn without the extension (canvas and SVG alike), the issue once, the version said;
    `viewer.extend` and its undo redraw (EX's key);
  - EC-108, EC-169, EC-171, EC-206, EC-210, EC-218, EC-225, EC-235 named by tests.
  - Performance, each with its path: a ring of 23 blocks and 2,000 links laid out (stems and ribbons) under 4 ms and
    drawn under 8 ms (quick representative); 12,000 links at the fit: layout under 16 ms once, a frame under 8 ms, a
    hover under 2 ms (heavy, `itWithin`); the draw and hover loops at 0 bytes warm (the allocation tier, the plan's
    "a ribbon frame"); the ring in the browser at DPR 2 (M5).
- **`seqx-core`**: `section.include` and its resolve (GC1's bump, GA's text); `section.center` (GC2); the `system`
  field's domain before `include` (GC2); the `extension` entry exports RR's `PairsReader`, the bins' types, the level
  rule (`pairsLevels` with `binPx` and `gate` arguments) and `coarsenBins`, GA's `placesOf`, and the polar warp (`PolarWarp`, already
  listed); the `LayerType` hooks `describe`, `summary`, `announce` (EX amended); `navigator.switch`, the square bound
  to one section on both axes, its lattice the assembly's (GC3); `rings.center` promoted to `section.center` (GC2).
- **`seqx-viewer`**: the `section.include` command (`patch`) and the label menu's "Show only" and "Show all" (GC3);
  extension bindings scoped to a layer type's hits (EX amended, Decided, 7); the hover card, ARIA and cursor calling the
  hooks.
- **`seqx-examples`** (GD's): the GRCh38 assembly by inline lengths, the seeded links generator (held and BEDPE),
  `syntheticBands`, the recorded cytobands once their licence is agreed; the showcase spec through `rondoRing`.
- **Storybook**: a Design page "Genome ring"; a story per default above, light and dark, canvas and SVG, each with its
  headless twin on a fake clock (bins landing, then items in chr4's ring); the same story with the extension removed.
- **`_spec` sections that change**: core Language (section `include`, `center`; navigator `switch`), Resolved view,
  Render (Rings: the title; Navigators: the switch; the extension entry); viewer Commands, Menus, Bindings.
- **Other designs amended**: [extension-openings](extension-openings.md) (the entry's exports;
  extension bindings scoped to their layer's hits; Part 1's "Enough for the genome ring?" replaced by the text below);
  [layouts](layouts.md) (Rings: `center` a spec param; the zoomed ring a section's `include`);
  [themes](themes.md) (an extension's unsafe scheme with its substitute, `rondo:hues` as the example);
  [features](features.md) (the genome ring row planned with RX); [navigators](navigators.md)
  through GA (the switch).

> **Enough for the genome ring** (the fifth slice's `rondo`): `rondo:ribbons` (a `scope: 'section'` layer over core's
> pairs level, Gate 0, 4, on `rings` through the warp), `rondo:cells` and the `rondo:square` picture on core's square,
> `rondo:hues`, two commands and a binding; the ring zoom is core's `section.include` (a cut), the square's switch
> core's. The entry exports `PairsReader`, the pairs bins' types, `pairsLevels`, `coarsenBins` and `placesOf` beside
> the warp. A `LayerType` may give `describe(hit, lay)` (the hover card's lines), `summary(lay)` (its group's ARIA
> description) and `announce(lay, system, position)` (what the keyboard cursor adds there); the viewer calls them, its
> own text after. An extension binding may be scoped to its own layer type's hits, and then wins over the grammar's
> default there. One `seq.view.layer.extension` a namespace and view says every part drawn as its fallback (layers,
> pictures, schemes).

**Built** (RX, merged 2026-11-23): the package (`rondo` 28 kB, its row 32 kB), `section.center` (language 0.10.0),
`PairsReader`, the hooks and scoped bindings, the fallback, the showcase data and stories; a hit slop for thin ribbons
(`rondo.hitSlopPx` 4); stub and block labels without overlap; a band's hover lighting only its links. Not built, RX2's:
`rondo:cells` and the `rondo:square` picture, element menus and the stem menu, `.` and `,`, the extension's config read
by its layer types, a failed source's status; and Decided 16 to 18.

**Built by RX2** (merged 2026-11-23): Decided 16 to 18 (the square a window: move, edge and corner resize, one patch
a gesture, the keyboard; Cmd/Ctrl+click sets outlined cell by cell; one click on a cell shows its pair), the square's
cells (`LayerType.square`, `rondo:cells` by count on a square-root ramp), element menus by layer type and the stem menu,
`.` and `,` walking partners, `settingOf` (an extension's settings, recomposed on a late extend), a failed source read
failed (the store's `unopened`), on one chromosome or a pair lighter stems (`focusStems` 0.35) and stubs under
`rondo.stubMin` folded unlabelled, unfolding on hover. The square, its gestures, cells and lattice load as one lazy
chunk. Not built: Shift+F10 at a stem; `rondo:cells` lazy (needs `withViewKind` on the extension entry).

**Built by RX3** (merged 2026-11-23): kinds of links as layers (`rondoRing({ kinds })`, the `(item :type …)` base
filter), a ribbons layer's colour (fixed, drawn opaque; or a field per link, exact per bin where held bins carry the
field, else ink and said), `opacity`, a legend entry per layer with its scale ("1 px = n links"), widths per layer by
default and `scale: 'section'` to share one (Decided 19), the card naming a link's fields and a bin's count per layer,
invented translocations named "(made up)". Core: filters applied to pairs read by region (items; bins read by region
drawn once and said), `strandKnown` for unstranded ends, settings composed on a late extend.

## Decided by the author

2026-11-21: every gap as recommended (the author on the drafts: "generally reasonable and we can fold in the reviews";
on the whole-designs review's twelve recommendations: "yes with the recommendations"; on the fold-in: "yes"). An item
the author decided by name carries its own date.

1. **A stem is a bin of a nested coarser level** (`2^j` × a listed level under `rondo.stemPx` in CSS px: exact,
   cached, the same on screen and in the SVG), not a clustering by gaps: **decided yes**.
2. **Width one px a count for the view, no floor inside a neck**: ribbons under a device px go to the stem's fan, so a
   neck is always the sum (the cap and additivity hold; the floor gives): **decided yes**.
3. **The ring zoom as `section.include`** (a core spec param and `patch` command, in GA), in place by default, a
   second ring below as the builder's option; back by undo: **decided yes**.
4. **The square bound to one section on both axes, with `navigator.switch: include`** (core, GA's), RN drawing
   only its cells: **decided yes**.
5. **The layer's default palette the safe `categories`** (repeating past 8); `rondo:hues` (24 distinct, not safe,
   substitute `categories`) named by the showcase spec: **decided yes**.
6. **Stubs inward**, one a stem a partner, in the partner's hue (the `system` domain before `include`): **decided yes**.
7. **`opens`**: a double-click on a ribbon, stem or stub opens it in the named section, through an extension binding
   scoped to its layer's hits that wins over the grammar's focus there (extension-openings amended): **decided yes** (else menus
   and the square only).
8. **A data filter on an item field rebuilds held and summarized pyramids**; a source's own levels filter items only,
   said (RR's): **decided yes**.
9. **`rings.center` promoted to a section spec param `center`** (`title`, `none`) so ribbons own the hole: **decided yes**.
10. **Cytobands**: ask UCSC's licence for `cytoBandIdeo.txt.gz`, then fetch it once into `seqx-examples` with its
    script and `sources.json` row (not the current Seq\*'s unlicensed copy); NCBI's public-domain ideogram if UCSC's
    terms refuse: **decided yes**.
11. **`rondo.maxRibbons` 1,024 as both the items gate and the drawn bound** (`pairsLevels`' `gate`, RR's),
    fans only for ribbons under a device px or past the bound at the bins level: **decided yes**.
12. **Two entries** (`@molstar/seqx-ext-rondo` headless, `/viewer` with commands and the binding), importing only the
    `extension` entry: **decided yes**.
13. **Ribbons on a `tracks` row too** (filled arcs between stems), text and grids refused: **decided yes**.
14. **Three optional `LayerType` hooks** (`describe`, `summary`, `announce`) for the hover card, ARIA and the cursor
    (EX amended): **decided yes**.
15. **The square's lattice the assembly's**, washed cells clickable under the switch (GA's): **decided yes**.
16. **The square over a ring is a real window** (the author, 2026-11-23, on trying it: "i dont understand how this
    control works ... its impossible to drag or resize the box. what works is clicking once to select a grid point, and
    then the 2nd click enlarges that section... very confusing"; then "a real window, as recommended"): its box is the
    blocks the ring shows (`section.include` as a run of blocks on the diagonal, or a pair off it); dragging the box
    moves the run along the diagonal, dragging an edge or corner grows or shrinks it (snapped to whole blocks), one
    click on a cell shows that pair at once (never a select first), its title or Home shows all; each gesture one patch
    and one undo step; the cells coloured by count (`rondo:cells`) so the window is read against the links. Replaces
    "idle over a ring" for zoom and pan in the matrix above.
17. **On one chromosome or a pair, the information first** (the author, 2026-11-23, "yes, fold the stem and stub idea
    into RX2"): stems drawn lighter so the ribbons read over them; stubs under `rondo.stubMin` (2) folded into one "other
    chromosomes · n" stub a stem, the card listing the partners, a neck still the sum of its parts.
18. **The square selects any set** (the author, 2026-11-23, "yes, fold it into RX2"): Cmd/Ctrl+click on a cell or a
    label adds or removes chromosomes from `section.include`, a label's click shows only it; a contiguous run is one
    draggable box, a scattered set outlined cell by cell (its diagonal and cross cells) without grips; the cells
    coloured by count with self-connections on the diagonal (`rondo:cells`).
19. **Widths per layer, said** (the author, 2026-11-23, "Per layer, said"): each kind's widths on its own scale, each
    legend entry stating it ("1 px = n links"), and `scale: 'section'` for a host that wants widths comparable across
    kinds.
