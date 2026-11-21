# Genome coordinates: assemblies, chromosomes as blocks, a direction per block, navigating across systems

Lane GA's design (the fifth slice, stage 0; agreed by the author 2026-11-21 and promoted from its draft in `execution/fifth-slice/drafts/`). It amends
[axis](axis.md), [navigators](navigators.md), [layouts](layouts.md),
[addressing](addressing.md) and [extension-openings](extension-openings.md) (below,
Consequences). Built on Gate 0's answers as recommended ([plan](../execution/fifth-slice/plan.md)): the navigator over several systems is
core's (8); a ring's zoom is a cut (5: `section.include`, or a second section); the ring extension is a package of its own on the
public surface (3), so what it needs from coordinates is here, its pictures are RN's.

## Intent

Version 1's fourth milestone, Genomes: "two assemblies on one axis", "a genome ring like Rondo's", synteny with an
inversion, a liftover, open on a gene by name ([version-1](version-1.md), [examples](examples.md),
Genomes). Four things the agreed designs left open, each by name:

- **A direction per piece** ([axis](axis.md), Tensions: "an inversion between genomes on one axis needs a
  direction per piece … Open until Genomes").
- **Two assemblies on one axis or ring** (EC-210; [features](features.md), "Two assemblies on one axis",
  deferred; [layouts](layouts.md), Blocks: "Two genomes on one ring wait for the Genomes milestone").
- **Switching chromosomes and a navigator over several systems** (EC-184; [navigators](navigators.md),
  "switching chromosomes on a genome-wide ideogram: refused for now", its self-check's standing residue;
  [extension-openings](extension-openings.md), Decided 7: "designed in the fifth slice").
- **Scaled segments and a zoomed ring** ([layouts](layouts.md), Tensions and Self-check (f): "the
  Genomes milestone decides").

And under them, what nothing yet names: an **assembly** (its sequences, their names and aliases, an HLA contig's name,
EC-186 to EC-189, EC-200), and what a position, a location, an address and a selection are when blocks run in opposite
directions, read by paths, hover, follow, history and export.

## Precedents

| Precedent | Taken | Not taken |
|---|---|---|
| The current Seq\*: blocks with a `direction` each (1 or −1), a layout map, never new data (EC-31, EC-32); the JBrowse reply's "keep `direction` on sections and blocks" | a direction per block; reversing never changes a position | complement tied to it (the axis design's `complement: auto` reads per block) |
| UCSC: a chromosome at a time, the position box (`chr17:7,668,402-7,687,550`, 1-based, commas), chromAlias tables, the "reverse" button, multi-region mode | the locus text, aliases by table, one chromosome at a time as a configuration (`navigate.acrossBlocks: stop`) | a navigation that reloads; a reverse of the whole window only |
| IGV: a genome dropdown and "All", the chromosome ideogram with its red box above the tracks, a click on "All" choosing a chromosome, `chr` prefixes matched loosely, the locus box with gene names | the whole genome as blocks, a click on a chromosome a jump, the ideogram of the current chromosome under the genome's (the zoom strip, below), the `chr` prefix rule | the whole-genome view as a separate mode with its own data rules (here the same section at the fit) |
| JBrowse 2: multi-region linear views (several displayed regions end to end, panning across them), per-region "reverse", refName aliases, synteny views with a view per assembly, "horizontally flip" | blocks end to end with panning across them (as built), a direction per block, an alias table per assembly, each assembly its own section in a synteny view | dynamic region lists edited by drag (here a block list is the spec's) |
| Circos and pyGenomeViz: ideograms as arcs with gaps, reversed ideograms (`chromosomes_reverse`), scaled segments (`chromosomes_scale`, "chr1:0.25"), a gap between karyotypes | a direction and a scale per block, a wider gap between assemblies | static only; scales by regex over names |
| Rondo (research-genome-ring, third slice): chromosomes on a ring, one chromosome as its own ring, the square navigator | the ring of blocks (built), one chromosome's ring as a cut (Gate 0, 5), the square's coordinates (a two-axis window over blocks × blocks) | the square's picture and the zoom's transition (RN's) |
| This codebase: blocks (one strip a block, each with its system and axis, `arrange`, `copies`), a section `of` a complex deriving a block a unit with overrides by id, the navigator's path-bound window and zoom strip, the `(sequence :name …)` form for names that are not ids, coverage aliases in the store | an assembly derives a block a sequence as a complex derives a block a unit; the navigator's window carried along a path; aliases resolved where the store already resolves them | a new layout or a new transform: none is needed |

### Interactions in the precedents

| What a user does there | Here |
|---|---|
| choose a chromosome from a list or the whole-genome ideogram (IGV, UCSC, JBrowse) | built: a click on a chromosome in a navigator over the assembly, a jump; `block.next`, `block.previous` (Alt+PageDown, Alt+PageUp); a host's list calls `view.focus` with the block's address |
| type a locus `chr17:7,668,402-7,687,550` or `17:7668402` (all) | built: the `(locus :text …)` form, aliases and the `chr` prefix matched; a host's box calls focus with it |
| type a gene name (IGV, UCSC) | built: the locus form falls back to GF's `findNames` (a GFF3 read whole, EC-228) when no sequence matches |
| see the whole genome, double-click a chromosome to open it (IGV "All") | built: the section at the fit; double-click on a block's label or a navigator's block fits it |
| pan past a chromosome's end into the next (JBrowse) or stop at it (IGV, UCSC) | built: `navigate.acrossBlocks: cross` (as built) or `stop` |
| reverse a region or a chromosome (JBrowse, Circos) | built: `block.direction`, `block.reverse` (a patch, one undo step) |
| reverse the whole view (UCSC, JBrowse flip) | built as before: `section.reverse` (each block flips; the block order is the spec's) |
| zoom one chromosome on a genome ring (Circos scale, Rondo zoom) | built: `block.scale` (Circos's scale); the zoom into one chromosome a cut by `section.include` (Gate 0, 5) |
| two genomes on one circle (Circos karyotypes, the JBrowse reply's peach and grape) | built: `section.of` a list of assemblies, `blocks.groupGap`, a label a group |
| a square navigator of chromosome pairs (Rondo) | built in core: a two-axis window over blocks × blocks; its picture RN's |
| drag a region across chromosomes to select it (JBrowse) | built: shown positions of every block between the ends |
| reorder chromosomes by drag (JBrowse's region list, Circos's order) | refused: the order is the spec's (`order`, a list by patch); a drag-reorder of blocks has no second use yet |
| split a chromosome into reversed pieces by gesture (JBrowse) | refused: the writer spells blocks (the builder's `colinearBlocks`, below); a gesture that edits the block list is a figure edit with no precedent in our grammar |

## Decision

### Whose each choice is

By [axis](axis.md)'s rule (a spec param, a `RenderConfig` preset, a `SpecBounds` or `ResolveConfig`
bound, a `ViewerConfig` item, a port or registry entry, an option of a call, or fixed with its reason). No new port:
assemblies, blocks and navigation are pure core over the resolved view; names and lengths arrive through the sources'
coverage, which exists.

### An assembly

- **A new data kind, `assembly`**: a genome's sequences, named, ordered, with their aliases. Params: `name` (the
  accession or the assembly's name, `GRCh38`; EC-187), `aliases` (other names of the assembly, `hg38`), `source` (a
  ref: the source that lists its sequences with their lengths: a chrom.sizes table (`chrom-sizes`, GF's BED entry), a
  FASTA read whole, a host's adapter; a 2bit or an indexed FASTA only through a host's adapter, their runtimes being
  designed only, GF), `sequences` (an inline table: `name`, `length`, `aliases`, `cyclic`; instead of or beside `source`),
  `alias_source` (a ref: a chromAlias-like table, a column a naming), `names` (the naming printed: a column of
  `alias_source`; absent: the names `source` gives), `alphabet` (`dna`), `order`, `include`, `match` (below). Why a
  kind: a genome's twenty-five to several hundred sequences written as `sequence` nodes is a manifest by hand (the JBrowse
  reply's point 2); a complex already derives its units this way. **The kind's name** is kept though mmCIF's source
  option is also called `assembly` (a biological assembly's id): that is a param's value on a structure source, not a
  kind, the language's kinds and params are separate names, and `source.assembly` (built) already means a genome
  assembly; the docs say "genome assembly" (a gap asks whether to rename it `genome`).
- **Its sequences are derived**, a system each, as a complex's units are: from `sequences`, else from its source's
  coverage (names, lengths, aliases), with the letters read from `source` when it serves them (a FASTA read whole, a
  host's adapter by region). A **`sequence` child whose `name` is a derived one's** overrides that one (`cyclic: true`
  on a mitochondrion, a `label`); the rest stay derived. Why: one chromosome's setting must not freeze the list (the
  complex's rule, [molecules](molecules.md)).
- **Every sequence has a length**: a derived sequence with none (a source listing names only) makes the assembly
  unavailable (`seq.resolve.assembly.lengths`, error, naming the sequences and the remedy: inline `sequences` or a
  `chrom-sizes` source), never drawn with a guessed length. Why: `include: { min_length }`, `order: length`, the gaps'
  share and the closing run all read lengths.
- **Static, after open**: the derived sequences enter the static model as the sources' lengths already do, an input of
  resolve (`ResolveOptions` gains the assemblies' coverage beside `lengths`, hashed with it), so the resolved view stays
  a pure function and never holds a dynamic node ([dynamic-nodes](dynamic-nodes.md)). Until the source
  has opened, a section of the assembly is pending (EC-112), drawn as the axis design's pending axis with the
  assembly's name. Why not dynamic nodes: blocks, scroll groups and follows are the resolved view's, which must not read
  the overlay.
- **Addresses**: a derived sequence is `/data/hg38/(sequence :name chr17)`, always the form, never an id step (an id
  step would be a second address of one system: no elision, [addressing](addressing.md)). `:name`
  reads any alias and prints the assembly's `names`; `HLA-A*01:01:01:01` is a quoted string there (EC-200), the id
  pattern unchanged.
- **A derived block is not a node and has no id**: its address is a keyed form under its section,
  `/views/main/genome/(block :name "HLA-A*01:01:01:01")` (the sequence's printed name), as an alignment's rows are
  `(row :id …)`; `ResolvedBlock.address` holds that form, never `${section}/${id}` (today's `blockAxes` builds it from
  `n.id`, which would throw on a name outside `IdPattern`, core `repr/tree/node.ts`). A `block` child that overrides a
  derived block is matched by its `of` (the sequence's address), not by its id, so its id is free (`hla_a`); a name in
  the id pattern may still be the override's id, as shorthand. Why: names that are not ids go through forms (EC-200,
  the addressing design's rule), and an override must reach every sequence. Amends addressing.md (the `block` form) and
  the molecules design's override-by-id rule for assemblies only (a unit's ids are its chain ids, always in the pattern
  after the structure reader's check).
- **Names match through aliases, then the prefix** (EC-188, EC-189): a name a source or an address gives is looked up
  in the assembly's names, its aliases (the table's, then each sequence's), then, with `match: chr` (the default), with
  a leading `chr` added or removed (`1` is `chr1`, `M` is `chrM`; `MT` needs an alias). An exact name always wins over
  the prefix rule; an alias naming two sequences is `seq.resolve.assembly.alias` (warning) and dropped, never guessed;
  so is an alias joining two sequences of different known lengths. Why `chr` by default: UCSC's and Ensembl's names for
  one assembly differ by the prefix in the common case, and a mismatch draws nothing with no word why; `exact` is the
  writer's for an assembly where `1` and `chr1` differ. **Not every pair is an alias**: hg19's `chrM` (NC_001807) and
  GRCh37's `MT` (the rCRS, NC_012920) are different sequences of different lengths; the prefix rule never joins `M` and
  `MT`, and the length check refuses an alias table that does. The builder's docs say so beside `alias_source`.
- **A source on an assembly** (`source.assembly`, built, EC-187): checked against the assembly's `name` and `aliases`;
  another is `seq.resolve.assembly.mismatch` (warning), the data still drawn (a reader may know better).
- **Annotations over an assembly**: `target` (ranges, values, pairs, grid) may name an assembly; each row is on the
  sequence its `target` field names, through the names rule. A row naming none is counted and said once
  (`seq.resolve.assembly.name`, info, with the first unknown names), never dropped silently. Why: a BED or a bigWig is
  one file over every chromosome.
- **`order`** (spec param): `natural` (the default), `length` (longest first: a draft assembly of scaffolds), `source`
  (as listed), or a list of names. `natural` ranks, on the name less a leading `chr` (case-insensitive), so UCSC's and
  Ensembl's unprefixed names sort alike: (1) a whole number, ascending (`1` … `22`; `2` before `10`); (2) `X`, then `Y`
  (then `W`, `Z`: birds); (3) `M`, then `MT`; (4) the rest in the source's order (unplaced, alternates, `HLA-…`).
  Why `natural`: the source's order differs by format for one assembly (a FASTA index, a 2bit, a chrom.sizes sorted by
  length), so the same genome would draw in three orders; a sort over names is a rule, not an inference from data.
- **`include`** (spec param): every sequence (the default), a list of names and patterns (`^chr[0-9XYM]+$`), or
  `{ names?, min_length? }` (patterns and a length floor: `{ min_length: 1000000 }` keeps a draft's chromosomes and drops
  its scaffolds). Why all by default: which contigs matter is the figure's (a library, not a viewer); the examples and
  the builder spell the primary set.
- **Bounded, for assemblies only**: `ResolveConfig.maxAssemblyBlocks` (4,096: GRCh38 with its alternates and decoys is
  about 3,400) bounds the blocks a section of assemblies derives after `include`; complexes keep their own
  `maxUnits` (65,536, refused whole). Past it the section draws the first 4,096 in `order` and **one closing run** for
  the rest, labelled with what is missing ("+96,312 sequences not shown"), never silently; `seq.view.section.blocks`
  (warning) names the count and the remedy (`include: { min_length }`, `order: length`). Why not refused whole: a draft
  assembly with 100,000 scaffolds still has a figure in its largest, and the run says what is left out.
- **Resolve's cost at this scale** (100,000 sequences listed, a few thousand derived): linear in the sequences, each
  once. `blockAxes`' copy groups go through a `Map` by key (today a linear `groups.find` per block: quadratic); an
  assembly's coverage enters the model's hash as one entry per source, its content hashed once when the source opens
  (not by source and version, which would miss a changed file at one version; not each length sorted and joined at each
  resolve as `lengthsKey` does today); the names rule is one `Map` of names and aliases per assembly, built once per
  coverage. A quick representative: an assembly of 100,000 sequences resolved with `include` to 25 under 50 ms.
- **Maps between two assemblies** (an alignment's runs through `map.source`: the `seq` runtime's sided pairs, a host's
  chain or PAF adapter, GF; LX's figures): a `map` whose `from` and `to` are
  assemblies is **one map over every sequence pair**, its runs carrying their two sequences; a path between two
  sequences of them takes that map's runs for the pair, as one step. Why: a chain file is one map, not one per
  chromosome pair; held so, EC-207's "several maps between one pair of systems" never makes a path ambiguous.
  By region, as any large map ([building-blocks](building-blocks.md)): a projection asks for the runs
  of the region it carries (EC-211). **This amends** [maps](maps.md) and building-blocks.md: a map's
  step is no longer only system to system; a map between assemblies is a step from any sequence of one to any of the
  other, its runs keyed by the pair, so path search counts it once and the shortest-path rule (EC-22) and the ambiguity
  rule (EC-23, EC-24) read it as one map. A liftover (EC-212) is projection through it (GF's `map.source`, LX's
  figure).

### Loaded when used (the author, 2026-11-21: "Lazy genome resolve")

Assembly resolve, the names rule, derived sequences and the `(locus …)` and `(block :name …)` forms are one chunk,
loaded the first time a spec holds an `assembly` (or a query names a locus or block), as navigators load in
`renderSvg`. The language's kinds and checks stay in the reading path (a reader checks every spec); what reads an
assembly does not. Until the chunk arrives the spec is pending (EC-169's state), never wrong. Why: a viewer that never
opens a genome pays nothing for it, and the minimal and headless entries keep their bounds.

### A section of an assembly: chromosomes as blocks

- **`section.of` takes an assembly**: a block a sequence (the layouts design's blocks, as a complex gives a block a
  unit), in the assembly's `order`, after `include`. Explicit `block` children still win; a `block` child whose id is a
  sequence's name overrides that block only (its `direction`, `scale`, `label`, `range`). Why blocks: they already
  carry a system, an axis and layers each, a strip each in every layout, rings included (built: the ring of blocks,
  `rings-fixtures`); an assembly adds no layout.
- **`section.include`** (spec param, on a section of an assembly or a list of them): which of the derived blocks the
  section shows, in the assembly's `include` grammar (names, patterns, `{ min_length }`), intersected with the
  assembly's own; absent, every block. A command `section.include` (`patch`, one undo step) sets or removes it; a
  block label's menu has "Show only <name>" and "Show all". So a ring of one chromosome or of a pair is the genome's
  section with one param changed, a cut (stories), and a section of one chromosome can switch. Why a section param and
  not new `block` nodes: one param to patch and undo, the assembly's grammar, the hue domain and `order` untouched (the
  ring's zoom, RN).
- **`section.of` takes a list of assemblies**: their blocks in turn, each assembly a **group** (EC-210). A wider gap
  between groups (`blocks.groupGap`, `RenderConfig`: 3 × the block gap in tracks, `rings.groupGap` 4° on a ring), and a
  group label (the assembly's name) on a bracket under a row or an arc outside a ring's labels. Why one section: the
  peach-and-grape ring and a genome-vs-genome ribbon figure need both on one axis; two sections (one an assembly) is
  the other figure, and both are the writer's.
- **The axis of a section of several systems** stays units on no system (as built); a position on it is no block's.
  Every hit, cursor, selection part, link and bookmark carries the block's system and a position on it.
- **Across blocks, navigation** (`navigate.acrossBlocks`, `ViewerConfig`): `cross` (the default, as built: one viewport
  over the blocks end to end, a pan carries into the next chromosome across the gap, JBrowse's) or `stop` (a pan, a
  zoom and a fit clamp to the block holding the window's centre, IGV's and UCSC's one chromosome at a time; a jump to
  another block moves there). Why both: two real uses (the precedents split); `cross` keeps the built behavior.
- **The whole genome at the fit view** is the same section, not a mode: every block at its share of what the gaps
  leave. Its layers read their coarse levels by block (each block a system the store reads by region; the navigator's
  rule: never a chromosome at level 1).
- **Gaps at scale, one rule for rows and rings**: the gaps in view share at most `gapShare` of the length (a quarter),
  each at most the layout's gap (tracks: 12 px, as built; rings: `rings.blockGap`, 1°), so `gap = min(max gap,
  gapShare × length / gaps in view)` with no floor: at 3,000 blocks a ring's gaps are 0.03° each and three quarters of
  the circle stays positions. Amends the rings layout: today each block takes a fixed 1° and the arc left is
  `max(2π − n·gap, 1e-6)` (`render/layout/rings.ts`), so past 360 blocks there is no ring. A group gap (two assemblies,
  below) takes its share first, fixed.
- **Small blocks run together**: consecutive blocks each narrower than `blocks.runPx` (2 CSS px, `RenderConfig`) at the
  current scale are drawn as one **run**: one band over their joined extent, no gaps inside, a hairline at each end,
  labelled when it has room ("2,841 contigs"). **The joint rule with RR** (one text in both designs): inside a
  run a block's layers read and draw its data at its place while the block is at least one device pixel wide
  (`blocks.readPx`, 1 device px); a block under one device pixel reads nothing at that zoom, and its stretch of the run
  shows the run's band alone, its hover naming the block and "not read at this zoom". Pairs with an end on an unread
  block are still read from their other end's tile and drawn at the block's place (`placesOf`); a pair with both ends
  on unread blocks is not drawn, and the run's hover counts it. Runs and the read set are recomputed on a zoom (a
  prefix sum over blocks, O(blocks)), never on a pan. A run is a drawing of gaps and labels, not a merge of data. A hit
  in a run names the block under the point and the run's count (its hover: "chrUn_KI270302v1, one of 2,841"); a click
  selects that block; a double-click fits the run, where its blocks part again. Why: below
  two px a block's gap and label are noise and a hairline each is a comb; data needs one device pixel to be seen, and
  under it a block's data is more requests than pixels (3,000 contigs × layers), while its pairs, the data a genome view
  reads across blocks, survive through their other end. The same rule on rings, measured at the outer radius.

### A direction per block (closing the axis design's tension)

- **`block.direction`** (spec param): `section` (the default: the section's resolved direction, as built) or
  `reversed` (against the section's). Relative, not absolute, so `section.reverse` flips one param and every block with
  it. **No `{ along }` on a block**: the maps between assemblies arrive by region (EC-211), so a direction read from
  their majority run would change as regions arrive and turn a block over under the reader, the reason
  `pieces: { along }` is refused too; `colinearBlocks` (below) reads the map once and writes `reversed` where a run is.
- **`section.reverse` on a section of blocks** flips every block as built and, as built, keeps the block order (chr1
  stays left). JBrowse's flip mirrors the whole view, its regions' order too; a gap asks whether to follow it.
- **A block is the piece that turns over**: an inversion between genomes is drawn by cutting B's chromosome into
  blocks at the inversion's edges and reversing the inverted one; each block one direction. The axis design's
  one-system axis (blocks on one system end to end with breaks) holds only while its blocks share a direction; blocks
  on one system with mixed directions lay out as blocks on no system (a gap between them, each its own strip). Why:
  every strip already has its own `AxisTransform` and `dir`; a direction per piece inside one transform would change
  the transform, every lookup and the pieces' shape for one figure, while blocks give it with no new code path.
- **Overlapping blocks on one system** (a duplication shown twice to line up with two partners): allowed; a Location
  on the system lights every block that shows it (the merged-copies rule), a hit names its block. Why: a synteny figure
  shows a duplicated segment once per partner.
- **The builder's `colinearBlocks(map, of, along, { min_length })`** writes the blocks of one sequence in another's
  order from a held map's runs (a block a run of at least `min_length`, 100 kb by default, shorter runs merged into
  their neighbour; its direction `reversed` where the run is), so a figure "B in A's order" is explicit data, exported
  as drawn. Why a helper and not `pieces: { along: map }`: an axis derived from a map read by region
  would change length as regions arrive (EC-211); the helper reads the map once, the spec keeps the result.
- **What reads a reversed block** (each the axis design's rule, per strip): the ruler counts down within it; its label
  carries the direction mark (`chr4 ←`, by `axis.directionMark`); `complement: auto` is the section's, times the
  block's, times the run's; on a ring it runs counterclockwise within its arc; navigation stays in display order.
- **What this closes, and what stays open**: closed: an inversion between genomes drawn on one axis or ring (EC-206, EC-210), by a
  reversed block. Open: a direction per piece inside one block, an inversion reversed in place on one system's axis
  with a zero-width break instead of a gap (a collapsed, exons-only axis turning over mid-way); the axis design's
  tension is narrowed to that, not closed. Many small inversions need neither: LX's ribbons draw them twisted in each
  genome's own order (below), with no blocks cut.
- **`block.reverse`** (command, `patch`, one undo step): flips one block's `direction`. Reached from the block label's menu (and Shift+F10 on
  it); no key of its own (Shift+R stays the section's). Why no key: Alt+Shift+R types a character on macOS (EC-155).

### Scaled segments and the zoomed ring (layouts' open tension, decided)

- **`block.scale`** (spec param, 1): a block's display units are its length times its scale, in every layout (rows and
  rings), so a region can take more of the row or the circle (Circos's `chromosomes_scale`). Its label says it
  (`chr17 ×10`), its ruler counts its own positions at its own density, and the hover and the ARIA text read positions,
  never units. Why: a block already owns its strip's transform, so a scale is one factor at its offset; unlabelled, a
  scaled block would lie about length.
- **A zoomed ring is a cut** (Gate 0, 5): the genome's ring with `section.include` patched to one chromosome or a pair
  (in place), or a second ring section of it (a section `of` its sequence, a block with a `range`, or the builder's
  ring below); reached by the label's menu, the square's switch, a story step or a host's layout; the ring's key never
  holds a viewport. An animated ring zoom waits for the use gate to ask for it (a story channel, RN). Why: the ring's
  laid-out labels stay placed once per layout, and its window stays the linear detail's.

### Positions, locations, addresses and selections when blocks reverse

- **A position is the system's**, 0-based, whatever a block's direction or scale; **display order** is the section's
  offsets. Nothing stored or exchanged (a Location, an address, a bookmark, a link, a saved viewport, an export's
  `data-` attributes) is in display order or in units.
- **A Location** is parts by system (built: `{ parts: { system, ranges }[] }`), its ranges ascending on the system.
  A selection across blocks of several systems is one Location of several parts.
- **An address** names a sequence by `:name` and its positions by the forms as before; **`(locus :text …)`** (a new
  form under an assembly or a section of one) reads the precedents' locus text: `chr17`, `chr17:7668402`,
  `chr17:7,668,402-7,687,550` (1-based, inclusive, commas), its name through the names rule; a name with colons (an HLA
  contig) is matched longest first (the whole text, then up to the last `:`). It prints as the canonical
  `(sequence :name …)` and range form; an unknown name is `seq.query.locus.name`. Why a form: a link, a bookmark and a
  host's box carry the same text the reader typed. Amends [locations-and-queries](locations-and-queries.md)
  (a generator form) and addressing.md (its table).
- **Selection in display order**: a drag, Shift+click and Shift+arrows extend over the display offsets between the
  anchor and the pointer: the shown positions of each block between (whole blocks in the middle, `selection.acrossBreaks:
  shown` as for collapsed introns), each block's part on its own system. The anchor stays a position on its system.
- **The keyboard cursor crosses blocks** (amends EC-154's built "does not cross"): at a block's last position in
  display order a step goes to the next block's first, announced with its name ("chr18, 1"); `block.next`,
  `block.previous` jump a block (Alt+PageDown, Alt+PageUp; a navigator's slider PageDown and PageUp, ARIA's large
  step). Why: a genome read with the keyboard must reach the next chromosome; the old rule kept a step from landing
  on another chain's positions, which crossing to its first position does not do.
- **Hover and paths**: a hit names its block, system and position; projection carries it along paths as before (a
  synteny partner lit through the assemblies' map, its run's direction its own). Direction never enters a path: maps
  carry their runs' directions (EC-206); blocks only lay them out.
- **Follow**: a section follows another whose block list is equal (systems, ranges, directions, scales in order) by
  identity in units, the scroll group's key; this lifts `seq.render.follow.blocks` for them (two sections of one
  assembly, a ring of the genome and its linear detail). A follow between unequal block lists (assembly B following A
  through their map) stays refused, its reason now in the issue: which block's map a window carries is a figure's
  choice no rule makes; synteny sections move on their own (LX), or through a navigator.
- **The view history and saved state**: a section of blocks saves and steps its viewport as the block's system and the
  position at each end (`{ system, position }`), never units, so a step survives a change of `include`, `order` or
  `scale`; a step over a changed block list is restored by position under the one rule (axis, history), clamped to what
  is drawn (an end on a block no longer shown goes to the nearest shown). Axis signatures gain the blocks' directions
  and scales.
- **Export**: each block a group with `data-system`, `data-direction` and `data-scale`; a group of an assembly
  `data-assembly`; `renderSvg`'s `viewport` in positions takes `{ system, position }` ends for a section of blocks.
- **Stories**: a step across blocks is a `view` tween in units (it passes the chromosomes between); a change of a
  block's direction, scale or the block list is a `cut` (stories, as for an axis change).
- **A hue per block** (EC-218, and RN's "a hue per block"): every item drawn in a block has the built-in field `block`
  (its block's label); a categorical theme over it takes the section's blocks in order as its domain (never first
  occurrence by region, EC-191), and `system` the same over an assembly's names, in the assembly's `order` before any
  `include` (the assembly's or a section's), so a sequence has one hue on every section of its assembly, a ring of one
  chromosome and a stub pointing at it included. A pair drawn between two places takes the field of its end on the
  first place's side: in a bridge, the bridge's `a` side (LX); on one section, the data's `a` end (RN). Why here: the
  domain is the block list or the assembly's names, which only the resolved view knows.

### What layers across blocks and sections read from here (for LX, RR and RN)

- **A ribbon's twist** (LX): a ribbon joining a run of a map (or a pairs row) between two blocks is drawn crossed when
  `run direction × display direction of block A × display direction of block B = −1` (a block's display direction is
  the section's times its own), straight otherwise. So an inversion is a crossed ribbon in each genome's own order, and
  reversing the block untwists it: many small inversions need no cut blocks.
- **From a system and a position to its places** (GA builds, GC2): each section of blocks keeps an index, by system,
  of the blocks that show it, sorted by their range's start (built once per resolved section); `placesOf(system,
  position)` returns every block and offset showing it. **Overlapping blocks**: a pairs end, a ribbon end, a projected
  hover or a selection is drawn at every block that shows its position (the merged-copies rule), each drawing counted
  against the layer's bound (`pairs.maxArcs`); a hit names the drawn run's block. A position no block shows has no
  place: the ribbon's end is a stub (in a bridge, core's: LX; on a ring, RN's outward stub) or the pair is not drawn,
  counted in the issue.
- **Hit tests are per drawn run**: a ribbon's hit is the run as drawn, its two blocks and ranges (one block, the
  blocks of one alignment merged under a pixel, or a bin pair), selected as its `(rect :a … :b …)`; one map's
  thousands of runs are hit one by one, through the layer's density cull as `pairs` today.
- **Pairs across blocks of one section** (a genome ring's 12,000 links, RN): the coordinates are GA's (the index above,
  ends on blocks of several systems: built for chords between blocks on a ring), reading them by region and by level is
  RR's (the pairs level, Gate 0, 4), and the bundled picture RN's (`rondo:ribbons`). Recommendation: no new owner; the
  built `pairs` layer stays core's, RR extends its reads, RN its drawing, each against this index.
- **Grids across blocks** (EC-223, a Hi-C map over a genome): a grid's x and y on sections of blocks read 2D tiles by
  block pair (RR's); its cells drawn per block pair's strip pair. Matrices by sample over positions (EC-220 to EC-222)
  are rows over one block's axis, nothing new from coordinates (RR).

### A navigator over several systems (core's, Gate 0, 8)

- **`navigator.system` may name an assembly** (or a list): the picture lays it out as a section of its blocks (a row,
  in the bound section's direction; the same `order`, `include`, directions and scales when the bound section is of it),
  each block's content drawn on its strip (`bands`: each chromosome's bands from one bands annotation over the
  assembly; `axis`: the blocks and their labels; `tracks`: the miniature by block). Why: a navigator already draws a
  section of its system through the layout interface; an assembly's section is blocks.
- **The window** is the bound section's viewport carried along its path, as before: the same block list, identity in
  units (a window across two chromosomes one window over the gap); a section of one chromosome, into that chromosome's
  block, the other blocks washed (unreachable, as built). So **switching chromosomes is a jump on a section of the
  assembly**: a click on another chromosome centres the view there (`navigator.center`), a sweep zooms to the region
  swept, Alt+Left returns; a section of one chromosome cannot switch (its `of` is the figure's), and its navigator says
  so by the wash.
- **The zoom strip over several systems**: on a navigator over blocks, when the window is narrower than
  `navigator.zoomStrip.minPx`, the zoom strip shows the block holding the window's centre (IGV's chromosome ideogram
  under the genome's), its real window, grips and gestures; a window across two blocks shows their hull. Why: the
  genome's and the chromosome's scales are two, as the built zoom strip's context and region are.
- **Double-click** on a navigator over blocks: fits the bound section to the block under the pointer, and from that
  fit back to the whole (`navigator.doubleClick: toggle-fit`, one level deeper); on a block's label in the view, the
  same. Why: "open this chromosome" is the precedents' first gesture.
- **Two axes over blocks (the square)**: a navigator with `axes: both` bound to a grid section whose x and y are of
  assemblies (or bound to two sections, `of: [x, y]` with `axes: both`) draws blocks × blocks, its window the x and y
  windows as one rectangle; a click centres both axes, a double-click fits both to the cell's chromosome pair, the
  grips resize each axis. This closes the navigators design's "no navigator for a grid's own y window". Its picture
  (Rondo's square) is a `navigator-picture` entry, RN's.
- **One section on both axes, and the switch**: a navigator may bind one section on both axes (`of: [s, s]`). Its
  lattice is the assembly's (after the assembly's `include`, never the section's), its cells the section shows drawn
  plain and the rest washed ("not shown"), every cell still clickable under the switch. **`navigator.switch`**
  (`none`, the default; `include`): a click on a block or a cell runs `section.include` on the bound section with that
  block or pair (a click on what it already shows: no change), a click on the navigator's title or Home on its focus
  removes it; the window, where the bound section has one, keeps its gestures. Why in core: a ring's window is idle
  unless followed, so a centring click does nothing there.
- **What the square reads, bounded**: each axis lays out by the rows' run rule, so at 3,000 contigs it is the primary
  chromosomes and a few runs, never 9 million cells. For pairs there is no 2D summary: the picture reads the pairs
  level's coarsest tiles **by block** (one want a block: a pairs tile holds every pair with an end in it, RR), sums the
  bins into one image of block × block cells and reduces it to its pixels by the shared reduction, so its cost is the
  blocks, not their square. For a grid source it reads the coarsest 2D tiles by block pair, the upper triangle when
  symmetric, at most `navigator.square.maxCells` (4,096 block or run pairs, `RenderConfig`) wants, past it the
  whole-assembly level (RR). The image is always drawn; past `maxCells` the lattice's lines are not, said once.
- **The readout and the slider text** name the block: "chr17:7,668,402" on hover; `aria-valuetext` "chr17:
  7,668,402–7,687,550 of GRCh38", across blocks "chr16:90,100,000 – chr17:1,200,000".


## Edge cases this design takes or assigns

- **Taken here**: EC-102 (gaps at scale, runs), EC-112 (a pending assembly), EC-154 (amended: the cursor crosses
  blocks), EC-184 (navigation across chromosomes, switching by `section.include` and `navigator.switch`; zoom levels
  at whole-chromosome views by the coarse reads per block), EC-186 to EC-189 (an assembly's source, name, names and
  lengths, a sequence without a length refused, aliases), EC-191 (the block-ordered domain; `system`'s in the
  assembly's order before `include`),
  EC-200 (names through forms, derived blocks keyed), EC-206 (direction on runs, a block's direction for display only),
  EC-207 (one map over sequence pairs), EC-210 (two assemblies on one section), EC-218 (the `block` field), EC-233
  (complement per block, the axis design's `auto`).
- **Its coordinate part here, the rest assigned**: EC-211 (projection through a map read by region here; no direction
  or axis derived from one; a follow through it refused), EC-212 (projection through the assemblies' map here; its
  runs GF's `map.source`, the liftover figure LX's).
- **Assigned**: EC-216 (points at a height over a genome: the Manhattan family, deferred to the close-out slice; the
  `block` field here), EC-220 to EC-222 (sample-by-position matrices, EC-221's gene-by-sample heatmap included: rows on
  one block's axis, RR's matrix), EC-223 (a contact map across blocks: RR's 2D tiles by block pair, LX's grid), EC-208
  and EC-209 (RR, LX), EC-213 (MAF: cut), EC-228 (open by gene name: GF's name index, the locus form's fallback).

## Defaults, with a picture

What a reader sees with no setting but `of` (the builder's `genome()` adds a ruler in `units: bases` and an ideogram
navigator; both spelled, never inferred).

1. **One assembly, tracks, the fit view**: blocks end to end, labels at each block's start (as built), each its own
   ruler in Mb, the gaps sharing at most a quarter of the width, 12 px each at most (EC-102); the default rail under
   the ruler, its window the whole row (idle at the fit).
   ```
   chr1             chr2            chr3          …  chr22 chrX      chrY chrM
   |0   100   200 Mb| |0   100   200 | |0  100  190|   |··| |0  100 | |··| |
   ▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬  (rail, idle)
   ```
2. **The closest zoom**: inside one block, its letters where the assembly's source serves them; the next chromosome only
   past the gap (`cross`) or never (`stop`).
3. **A ring**: arcs clockwise from the top, the gaps by the same share rule (1° each at most), labels outside, the hole
   titled with the assembly's name and total length.
   ```
                 chr1
            chrM  ╭──────╮  chr2
           chrY ╭─╯      ╰─╮
          chrX  │ GRCh38   │ chr3
                │  3.1 Gb  │
           chr22╰─╮      ╭─╯ chr4
                  ╰──────╯   …
   ```
4. **Two assemblies**: a group gap (3 block gaps; 4° on a ring, taken before the shared gaps) and the assembly's name
   on a bracket under a row or an arc outside a ring's labels.
   ```
   ┌──────── peach ────────┐   ┌──────────── grape ────────────┐
   |1  | |2 | |3  | … |8 |      |1 | |2 | |3  | … |19 |
   ```
5. **A reversed block and a scaled block**: the direction mark in the label, the ruler counting down, letters
   complemented under `auto`; a scaled block's label `×N` and its ruler at its own density.
   ```
   chr3                 chr4 ←                 chr17 ×10
   |0   100   190 Mb|   |190   100    0 Mb|    |0   2   4   6   8 Mb … 83|
   ```
6. **Pending**: before the assembly's source opens, the section's axis drawn pending (the hatched axis of EC-112)
   with the assembly's name; no blocks, no labels, no window, the rail idle.
   ```
   GRCh38 …
   ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░
   ```
7. **Crowded** (GRCh38 with alternates and decoys, about 3,400 sequences, every one included): the primary chromosomes
   as above, then one run per stretch of blocks under 2 px, labelled when it has room; on a ring the same at the
   outer radius. Past `maxAssemblyBlocks`, a closing run "+N sequences not shown".
   ```
   chr1 … chr22 chrX chrY chrM ▕▒▒▒▒ 3,366 contigs ▒▒▒▒▏
   ```
8. **A navigator over the assembly** (the builder's ideogram): every chromosome's bands, the window floored and centred
   (16 px), and under it the zoom strip of the current chromosome with the real window and connector lines.
   ```
   ▐▌▐▐▌▌ ▐▌▐▌▐ ▐▌▐▐▌ … (genome, window floored on chr17)
             ╲      ╱
   ▐▐▌▌▌▐▐▌[▐▌]▐▌▐▐▌▌▐   (chr17, the real window)
   ```
9. **Dark**: tokens only (`seq-window`, `seq-wash`, the `giemsa` grays from `msui-text`, the group bracket and a run's
   label in `msui-text-muted`, a run's band in `seq-wash`); nothing hard-coded.

## The interaction matrix's cells

Cells this design fills or changes; the others stand as the axis, layouts and navigators designs give them.

| | select | extend | add / remove | focus | back | zoom | pan | hide | edit |
|---|---|---|---|---|---|---|---|---|---|
| **a section of an assembly** | no: a section is not selected whole (select its blocks or a range; Cmd/Ctrl+A selects all, as built) | — | — | `(locus :text …)` or a block's address through `view.focus`; `block.next`, `block.previous` (Alt+PageDown, Alt+PageUp) | jumps across chromosomes are steps (Alt+Left) | about the pointer; with `stop` clamped to the block | across blocks (`cross`) or stopped (`stop`) | no: `include` by patch | `order`, `include`, `of` by patch; Shift+R reverses every block |
| **an assembly's group** (one of two on a section) | click on its label: every block of it | Shift+click on the other's label: both | Cmd/Ctrl+click | double-click on its label: fit to it, again: the whole | as built | — | — | no: `of` by patch | no: its params are the data's |
| **a block** (a chromosome) | click on its label: the whole block on its system | Shift+click on another block's label: every block between | Cmd/Ctrl+click on a label | double-click on its label: fit to it, again: the whole | as built | stacked: about the block (EC-115) | — | its label's menu: "Show only" (`section.include`, one undo step), "Show all" | `block.reverse` (its label's menu, Shift+F10), `scale` and `direction` by patch: one undo step each |
| **a run of small blocks** | click: the block under the point | Shift+click: every block between | Cmd/Ctrl+click | double-click: fit to the run | as built | — | — | no | no |
| **an element** across a block's edge | as built | Shift+arrows cross into the next block, announced | as built | as built | as built | — | — | — | — |
| **a range** across blocks | drag: shown positions of each block between, a part a system | Shift+drag, the same | Cmd/Ctrl+drag, Alt+drag, the same | Enter: its hull in display offsets | as built | — | — | — | — |
| **a bookmark or link** | — | — | — | a locus or `(sequence :name …)` address read through aliases; an alias naming two is said, never guessed | as built | — | — | — | — |
| **a navigator over several systems** | no (as navigators) | no | no | click on a chromosome: centres there (a jump); with `switch: include`, a click on a block or a cell shows only it, its title or Home shows all; double-click: the block, then the whole; PageUp, PageDown on its slider: the previous, next block | Alt+Left; undo after a switch | grips; a sweep across blocks; on two axes each axis's grips | drag the window across blocks; on two axes both | as navigators | its params by patch |
| **the view** | — | — | — | 0 fits every block (the whole genome) | history by system and position | — | — | — | undo returns a block's direction or scale |

## Configuration

| Item | Whose | Default | Why the default; the alternatives |
|---|---|---|---|
| `assembly.name`, `aliases`, `source`, `sequences`, `alias_source`, `names` | spec params | — (one of `source`, `sequences` required; every sequence with a length) | what the genome is (EC-186 to EC-189); lengths inline or from a `chrom-sizes` source |
| `assembly.order` | spec param | `natural` (numbers, X, Y, M, MT, the rest) | one genome draws in one order whatever the file; `length`, `source`, a list |
| `assembly.include` | spec param | every sequence | which contigs matter is the figure's; names, patterns, `{ min_length }` |
| `assembly.match` | spec param | `chr` | UCSC and Ensembl names differ by the prefix; `exact` |
| `section.of` | spec param | — | a system, a complex, an assembly, or a list of assemblies (new) |
| `section.include` | spec param | absent (every block of `of`) | the chromosome or pair a section shows; names, patterns, `{ min_length }` |
| `navigator.switch` | spec param | `none` | a click centres; `include` sets the bound section's `include` (the square over a ring) |
| `block.direction` | spec param | `section` | the section reads one way unless told; `reversed` |
| `block.scale` | spec param | 1 | lengths true unless the figure says; any positive number, labelled |
| `axis.directionMark` | `RenderConfig` preset (built) | `label` | now drawn per reversed block too; `none` |
| `blocks.gapShare`, `rings.gapShare` | `RenderConfig` preset | 0.25 (tracks: built; rings: new) | positions keep three quarters at any count; any share |
| `blocks.runPx` | `RenderConfig` preset | 2 CSS px | below it a gap and a label are noise; any px, 0 for never |
| `blocks.readPx` (joint with RR) | `RenderConfig` preset | 1 device px | a block's data read while it can be seen; 0 reads every block |
| `blocks.groupGap`, `rings.groupGap` | `RenderConfig` preset | 3 block gaps; 4° | a boundary between genomes reads wider than between chromosomes (Circos); any |
| `navigator.square.maxCells` | `RenderConfig` preset | 4,096 | the lattice's lines stay readable and a grid source's block-pair wants bounded; the image is always drawn |
| `navigate.acrossBlocks` | `ViewerConfig` | `cross` | as built, JBrowse; `stop` (IGV, UCSC) |
| keys of `block.next`, `block.previous` | `ViewerConfig` (`view.keys`) | Alt+PageDown, Alt+PageUp | a keyboard way to the next chromosome; any |
| `ResolveConfig.maxAssemblyBlocks` | `ResolveConfig` bound, sections of assemblies only | 4,096 | GRCh38 with alternates and decoys; past it a closing run and an issue; complexes keep `maxUnits` |
| `colinearBlocks`' `min_length` | option of a builder call | 100 kb | shorter runs merged into a neighbour, so a chain's thousands of runs give tens of blocks; any |
| the locus text's grammar (1-based inclusive, commas, longest name first) | fixed | — | the precedents' text, typed by readers as they know it |
| a block's ruler at its own density; overlapping blocks lit together; a ribbon's twist by the three directions | fixed | — | the reasons in the Decision |
| the hue per block's domain (the blocks in order) | fixed | — | known up front, never by region (EC-191) |

Not built: the navigator rows, `navigate.acrossBlocks`, the keys, `navigator.switch` (GC3); data over an assembly read by sequence (its unread info stays).

## Tensions

- **A direction per block, not per piece**: an inversion inside a collapsed (exons-only) axis on one system, reversed
  in place with breaks and no gap, is not drawn; blocks give it with a gap. The axis design's tension is narrowed to
  this, not closed. Standing.
- **Follow across unequal block lists** (B following A through a chain): refused with its reason; synteny sections move
  apart, a navigator or the selection links them. Open until a figure asks; the options stay layouts' (a follow names
  a block, or carries each block along its own map).
- **`natural` order is a rule over names**: an assembly named in words (a plant's `Pp01` … `Pp08`) sorts by the source's
  order, which may differ by format. Recommendation: accepted; `order` as a list or `length` fixes any figure.
- **`match: chr` by default** guesses less than IGV and more than an exact match: `MT` and `chrM` need an alias, and
  must not get one where they are different sequences (hg19, GRCh37). Standing; the issue names the unmatched names.
- **Reads at the whole genome**: 25 blocks × layers × a tile or two is cheap; 3,400 blocks × layers is more requests
  than pixels. Taken by the joint rule (a block read while a device pixel wide, pairs through their other end); its
  cost is that a run's values under a pixel are not drawn at the fit. Open until GC2's measure on the 3,000-contig
  fixture.
- **Units are not positions across a scaled block**: a pan's speed changes as the window crosses a scaled block.
  Accepted: the scale is the figure's, labelled.
- **The closing run past `maxAssemblyBlocks`** shows a count, not the sequences; a reader cannot reach them without a
  patch to `include`. Accepted: the bound protects resolve and layout; the issue names the remedy.

## Self-check

The claim: **this design covers what a user does with genome coordinates.**

1. **Assume it holds.** An assembly gives its chromosomes as blocks of one section, linear or ring, one assembly or
   two; each block turns over or scales by its own param; positions, locations and addresses stay the systems'; the
   navigator over the assembly moves the view across chromosomes; the zoomed ring is a cut.
2. **What it rests on, negated.** (a) A genome is seen whole, chromosomes end to end; negated: one chromosome at a
   time. (b) One assembly a section, of a few dozen sequences; plural: two, a pangenome's ten, a draft's 100,000
   scaffolds. (c) A block reads one way; negated: an inversion inside it, many small ones. (d) Names are known, unique,
   ids, spelled one way; negated: names arriving after open, aliases, an HLA name, two assemblies sharing `chr1`, a
   source naming a sequence not there, two names that look alike and are not (chrM, MT). (e) A unit is a position;
   negated: a scaled block, a contig under a pixel. (f) The user goes forward; inverse: back across a switch, a saved
   view restored after the block list changed. (g) A navigator moves one axis of one system; negated: several systems,
   two axes, thousands squared. (h) The section leads; symmetric: a chromosome section following the genome, the
   genome following it. (i) Data is held when it is drawn; negated: a map arriving by region decides a direction.
3. **The negated world.** (a) the UCSC reader needing a stop at the chromosome's end; (b) group gaps, a ring with no
   room for fixed gaps, resolve quadratic in the blocks, a bound that must say what it hides; (c) the inversion needs
   its own block, and many small ones need twisted ribbons, not blocks; (d) pending sections, alias tables, the prefix
   rule and its limit, a quoted name, a block with no id, a longest-first locus, group labels, a counted unknown name;
   (e) labels and rulers that say the scale, runs for contigs; (f) saved ends by system and position; (g) a window over
   the gap, the zoom strip of the current chromosome, the square bounded by runs; (h) a follow between unequal block
   lists has no one map to carry; (i) a block turning over as regions arrive.
4. **Both held.** Integrated: `navigate.acrossBlocks: stop` (a); `of` a list, groups, the share rule for gaps on rows
   and rings, `include: { min_length }`, `order: length`, `maxAssemblyBlocks` with its closing run, linear resolve (b);
   `colinearBlocks`, a reversed block, the ribbon's twist rule and the places index for LX (c); resolve's coverage
   input and the pending section, `alias_source`, `match`, the length check, `(block :name …)`, `(locus …)`, group
   labels, `seq.resolve.assembly.name` (d); `×N` labels, runs (e); `{ system, position }` ends and the one history rule
   (f); the navigator over an assembly, its zoom strip, the two-axis window, `navigator.square.maxCells` (g); follow
   lifted for equal block lists (h); `{ along }` dropped from blocks (i). Standing: a direction per piece inside one
   block (c); follow across unequal block lists (h); reads inside runs (e, until measured); `natural` order on worded
   names (d); the closing run's hidden sequences (b).

## Consequences

- **`seqx-core`** (GC1: language and resolve; GC2: axis and layouts):
  - language (the next minor version): the `assembly` kind and its params; `section.of` an assembly or a list;
    `section.include`; `block.direction`, `block.scale`; `navigator.switch`; `target` an assembly; the
    `(locus :text …)` and `(block :name …)` forms; the builder's `genome()`, `assembly()`, `colinearBlocks()`; issues
    `seq.resolve.assembly.alias`, `.mismatch`, `.name`, `.lengths`, `seq.view.section.blocks`, `seq.query.locus.name`.
  - resolve: derived sequences from `sequences` or the coverage (`ResolveOptions`, one content hash per source); the
    names rule (aliases, `chr`, the length check); a map between assemblies as one map over sequence pairs (EC-207);
    linear in the sequences (`blockAxes`' groups by `Map`).
  - view: `blockAxes` with a block's direction and scale (blocks on one system merge only on one direction), derived
    blocks keyed by name, overrides by `of`; groups; `maxAssemblyBlocks`; equal block lists follow;
    `seq.render.follow.blocks`' message; the places index (`placesOf`).
  - render: block strips with their own direction and scale (rulers, the direction mark, `×N`); the gap share on rings
    (`rings.ts` amended) and runs in both layouts; group gaps and labels; the `block` field for themes; the SVG's
    `data-system`, `data-direction`, `data-scale`, `data-assembly`; `renderSvg`'s viewport by `{ system, position }`.
  - navigators (GC3): `navigator.system` an assembly; windows over blocks; the zoom strip of the current block; two
    axes over blocks × blocks with its bound, one section on both axes, its lattice the assembly's; `navigator.switch`;
    the square's reads by block (pairs) and block pair (grids); readouts and slider text with block names.
- **`seqx-viewer`** (GC3): `navigate.acrossBlocks`; the cursor across blocks (EC-154 amended); `block.next`,
  `block.previous`, `block.reverse`, `section.include` (change kinds `navigate`, `navigate`, `patch`, `patch`) with keys
  in `view.keys`; saved and history viewports by `{ system, position }`; the navigator's double-click by block; the
  label's menu entries ("Reverse", "Show only", "Show all"; lazy with the menu chunk); a run's hover.
- **`seqx-examples`** (GD's): the two toy assemblies (an inversion, a translocation, a gap, aliases, an HLA contig)
  as `assembly` nodes; the whole genome in a row and on a ring; peach and grape on one ring; B in A's order with the
  inverted block reversed; a chromosome at a time with the ideogram navigator; a generated assembly of 3,000 contigs
  for the runs.
- **Storybook**: a Design page "Genome coordinates"; stories per default above, light and dark, canvas and SVG, each
  with its headless twin.
- **Tests**: property tests (a block's direction reversed twice is the identity; hit-testing the inverse of placement
  on blocks of mixed direction and scale, runs included; offsets round-trip through `{ system, position }`; a drag
  across blocks selects exactly the shown positions; switching blocks keeps the selection's addresses; the names rule
  never maps an alias to two sequences or two lengths; a ribbon's twist flips with either block's direction); drawing
  against the data at the fit view, across a block gap, inside a run and across a tile edge; a ring of 3,000 blocks
  keeps three quarters of its circle for positions; the locus grammar's conformance table (aliases, prefixes, commas,
  an HLA name, unknown names); EC-102, EC-154, EC-184, EC-186 to EC-189, EC-200, EC-206, EC-207, EC-210, EC-218 named by
  tests. Performance (each with its path): resolve of 100,000 sequences under 50 ms (quick representative, linear);
  a whole-genome layout of 25 blocks × 3 layers at the fit and of 3,400 blocks with runs (under `itWithin`); a pan
  across a block gap at 0 bytes warm; the window's move on a navigator over 25 blocks in 0.5 ms; the heavy lane's whole
  genome by region.
- **`_spec` sections that change**: core Language (assembly, section, block, locus), Resolve (assemblies, names,
  maps between assemblies), Resolved view (blocks: direction, scale, groups, keys, follow), Render (Blocks, Rings,
  Navigators, Export, Configuration); viewer Models (history by system), Input, Commands, Bindings, Configuration.
- **Other designs amended**: [axis](axis.md) (the tension "One direction per axis" narrowed to a
  direction per piece inside a block; "Not taken: a per-block direction" reversed); [layouts](layouts.md)
  (scaled segments and the zoomed ring decided, a ring's zoom a section's `include`; the gap share on rings and runs;
  "Two genomes on one ring" built; follow for equal block lists); [navigators](navigators.md) (the
  refused row "switching chromosomes" built, `navigator.switch`; the grid's y window closed; Self-check's standing
  residue closed);
  [extension-openings](extension-openings.md) (Part 2's refusal: answered here);
  [addressing](addressing.md) (the `(locus …)` and `(block :name …)` forms; derived sequences);
  [locations-and-queries](locations-and-queries.md) (the locus generator);
  [maps](maps.md) and [building-blocks](building-blocks.md) (a map between assemblies:
  one map, a step between any of their sequences); [molecules](molecules.md) (overrides by `of` for
  assemblies); [features](features.md) ("Two assemblies on one axis" planned with GC1 to GC3);
  [edge-cases](edge-cases.md) (EC-154 amended).

## Decided by the author

2026-11-21: every gap as recommended (the author on the drafts: "generally reasonable and we can fold in the reviews";
on the whole-designs review's twelve recommendations: "yes with the recommendations"; on the fold-in: "yes"). An item
the author decided by name carries its own date.

1. **A direction per block, not per piece** (an inversion drawn as its own block, with a gap; per piece inside a block
   stays open): **decided yes**.
2. **No `{ along }` on a block** (a map read by region would turn it over; `colinearBlocks` writes `reversed`): **decided yes**.
3. **An assembly's sequences enter resolve from the coverage** (as lengths do), not as dynamic nodes: **decided yes**.
4. **`order: natural` by default** (numbers, X, Y, M, MT, the rest; not the source's order): **decided yes**.
5. **`match: chr` by default** (`1` is `chr1`; never `M` to `MT`): **decided yes**.
6. **Keep the kind's name `assembly`** (not `genome`), despite mmCIF's `assembly` option: **decided yes**.
7. **Gaps at scale by a share on rows and rings** (rings amended from a fixed 1°), small blocks as runs below 2 px: **decided yes**.
8. **`maxAssemblyBlocks` (4,096), sections of assemblies only, past it a closing run and an issue** (not refused whole
   as complexes are): **decided yes**.
9. **`block.scale` built** (Circos's scaled segments, labelled `×N`), the zoomed ring a cut only (`section.include` or a second section): **decided yes**.
10. **Pan across chromosomes by default** (`navigate.acrossBlocks: cross`, `stop` the option): **decided yes**.
11. **The keyboard cursor crosses into the next block**, amending EC-154's built rule (chains in a complex too): **decided yes**.
12. **Follow lifted for equal block lists only**; B following A through a chain stays refused: **decided yes**.
13. **`section.reverse` keeps the block order** (as built; JBrowse mirrors the order too): **decided yes**, keep (chr1 stays left,
    every block turned over; a mirror of the order is `order` by patch).
14. **`(locus :text …)` and `(block :name …)` as address forms**: **decided yes**.
15. **Pairs across blocks of one section**: no new owner (core's `pairs` on GA's places index, RR's reads, RN's
    picture): **decided yes**.
16. **The `block` field and its block-ordered domain** here (EC-218), not in RN: **decided yes**.
