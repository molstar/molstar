# Precedents

Other tools and formats Seq\* learns from (2026-11-17). Each: what it is, what we take, what we avoid, which designs it
informs. Three claims were checked against their docs (pretext's README, MDN's `corner-shape`, Gosling's semantic
zoom); the rest is from general knowledge of each project, details marked *(unverified)* where uncertain. Sources point
at the designs, never the other way (§self-contained).

## Grammars and specs

**Gosling**, a JSON grammar for genomic visualization rendered through HiGlass (gosling-lang.org).
- Take: semantic zoom as data (`visibility`: an operation, a measure of width, height or zoom level, a threshold, a
  target, a transition padding); `alignment: "overlay"` to stack tracks over one data; linear and circular as one
  switch; views linked by an id; a schema that drives editor completion.
- Avoid: verbose overlays; every view tied to one runtime; interaction limited to brush and zoom.
- Informs: [composition](../../design/composition.md), [performance](../../design/performance.md) (density steps as
  data), [ux-foundations](../../design/ux-foundations.md).

**Vega and Vega-Lite**, declarative visualization grammars; Vega-Lite compiles to Vega's dataflow.
- Take: scales as named objects views share, and `resolve` for shared or independent scales (our domain groups);
  selections as data predicates; event streams as an interaction grammar; a high-level spec compiled to a low-level one.
- Avoid: signal expressions as a string language hard to validate; opaque dataflow to debug; a scenegraph that slows
  past about 100k marks.
- Informs: [composition](../../design/composition.md), [ux-foundations](../../design/ux-foundations.md) (the
  grammar), [addressing](../../design/addressing.md) (selections).

**Observable Plot**, a concise layered-mark plotting library.
- Take: marks as functions over data; transforms (bin, group, stack, window) composed as option rewriters; good default
  scales; facets; a pointer transform for hover.
- Avoid: specs are code, not data; SVG only.
- Informs: [composition](../../design/composition.md), [examples](../../design/examples.md).

**MolViewSpec (Mol\*)**, a declarative tree for molecular scenes with builders and a Mol\* runtime: our closest sibling.
- Take: typed nodes with param schemas; builders apart from the runtime; annotations from files referenced by field;
  selectors over `label_*` and `auth_*`; snapshots and stories; a versioned, validated format.
- Avoid: selectors written twice (inline and in annotation files); builder and runtime versions drifting.
- Informs: [format](../../design/format.md), [addressing](../../design/addressing.md), [maps](../../design/maps.md).

**LocusZoom.js**, regional association plots.
- Take: data layers declared apart from data adapters, with namespaced fields joined across sources; layouts as JSON
  with overridable defaults; adapter caches by region; a plugin registry for adapters, layers and widgets.
- Avoid: surprising deep merges of layouts; an ad hoc field-namespace syntax *(its rework: unverified)*.
- Informs: [sources](../../design/sources.md), [extensibility](../../design/extensibility.md); the LD and Manhattan views
  ([examples](../../design/examples.md)).

## Genome browsers and data formats

**HiGlass and Clodius**, a multiscale tile-based viewer and its tiler.
- Take: tile pyramids (level z: tiles of width max / 2^z, ids `uuid.z.x`), so a frame costs what is on screen; an
  aggregation chosen per tileset; a tileset-info answer (bounds, resolutions); zoom and location locks across views;
  multivec tiles for stacked rows.
- Avoid: a preprocessing step or a server before anything shows; a view config mixing layout, data and view state; a
  track plugin API through globals.
- Informs: [performance](../../design/performance.md), [sources](../../design/sources.md) (the protocol's levels).

**JBrowse 2**, a modular web and desktop genome browser (and the source of the reply).
- Take: the adapter boundary (`getFeatures(region)` as a stream), adapters swappable by configuration; tracks as an
  adapter plus displays plus a renderer; a plugin manager registering element types; serializable sessions; renderers
  in workers; comparative views.
- Avoid: heavy configuration and session layering; a steep plugin learning curve; a large bundle; serialization cost
  across workers.
- Informs: [ports](../../design/ports.md), [extensibility](../../design/extensibility.md),
  [sources](../../design/sources.md).

**IGV and igv.js**, the desktop and embeddable genome browser.
- Take: small track configurations (`{type, format, url, indexURL}`); index-driven range reads over HTTP; caching by
  loaded interval with padding; expanded, squished and collapsed modes; locus strings in URLs.
- Avoid: display-mode logic repeated per track type; vertical stacking only; little declarative interaction.
- Informs: [sources](../../design/sources.md), [addressing](../../design/addressing.md).

**UCSC track hubs**, remotely hosted track collections (`hub.txt`, `genomes.txt`, `trackDb.txt`).
- Take: decentralized, URL-addressed distribution; composite and container tracks; display settings inherited from a
  parent; hubs for non-standard assemblies.
- Avoid: a loosely validated stanza format with hundreds of settings; renderer-specific semantics.
- Informs: [sources](../../design/sources.md), definitions and uses in [composition](../../design/composition.md).

**bigWig, bigBed, tabix and BGZF**, indexed binary formats for range queries.
- Take: bigWig's built-in zoom-level summaries (min, max, mean, sum of squares a bin) and R-tree index, the template
  for "a region at a resolution"; BGZF with `.tbi`/`.csi` for remote ranges in few requests; the GMOD readers as
  references.
- Avoid: bigBed's untyped extra fields; C-centric writing tools.
- Informs: [performance](../../design/performance.md), [sources](../../design/sources.md) (bigWig, the proof format).

**SnapGene and Benchling**, plasmid and sequence editors (SeqViz an open counterpart *(unverified as the best one)*).
- Take: circular map, linear map and sequence views in sync over one model; directed features with arrows; translation
  under the DNA in frame; restriction sites and primers as derived overlays; edits shifting downstream features; one
  selection across views; label layout on circular maps.
- Avoid: proprietary formats; edit semantics hidden in the UI.
- Informs: [ux-foundations](../../design/ux-foundations.md), the plasmid views ([examples](../../design/examples.md)).

## Protein and alignment viewers

**Nightingale (EBI), after ProtVista**, web components for protein feature tracks (UniProt, InterPro).
- Take: one custom element a track; a manager syncing the display range and highlight; a shared residue axis;
  framework-agnostic delivery.
- Avoid: synchronizing through DOM attributes (slow, stringly typed); scaling reimplemented in each component; spec and
  data inseparable.
- Informs: [rendering](../../design/rendering.md) (the reference adapter), [addressing](../../design/addressing.md).

**Jalview**, the desktop alignment editor.
- Take: the color scheme catalog as rules (Clustal, BLOSUM62, hydrophobicity, conservation thresholds); annotation
  rows under the alignment (conservation, quality, consensus, secondary structure); groups; structure links through
  residue maps with highlight both ways; hidden rows and columns that keep coordinates.
- Avoid: state scattered across panels; global mutable selection; quirky file grammars.
- Informs: [maps](../../design/maps.md), [ux-foundations](../../design/ux-foundations.md), the globins views.

**Web alignment viewers** (react-msa-viewer, AlignmentViewer 2.0, Nextclade).
- Take: tiled canvas rendering with cached glyph tiles; an overview of the whole alignment beside the detail *(Alignment
  Viewer's details: unverified)*; differences from a reference; codon-aware genome → gene → peptide coordinates
  (Nextclade); virtualized rows and columns.
- Avoid: redraws driven by a framework's re-renders; each viewer's own position mapping code.
- Informs: [performance](../../design/performance.md), [maps](../../design/maps.md).

**Skylign**, profile HMM logos.
- Take: letter heights by information content; insert probability, length and occupancy as their own rows; logos at
  fixed column width, scrollable.
- Avoid: a logo as a service apart from the view; ambiguity between entropy and probability scaling: the spec names it.
- Informs: the logo layer ([features](../../design/features.md)).

## Rendering, layout and text

**deck.gl**, a GPU layer framework.
- Take: layers as immutable descriptors diffed by id, with explicit triggers for what to recompute; data in binary
  buffers; composite layers from primitive ones; picking that returns `{layer, index, object}`; a tile layer by
  `(x, y, z)`.
- Avoid: stale renders when a trigger is forgotten: prefer versions from input hashes (our IR keys).
- Informs: [composition](../../design/composition.md) (keys and versions), [addressing](../../design/addressing.md).

**Pretext**, multiline text measurement and layout without DOM reflow (checked: its README).
- Take: `prepare` (segment, measure with a canvas, cache) then `layout` as arithmetic; line-at-a-time layout for text
  flowing around obstacles; output for canvas, SVG and DOM alike. Already adapted in `@molstar/core/util/text-measure`.
- Avoid: it throws without a canvas (Node, jsdom); px fonts only; no bidi; widths to clear after web fonts load.
- Informs: [rendering](../../design/rendering.md) (the text metrics port), [base-library](../../design/base-library.md).

**Apache Arrow**, a columnar typed-array format.
- Take: typed columns with validity; dictionary-coded categories; record batches as region chunks; zero-copy to workers.
- Avoid: the JS library's size and weight: adopt the conventions, not the dependency.
- Informs: [sources](../../design/sources.md) (the protocol's tables), [format](../../design/format.md).

**CSS `corner-shape`** (checked: MDN): `squircle` is `superellipse(2)`; limited availability, needs a
`border-radius`, and does nothing in exported SVG.
- Take: chrome as progressive enhancement; the same superellipse as a path in draw lists for data marks and exports.
- Informs: [ux-foundations](../../design/ux-foundations.md).

## Mappings and nomenclature

**SIFTS and PDBe-KB**, residue mappings between UniProt, PDB (`label_seq_id`, `auth_seq_id` with insertion codes) and
domain databases.
- Take: mappings as segments, not residue tables; gaps, insertion codes, many chains a UniProt and isoforms explicit;
  the release recorded.
- Avoid: assuming author numbering is monotonic; ignoring microheterogeneity.
- Informs: [maps](../../design/maps.md) (structure maps).

**HGVS nomenclature, GFF3 and the Sequence Ontology**.
- Take: coordinate prefixes (`g.`, `c.`, `p.`, `n.`) and intronic offsets (`c.88+2`) as addresses over maps; gene →
  transcript → CDS hierarchies typed by the ontology; one stated convention a coordinate system.
- Avoid: the off-by-one tangle across BED, GFF and VCF; untyped attributes.
- Informs: [addressing](../../design/addressing.md) (an HGVS-style form a candidate), [maps](../../design/maps.md).

**Ensembl and NCBI mapping services**.
- Take: mapping results as lists of ranges with gaps where exons split codons; liftover as one more map.
- Avoid: a network call a mapping: fetch the model once, map locally.
- Informs: [maps](../../design/maps.md) (computed maps).

## Lessons across them

1. Multi-resolution by pixels per unit, caches keyed by (source, level, tile): HiGlass, bigWig, deck.gl.
2. A narrow adapter boundary returning typed columns for a region, adapters registered by type: JBrowse, LocusZoom,
   Arrow.
3. Interaction as data, without a string expression language: Vega's selections without its expressions.
4. Shared scales and semantic zoom as spec, not code: Vega-Lite's `resolve`, Gosling's `visibility`.
5. Maps as invertible, composable segments, addressed in each system's own terms: SIFTS, HGVS, Nextclade.
6. Everything addressable by stable ids, never by DOM nodes: deck.gl picking, MolViewSpec selectors.
7. Avoid: sync through DOM attributes, heavy session trees, loosely validated settings, implicit off-by-one rules.
