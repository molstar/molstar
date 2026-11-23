# Version 1: the architecture reconciled

## Intent

After the independent review (2026-11-17, an advisor reading the designs without this session's conclusions), the
author: "so this is not just about scope, but about reconciling the architecture, fixing things from the "reply",
etc", with one version, a codec-independent tree core, and the review's other recommendations accepted ("rest as
recommended").

This note is the one place that says what version 1 is: the architecture as reconciled, what it fixes from the
previous implementation and the JBrowse reply, what it defers, and the order of work.

## The architecture, in one page

- **One language, one version**, four levels: `data`, `views`, `refs`, `snapshots` ([composition.md](composition.md)).
- **A codec-independent tree core** in `@molstar/core/repr` (tree, schema, patch, builder, codec interface); JSON and
  S-expression codecs in core, MessagePack in core's `data/msgpack` with the tree's codec on it in `repr/codecs/msgpack`, io re-exporting both ([tree-core.md](tree-core.md)).
- **JSON by default, MessagePack for large inline data**; typed columns declared in the description, hashed by bytes,
  printed elided; large data in sources ([format.md](format.md)).
- **Three compiled forms in three tiers**: the static model and the resolved view (from the spec), the overlay (dynamic
  nodes from the store), the draw lists (per frame) ([composition.md](composition.md), [rendering.md](rendering.md)).
- **A data store** with keyed, versioned entries, pending vs empty, settle, open domains at settle
  ([data-store.md](data-store.md)); answers carry nodes as well as tables, the resolved model's dynamic overlay
  ([dynamic-nodes.md](dynamic-nodes.md)).
- **Maps** by kind, sourcing and encoding, structure maps included ([maps.md](maps.md)).
- **Addresses** as `path/(form :key value)`, with strings, lists and URL encoding specified; what data makes reached
  by forms ([addressing.md](addressing.md)).
- **A viewer composed from a spec**, as Mol\*'s plugins: registry entries, behaviors, bindings, one typed
  configuration, three ports (assets, text metrics, structure) ([viewer-composition.md](viewer-composition.md),
  [ports.md](ports.md)).
- **Cost follows pixels**: multi-resolution structures per layer; numbers measured on the first slice
  ([performance.md](performance.md)).
- **Tests within budgets**: under 30 s for the quick lane, 2 s a folder, few records ([testing.md](testing.md)).

## What it fixes

From the JBrowse reply ([sources/jbrowse-reply](../sources/jbrowse-reply/README.md)), in version 1:

- self-describing sources: format, location, index, assembly, aliases, fields, in the spec ([sources.md](sources.md));
- one source protocol, formats converted onto it by runtimes; bigWig as the proof;
- nothing depends on a whole-data manifest: open domains at settle, unlisted categories by a stable hash;
- strand (direction) on every map run and on pairs;
- unknown kinds kept and reported, never failing the spec; `diff` writing defaults back correctly;
- JSON readable by any tool by default (the reply's portability point; a schema stays deferred by the author).

From the previous implementation ([sources/seqstar](../sources/seqstar/README.md)): one source of truth for the
language (no tables beside it), no rxjs, no process citations in code, test tiers cut to the repository's lanes, fewer
surfaces (one UI, internal layouts and chrome parts), and its edge cases kept, changed or dropped by decision
([edge-cases.md](edge-cases.md)).

## Deferred, with their triggers

- Layers across sections, synteny, contact maps, MAF over a genome, and the **Genomes dataset** with them: a later
  milestone, after version 1's views run: the fifth slice (the author, 2026-11-19). One chromosome by region is in the
  third slice (p53's chromosome, [examples.md](examples.md)).
- The genome formats beyond bigWig (2bit, BigBed, VCF, PAF, GFF, chain files): with the Genomes milestone.
- Definitions and uses with arguments: when an example needs them (track templates over dynamic nodes are in version 1,
  [dynamic-nodes.md](dynamic-nodes.md)).
- Layouts and chrome parts as public extension points: when an outside application needs one.
- Settings layers, effective values with origin, a generated settings form: after the typed configuration has users.
- The compute port (workers): when a profile shows the main thread is the limit; its boundary decided then.
- Clipboard, storage, links as ports: when a host needs a second implementation.
- A published JSON Schema, agent tools, a React adapter, archives' writer beyond zip.

## Order of work: the first slice

The smallest end-to-end slice that tests the riskiest decisions first (the first review's, refined by the second):

1. **Tree core, minimal** ([tree-core.md](tree-core.md), [language-description.md](language-description.md)): the
   delta from meso\* (typed columns, `custom`, severity, namespaced kinds), `repr/tree`, `repr/schema` with `column`
   and `table`, `repr/builder`, the codec interface, JSON and S-expression codecs, normalization on read, the
   form-resolver table for addresses.
2. **The language, first slice** ([language-description.md](language-description.md),
   [locations-and-queries.md](locations-and-queries.md)): `SeqLanguage` with the slice's kinds, the builder, the issue
   catalog; the Location and element types; the first forms (`range`, `at`, `all`, `item`, `node`, `extent`).
3. **End to end** ([composition.md](composition.md), [protocol.md](protocol.md), [layout-and-draw-lists.md](layout-and-draw-lists.md),
   [viewer-models.md](viewer-models.md), [examples-and-stories.md](examples-and-stories.md)): inline data as `seq`
   sources through the store; the static model and the resolved view; draw lists with the `letters`, `ranges` and
   `values` layer types, the tracks layout, the passes, the token table; the SVG, recording and Canvas2D backends; the
   `SpecModel`, `ViewModel` (a scroll group, one follow), `HoverModel` with projection, `StoreModel`; the examples
   package. Two views: p53 tracks, and p53's protein on its gene (following through the codon map, the gene's full range;
   collapsed introns later).
4. **By region**, in three steps ([protocol.md](protocol.md), [dynamic-nodes.md](dynamic-nodes.md)):
   - **4a, tables by region**: a synthetic asynchronous `seq` source with a fake clock: coverage, levels, tiles, pending
     vs empty, cancellation, errors, eviction, summaries; then FASTA with its index;
   - **4b, nodes**: the source sends a transcript as a dynamic node with its item coordinate system and codon map; the
     overlay; planning through a dynamic map; node lifetime against eviction;
   - **4c, nested sources and templates**: the transcript's domains as a nested source, shown by a track template.
5. **Measure**: frame time, allocation, test time, bundle size; confirm or move the targets
   ([performance.md](performance.md)).

## Milestones after the slice

Version 1's planned features ([features.md](features.md)) come in milestones, each shippable:

1. **Core editing and addressing**: patches and `diff`, undo through the history, previews, snapshots, refs and
   bookmarks, addresses and links, selection and the grammar's defaults, the security bounds, the MessagePack codec.
2. **Figures**: the text layout, alignment figures, rings, dot plots, the remaining layers and themes, the UX
   foundations' stories and presets, SVG and PNG export, the render CLI.
3. **Integration**: structure maps and the Mol\* showcase, services, bigWig and range reads, a second consumer of the resolved view (tests only),
   archives.
4. **Genomes** (the fifth slice): the deferred list above. Its showcase: a genome ring like Rondo's, bundled and count-weighted, as an
   extension package (the author, 2026-11-19), so `ext:` layer kinds (extensibility.md) land before it.

## Open, to decide as the work reaches them

- Worker boundaries (with the compute port).
- The asset manager's range reads (with bigWig; a public API change to `@molstar/core`).
- The runtimes for structures and computed maps (with milestone 3; [maps.md](maps.md)).
- Accessibility beyond the chrome (the cursor and hover announced, forced colors, titles in exported SVG) and the
  wording of issue messages for translation: with milestone 2.
- Versioning under 0.x (every 0.x minor may break; migrations start at 1.0) and extension versions: before the first
  release.
