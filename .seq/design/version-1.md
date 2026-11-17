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
  S-expression codecs in core, MessagePack in io ([tree-core.md](tree-core.md)).
- **JSON by default, MessagePack for large inline data**; typed columns declared in the description, hashed by bytes,
  printed elided; large data in sources ([format.md](format.md)).
- **Three compiled forms**: resolved model, resolved view (public, no data, no target), draw lists (for a size, a
  viewport, a target and the data store) ([composition.md](composition.md), [rendering.md](rendering.md)).
- **A data store** with keyed, versioned entries, pending vs empty, settle, open domains at settle
  ([data-store.md](data-store.md)).
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
  milestone, after version 1's views run.
- The genome formats beyond bigWig (2bit, BigBed, VCF, PAF, GFF, chain files): with the Genomes milestone.
- Definitions and uses with arguments, templates per member of a set: when an example needs them.
- Layouts and chrome parts as public extension points: when an outside application needs one.
- Settings layers, effective values with origin, a generated settings form: after the typed configuration has users.
- The compute port (workers): when a profile shows the main thread is the limit; its boundary decided then.
- Clipboard, storage, links as ports: when a host needs a second implementation.
- A published JSON Schema, agent tools, a React adapter, archives' writer beyond zip.

## Order of work: the first slice

The smallest end-to-end slice that tests the riskiest decisions first (the review's, accepted):

1. **Tree core, minimal**: `repr/tree` (node, values with typed columns, hashes, addresses, issues), `repr/schema`,
   `repr/builder`, the codec interface, the JSON and S-expression codecs, into `@molstar/core`.
2. **The language description, first slice**, with the query language's first forms (it is designed here, with
   addressing, not apart): `spec`, `data` with a sequence, `ranges` and `values` annotations, an `offset` map and a
   `codon` map; `views` with sections, tracks, the `letters`, `ranges` and `values` layers.
3. **End to end**: inline data through the data store → resolved model → resolved view → draw lists → the SVG and
   recording backends, and a thin Canvas2D viewer with projected hover. Two views: p53 tracks, and p53's protein on
   its gene.
4. **One source by region** through the data store: FASTA with its index (or a synthetic asynchronous source first).
5. **Measure**: frame time, allocation, test time, bundle size; confirm or move the targets.

Then: patches and `diff`, the MessagePack codec, structure maps and the Mol\* showcase, addressing and links, the
remaining layers and layouts, the typed configuration, the other datasets.

## Open, to decide as the work reaches them

- Session state's shape and whether it is saved (with the viewer's models).
- Worker boundaries (with the compute port).
- The asset manager's range reads (with bigWig; a public API change to `@molstar/core`).
