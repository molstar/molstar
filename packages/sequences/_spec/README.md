# Seq\*: the specification

Seq\* is SeqViewSpec, a format for sequence views built on one coordinate model: positions on anything from a single
residue to a whole genome, and the maps between them (a protein's residues onto its gene's codons across introns, an
alignment's rows onto its columns, a design onto its wild type, a chain onto its structure); data placed on any of these
shows on the others. Around it: its data model, SeqQL and the builder. The Seq\* viewer renders it, and so does
Nightingale. This spec says what the code cannot say for itself: what a document means, the concepts the packages share,
what a host can rely on, SeqQL's semantics, and the streaming protocol. The code is the truth for how. Where the code
holds a table (parameters, symbols, commands, config), the reference is generated from it into [reference/](reference/).

## What Seq\* is

- **SeqViewSpec** (SVS), a format: a document says what a sequence view shows (data defined once, views over it) and
  never how a viewer draws it. Patches change it by ref; snapshots make it a story.
- **A data model**: molecules, units, coordinate systems and the maps between them, alignments, annotations, locations.
  Platform-free and allocation-free on its hot paths.
- **The viewer**, the reference renderer: a headless context (state, managers, behaviors, commands) with a canvas
  engine, an SVG backend and two default UIs, React's and the DOM's. Nightingale is a second renderer
  (`@molstar/seq-nightingale`), reporting what it cannot draw.
- **SeqQL**, a query language over the model, written as S-expressions, stored as JSON.
- **Streaming**: data served by region and level through a host's provider, never whole.

## Packages

Two packages and two adapters. Each folder of a package is an entry (a subpath export); tree shaking keeps what a host
does not import out of its bundle. An entry names what it exports, the host's contract (what a host or an extension
uses, the registry's entry types included), and every type a public signature names is exported by some entry, so a host
can name it; its internals are not exported, and tests reach them by path.

**`@molstar/seq-core`**: platform-free; of other packages only Mol\*'s MessagePack codec (`@molstar/io`), zip,
synchronous task context, the Needleman–Wunsch DP that `@molstar/model` shares and the value palettes and colors that
MVS's runtime evaluates too (`@molstar/core`); Node's built-ins only in its CLIs and the tile server.

| Entry                           | Holds                                                                                          | Imports (of the core)                                   |
| ------------------------------- | ---------------------------------------------------------------------------------------------- | ------------------------------------------------------- |
| `/model`                        | the data model; `/model/query`, the SeqQL runtime                                              | `query`, `algorithms/search` (from `/model/query` only) |
| `/query`                        | SeqQL: the symbol table, reader, binder, printer, checker, skill                               | nothing                                                 |
| `/algorithms`                   | aligners (biomsa, Needleman–Wunsch, Myers); `/algorithms/search`, the matchers                 | nothing                                                 |
| `/translate`                    | translation between coordinate systems, a host's included; `Coordinates`, the system registry  | `model`                                                 |
| `/format`                       | the format: tables, validation, loading, resolving, patches, `.svsx` reading; the `svs-*` CLIs | `model`, `query`, `algorithms`                          |
| `/builder`                      | the builder; `/builder/svsx`, writing and bundling archives                                    | `format` (types only, but `/builder/svsx`)              |
| `/links`                        | command links for a host's markdown                                                            | `format` (types)                                        |
| `/stream/tiles`, `/stream/gmod` | stream providers: `seqstar-tiles/1` (and `/stream/tiles/server`); indexed genome files         | `format`                                                |

**`@molstar/seq-viewer`**: the context (`.`), the renderer API (`/draw`), backends (`/canvas`, `/svg`), renderers
(`/renderers/*`), UIs (`/react`, `/dom`), `/testing`, `/skin.css`, the `seqv-render` CLI; and the interactions extension
(`/interactions`, its format half, and `/interactions/viewer`, its renderers). It imports the core, rxjs and, from
`/react` only, React.

**Adapters**: `@molstar/seq-nightingale` (a resolved view as EBI's Nightingale components: the format's second viewer;
imports the core). The Mol\* adapter is Mol\*'s, in the Prism viewer (`@molstar/prism-viewer/sequence/molstar`;
[libraries.md](libraries.md#the-mol-adapter)). `@molstar/seq-storybook` holds the stories and docs pages (private; its
Mol\* demo imports the adapter and Mol\*'s plugin).

No Seq\* code imports application code, `/links` or the stream providers (a host's). Of Mol\*, the libraries import the
core's MessagePack codec and zip (from Mol\*'s io and core packages); the Storybook's Mol\* demo, a host, imports
Mol\*'s plugin and Prism's adapter. The core's folder graph is `core/_test/layers.test.ts`; the viewer's edge
`viewer/_test/layers.test.ts`; Biome (the repository's overrides on `packages/sequences`, and `layers.grit`) enforces
the rest.

## Principles

- **What, not how.** A document says what to show; navigators, skins, sizes and gestures are the viewer's. A drawing
  preference in a document is a hint a viewer may ignore.
- **Data once, views over it.** Views place data through coordinate maps; switching views reloads nothing.
- **Indices are canonical.** Positions are 0-based indices, ranges half-open; labels in a named numbering are for people
  and are converted once, at load.
- **Identity is the ref.** Never content: two sequences with equal residues are two sequences.
- **Replaced whole, or patched.** A viewer proposes patches; the document's owner commits them.
- **Themes are data, skins are the viewer's.**
- **Headless first.** All state is in one context that runs in Node; drawing reads it.
- **Mechanism in the context, policy in behaviors.** What a click does is a behavior with params.
- **Never silent.** Every issue has a stable code; fallbacks and omissions are reported; ambiguity is a status, never a
  drop; a quiet wait is a note. Missing pieces (an extension, a provider) never fail a document.
- **Domain-free core.** Biology lives in tables, presets, extensions and applications.
- **Performance is a contract.** No allocation per frame, per pointer event or per run (`seqspec:viewer#performance`).

## Glossary

| Term                 | Meaning                                                                                                  |
| -------------------- | -------------------------------------------------------------------------------------------------------- |
| molecule             | what units instantiate: a sequence or a ligand                                                           |
| sequence             | a polymer's full residue sequence, unmodeled residues included                                           |
| ligand               | a molecule indexed by its heavy atoms                                                                    |
| unit                 | one instance of a molecule over a range: a chain, a copy, a pose                                         |
| complex              | units shown together; also a composite coordinate system                                                 |
| copy group           | the units over one molecule with an equal range                                                          |
| index, label         | a 0-based position; a position's name in a numbering                                                     |
| numbering            | labels for positions (`ordinal`, `row`, `uniprot:<acc>`, an author's)                                    |
| coordinate system    | what an index counts in, named by a ref                                                                  |
| coordinate map       | a relation from one system's elements onto another's positions, possibly partial and one-to-many         |
| span                 | `to` positions per `from` element: 1, or 3 through a codon step                                          |
| extent               | disjoint half-open runs in one coordinate system                                                         |
| piece                | what a block shows of its coordinate system: a range, a direction, an origin                             |
| layout map           | the map from a block's piece to the section's shown positions                                            |
| projection           | the image of an extent through coordinate maps                                                           |
| closure              | adding the elements of other systems that contain an element (a nucleotide's residue)                    |
| alignment            | rows over columns; each row a coordinate system with a gapped map onto the columns                       |
| annotation, item     | a table of items over targets, by geometry; one feature, possibly several rows                           |
| field                | a column the kernel does not define                                                                      |
| location, loci       | a pointer to one element; many locations of one kind, stored compactly                                   |
| view, section, block | one way to show the data; one axis of it; one unit or alignment on that axis                             |
| track, layer         | a row of a section; what it shows, by type                                                               |
| channel, theme       | a visual variable (color, size, glyph); a named mapping from data to a channel                           |
| highlight channel    | a named set of lit elements, owned by its caller (`links:<link>`, `set:<section>/<n>`, `results:<name>`) |
| patch, snapshot      | ref-addressed operations over a document; a titled patch over the base                                   |
| invocation           | a request to a viewer as plain data, `{command, params}`                                                 |
| streamed source      | data served by region and level through the host's stream provider                                       |
| dynamic node         | a data node a stream patch adds, keyed by (source, key)                                                  |
| provider             | computes data from data: complex, annotation, node and field providers                                   |
| extension            | a package that owns a namespace and ships registry entries                                               |

## Files

| File                             | Says                                                                                        |
| -------------------------------- | ------------------------------------------------------------------------------------------- |
| [model.md](model.md)             | the concepts every package shares: positions, coordinate systems, maps, projection, closure |
| [format.md](format.md)           | what a SeqViewSpec document means                                                           |
| [sources.md](sources.md)         | data by reference: loaders, the streaming protocol, `seqstar-tiles/1`, `.svsx`              |
| [seqql.md](seqql.md)             | what a SeqQL query means                                                                    |
| [viewer.md](viewer.md)           | the host contract: context, state, commands, behaviors, config, extension points            |
| [interaction.md](interaction.md) | hover, highlight, selection, focus, viewports, keyboard, the default UIs                    |
| [libraries.md](libraries.md)     | the packages on the core: interactions, links, search, translate, Mol\*, Nightingale        |
| [dev-phases.md](dev-phases.md)   | what the phase ids in code comments name                                                    |

## Generated references

Written by tests as file snapshots: a change to a table fails the package's `test` until the reference is regenerated
(`vitest run -u` in that package). Never edit them by hand.

| File                                                 | From                                                                                                                                                                              | Written by                                                |
| ---------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------- |
| [reference/format.md](reference/format.md)           | the parameter tables (`core/src/format/tables/`) and the kernel's and `mvs` registries: nodes, types, layer types, themes, letter maps, color schemes, custom tables, issue codes | `core/_test/format/reference.test.ts` (`printReference`)  |
| [reference/seqql.md](reference/seqql.md)             | SeqQL's symbol table (`core/src/query/table.ts`)                                                                                                                                  | `core/_test/query/skill.test.ts` (`referenceText`)        |
| [reference/seqql-skill.md](reference/seqql-skill.md) | the same table and the fixed guidance: the agent's skill                                                                                                                          | `core/_test/query/skill.test.ts` (`skillText`)            |
| [reference/viewer.md](reference/viewer.md)           | `SeqCommands`, the default behaviors, `SeqViewerConfig`, `DefaultRegistry`                                                                                                        | `viewer/_test/reference.test.ts` (`printViewerReference`) |

## Citing the spec

Code names a section as `seqspec:<file>#<slug>`, the slug being GitHub's anchor of the heading: `seqspec:model#closure`,
`seqspec:format#patches`. Cite where code implements a stated contract: once in a file's header, or on the one function
that enforces an invariant. Never cite a generated reference (the code is the table), and never cite for how something
is built. A renamed heading means its citations are rewritten; the spec check (`scripts/spec.ts check`, run before every
commit) reports any citation that names no heading.

A change to a contract updates this spec in the same commit; where spec and code disagree, fix whichever is wrong, often
the spec.
