# Seq data ports: a host's data and Seq\*'s own reading, one path

A fourth-slice draft (lane SP, design only). Links are written as they will read once the file moves to `design/`.

## Intent

The author, 2026-11-19 (molecules gap 11, superseded): "we should also be able to share data representation with mol*
... maybe we could have generic representation ports, similar to mol io, we could write seq io"; "the same principle
not the same functionality". Then, on this design: "i suppose we can provide default ports which 'read data
independently' so we have the same mechanism?".

Seq\*'s nodes (sequences, molecules, complexes, units, their numbering and observed ranges) and tables (per-residue
values, ranges) can come from a host that already holds them in memory (a structure viewer's model), without a second
parse; and Seq\*'s annotations go back to the host in the host's own terms. Seq\*'s own reading (its format runtimes
over the assets port, [files.md](files.md), [sources.md](sources.md)) is the **default adapter of the same path**: a
host replaces an adapter, nothing else changes. Gate 0 (9): the design and a fake host's adapter in Storybook; no Mol\*
adapter.

## Precedents

- **The source protocol as built** ([sources.md](sources.md), [protocol.md](protocol.md), `data/protocol.ts`): a
  `FormatRuntime` opens a `SourceEntry` into an `OpenedSource` (coverage with nodes up front, `read(region)` async and
  cancellable, `issues`, typed failures); the store keys entries by address and holds nodes by (source, key). Taken
  whole: it already is the port's surface (Decision). Not taken: nothing, two additions only.
- **Mol\*'s own layering**: `@molstar/io` reads files into frames; `@molstar/model` builds a `Model` from them, PDB and
  GRO converted onto mmCIF's tables first. Taken: convert at the edge onto one shape named in mmCIF's words. Not taken:
  a CIF frame as that shape (a host holding a `Model` would have to keep or rebuild its frame).
- **JBrowse's adapters**: one interface for every format, and `FromConfigAdapter` serving features held in memory
  through it. Taken: in-memory data through the same adapter as a file. Not taken: adapters named in the document by
  type (a spec names the data's format, never a host).
- **igv.js**: a track's `features` array or a `reader` object in place of a URL, through one feature source. Taken: the
  host's reader stands where the file reader stood.
- **The previous Seq\***: opaque `stream` ids served by the host (EC-185: stripped of its host block, every example drew
  nothing). Refused: a spec opens with no host code; a host serves what the spec names in the data's own terms.
- **MolViewSpec**: data by URI only, no in-memory references. Taken: the portable name stays the URI.
- **The assets port** ([ports.md](ports.md)): a host passes its own `AssetManager`, so downloads are shared. Taken:
  the pattern (a host's adapter, the default otherwise). Not enough: shared bytes still parse twice.
- **This codebase's structure syncing** (fourth slice): `structureFacts(model, complex)` publishes the unit facts
  outward; the fake viewer holds its own recorded chain table, independent of Seq\*'s nodes. Taken: the inward
  description below is the same facts in the other direction, and the fake's chain table is its first host.

### Interactions (what a user does with shared data in the precedents)

| Interaction | Precedent | Here |
|---|---|---|
| open a structure in the 3D viewer, its chains appear as sequences | Mol\*'s sequence panel, ChimeraX | built: the host's adapter serves the complex; no second download or parse |
| the same spec without the host | MVS, JBrowse | built: the default adapter reads the file by `uri` (MC's runtime) |
| B-factors or pLDDT as a track | Mol\*'s sequence colors, AlphaFold DB's page | built: `values` on unit systems from the same source (`field`) |
| the host edits its structure (a chain removed, a residue mutated) and the sequences follow | Mol\* (the panel rebuilds) | built: the source reopened in place, node identity kept by key |
| a Seq\* track's rows as the host's annotation (domains as a 3D property) | ProtVista to LiteMol, Jalview's "color by" | built (core function); the host writes its own representation |
| load a file the host already read, by URL | Mol\* shares assets | built: the host's adapter matches the `uri`; else it delegates |

## Decision

### The port is the source protocol; adapters are registry entries by format

- **One path.** A data adapter is a `FormatRuntime` (`name`, `formats`, `open(entry, ctx) → OpenedSource`), a registry
  entry keyed by the formats it reads ([ports.md](ports.md)'s rule: several, looked up by name, so registry entries,
  not ports; "seq data ports" is the family's name). The core's runtimes (`seq`, FASTA, the alignments, the records;
  `runtimes/mmcif` and `runtimes/bigwig` opt-in) are the **default adapters**. A host lists its own adapter for a format
  instead of the default. Why: reading by region, streaming, waiting, issues and node identity are already the
  protocol's; a second mechanism would duplicate each.
- **A host adapter wraps the default** (`hostAdapter(fallback, serve)`, pattern, not a registry feature): for a source
  it holds it answers from memory; one it is loading it awaits (cancellable by `signal`); one it does not hold it hands
  to `fallback` (the default adapter, reading through the assets port). Why: a host holds some entries, not a format;
  the choice is explicit in the adapter, never an order of registration.
- **One adapter a format, or an error.** Registration refuses two entries claiming one format under different names,
  listing them (today `runtimesOf` checks names only, and `bindSource` and the overlay take the first that reads the
  format: an order). A host replaces the default through `withAdapter(registry, adapter)` (`seqx-core/context`), which
  removes the entries for the adapter's formats and adds it. Why: extensibility's "never resolved by order".
- **The spec names the data in its own terms, never a host** (EC-185): `format` and `uri`, as today. A host matches a
  source by its `uri` (normalized as its asset manager keys it); a source may add `ref` (the host's name for data with
  no file: an edited or predicted model), read by host adapters only (Gaps 2). Why: the same spec opens with the host
  and without it.
- **Inline data is never served by a host**: it reads through `seq`, as built. Why: it is the spec's own.
- **What crosses is the protocol's**: coverage (systems, shapes, levels, nodes up front), tables, typed failures,
  `issues` at the source's address, `AbortError` on cancel. A host's nodes are checked by the language as every dynamic
  node is ([dynamic-nodes.md](dynamic-nodes.md)); a host's bug is an issue, never a crash. Why: one validation, and a
  host is code Seq\* does not own.
- **Asynchronous always**: a host in memory answers with a resolved promise; there is no synchronous path. Why: one
  store path; hover and projection never read sources, so a microtask costs nothing visible.

### Seq io: one converter from a host-neutral description

- **`StructureDescription`** (`@molstar/seqx-core/io/structure`, pure, `@molstar/core` only): the columns every adapter
  of a structure fills, named in mmCIF's words, held as columns (SoA, typed arrays):
  - entities: `id`, `type` (`polymer`, `non-polymer`, `branched`, `water`), the polymer's type, `name`, canonical
    letters, components a position (with `hetero` alternatives), identifiers (CCD code, WURCS, LINUCS, PRD);
  - chain instances of the assembly asked: `label_asym_id`, `auth_asym_id`, entity, `operator`;
  - residues a chain instance: position in its entity, `auth_seq_id`, `pdbx_PDB_ins_code`, component, `observed`;
  - links (`struct_conn`-like, branch links), compounds (code, name, type, parents), heavy atoms a molecule (optional);
  - named per-residue values (`b_factor`, `occupancy`, `plddt`; any name), optional.
- **`structureNodes(description, options)`** turns it into the dynamic nodes [molecules.md](molecules.md) specifies
  (sequences, molecules, `compounds`, a `complex` of units with `auth` numbering, `observed`, `links`), with its issues.
  **Keys come from the data's identifiers, never indices** (an entity by its id, a unit by `label_asym_id` and
  `operator`), so the default adapter and every host give equal keys: a selection, a bookmark or a ref made with the file
  holds with the host, and the reverse.
- **MC's converter is shared** (the question Gate 0 left to SP): MC's runtime reads CIF categories through
  `@molstar/io` into a `StructureDescription`, then `structureNodes`; the fake host reads its chain table into one; a
  Mol\* adapter later reads a `Model` into one. Why at the description and not at CIF categories (molecules' option
  (b)): a host holding a `Model` need not keep a CIF frame; the description is format-neutral (PDB, BinaryCIF, a host's
  own); one conformance suite runs on descriptions, without CIF. Why not each adapter writing nodes: keys and the
  dictionary's reading would drift between hosts.
- **The inverse of `structureFacts`**: facts go out on units, the description comes in on entities and chains; both
  name the same identifiers (`chain`, `auth_chain`, `operator`, `auth` numbering), so a host's two adapters (data and
  structure port) translate with one vocabulary.

### What a source serves for a structure

- **Up front**: the nodes, in its coverage (held whole, small: one complex node with its units as children; bounded by
  `maxNodeBytes` and `ResolveConfig.maxUnits`).
- **By region**: `values` and `ranges` on the unit systems it brought (B-factors, pLDDT, a host's secondary
  structure), named in coverage by the units' keys. Why by region: the unit's positions are what a layer zooms over, and
  a host's computed property can arrive later than the units.
- **A values channel by `field`**: a `values` node reading a source names its column (`field: 'plddt'`), passed as
  `RegionRequest.field` (absent: the source's default, a records source's `fields.value` as today). Why: one source,
  one open and one parse serve every channel, and its units are named in the same source (no cross-source system names).

### A host's updates: the source reopened in place, identity kept

- **`OpenedSource.changes?(fn): unsubscribe`**, optional: the adapter says its data changed. The viewer's source sync
  reopens the source in place, as a transform's edit does ([files.md](files.md)): tiles dropped, opened again, attached.
- **The store's `reattach(address, opened)`**: nodes of the same key keep their identity and coordinate system id and
  take the new content (version moved); keys gone are released; new keys held. Why: EC-176's "a conflict keeps the
  held node" is right for two answers of one reading, wrong for a new reading; a host's edit is a new reading.
- **What referred to what went** (a selection part on a removed unit, a bookmark) is treated as when a patch removes
  that node, as built. A default adapter never fires `changes` (a file is read once a session, [files.md](files.md)).
- **Coalesced**: changes inside one frame reopen once (`DomHost`'s frame, as the store's other work). Why: a host edits
  in bursts.

### Annotations given back as the host's

- **Outward is a pure function, not an adapter**: `unitRows(model, store, annotation, units) → { complete, pending,
  rows }` (`seqx-core`, `structure/`, beside `colorRuns`): the annotation's rows projected onto bound units through the
  binding's map choice, unobserved cut, each row a `StructureTarget` with its fields. The host writes them as its own
  representation (a Mol\* custom property, an MVS annotation, the fake viewer's residue property). Why: the vocabulary is
  the host's (structure syncing, decided), and colors already go out this way.
- **No writer in Seq\*** for a host representation, as for MolViewSpec (decided). A host's writer may be a format
  writer registry entry of its own package.

### Layering

- **No package of its own.** `io/structure` (the description and `structureNodes`) is an entry of `seqx-core` with no
  dependency beyond `@molstar/core`; `runtimes/mmcif` stays the opt-in entry that imports `@molstar/io` (molecules gap
  12, decided); neither imports `@molstar/model`. Why: the entries already split the cost (`sideEffects: false`), and
  a package would add an inventory row, exports and a boundary for no second consumer.
- **A Mol\* adapter never lives in `seqx-core`**: it imports `@molstar/model` (outside the allowlist,
  [packages.md](packages.md)). It lives with its host: the new mol viewer, or a bridge package
  (`@molstar/seqx-molstar`: `seqx-core` + `@molstar/model`) when one is built. Not in this slice (Gate 0, 9).
- **Direction**: host → its adapter → `seqx-core/io/structure` → `@molstar/core`; Seq\* imports no host.

### The fake host's adapter (Storybook)

- Beside the fake structure viewer (`packages/seq/storybook`, `_stories/support/`): `fakeStructureAdapter(viewer,
  fallback)` reads the viewer's recorded chain table into a `StructureDescription` for the sources whose `uri` names an
  entry it holds (2OCJ, 1FIN), delegates the rest to MC's runtime, serves `values` by `field` from the table's
  per-residue columns, and fires `changes` when the viewer's "drop chain" or "mutate" control edits its table.
- **The proof of one path**: the same spec, composed with the fake adapter and with the default, gives equal nodes
  (keys and content) for 2OCJ and 1FIN, and the fake's run reads no asset.

## The interaction matrix's cells

A data path adds few gestures; rows it changes:

| Row | select | extend | add/remove | focus | back | zoom | pan | hide | edit |
|---|---|---|---|---|---|---|---|---|---|
| **an element** on a host-served unit | as the grammar | as the grammar | as the grammar | as the grammar | as the grammar | as the grammar | as the grammar | no: data | no: the host's (its edit fires `changes`) |
| **a track or layer** over a host's values | as the grammar | — | — | — | — | — | — | `track.hide` | `field` by patch (another channel, one undo step) |
| **the selection** across a host's update | kept by key | — | parts on removed units dropped, said once | — | Alt+Left: history as built | — | — | — | — |
| **a section of a complex** after an update | as [molecules.md](molecules.md) | — | — | — | — | — | — | — | reflows; no history entry (data, not an edit of the spec) |

Every other cell: no, because a data adapter carries no gesture.

## Configuration

| Item | Kind | Default | Why | Alternatives |
|---|---|---|---|---|
| a format's adapter | registry entry (`withAdapter`) | the core's runtime for it (mmCIF, bigWig opt-in) | reads the file on its own | a host's adapter wrapping it |
| a host adapter's `fallback` | adapter param | the default adapter | a spec naming data the host lacks still opens | `none`: unavailable with an issue (a host that must not fetch) |
| `field` | `values` node param | the source's default value | one source, many channels | — |
| `ref` | `source` param (Gaps 2) | absent | portable by `uri` | a host's name for data with no file |
| reopen on `changes` | viewer source sync | on, coalesced a frame | the sequences follow the host | — (an adapter that should not reopen fires nothing) |
| `StructureDescription` options | `structureNodes` options | as MC's (`assembly` absent: the deposited model; `water` off; `name: description`) | one reading for file and host | per adapter |

## Tensions

- **A portable name against data with no file.** `uri` keeps specs portable; a host's edited model has none. `ref`
  alone fails EC-185's stripped test. Recommended: `ref` beside a `uri` where one exists; `ref` alone allowed, the
  source unavailable without the host, with `seq.source.host_only` (info) so the stripped test names it (Gaps 2). Open.
- **Reopen against a diff.** Reopening drops every tile of the source on any change; a host mutating one residue
  re-sends the whole complex. Fine at structure sizes (a complex is held whole, its tiles one each); a genome-sized host
  would want region invalidation (`changes` carrying regions). Standing until such a host.
- **Wrapping against a registry chain.** A host adapter that delegates is code a host writes (`hostAdapter` helps); a
  registry with "first that accepts" would be shorter and an order. Recommended: wrapping. Standing.
- **One description against everything a host knows.** Mol\*'s `Model` holds secondary structure, coordinates,
  trajectories; the description carries what nodes need plus named values. More comes as named values or ranges by
  region, never coordinates. Standing.
- **Code names.** The protocol keeps `FormatRuntime`; the designs call it the default adapter of a data port. A rename
  (`DataAdapter`, kind `data-adapter`) would say it, at the cost of churn in every runtime and the registry (Gaps 1).

## Self-check

The claim: "one path serves Seq\*'s own reading and a host's data; a host replaces an adapter and nothing else
changes".

1. **Assume it holds.** Every node and table reaches the store through `open` and `read`; the default adapter reads
   files; a host's adapter reads its memory; keys are the data's.
2. **What it rests on, negated.** (a) *The data has a portable name.* Negated: an edited model, no file. (b) *Data is
   read once.* Negated: the host edits it. (c) *An adapter owns a whole format.* Negated: the host holds two entries of
   ten. (d) *A source serves one value.* Negated: B-factor, occupancy and pLDDT of one structure. (e) *Data flows in.*
   Inverse: Seq\*'s annotations out to the host. (f) *Keys are stable.* Negated: a host indexing by array position, so
   a removed chain renames the rest. (g) *The host's data is well formed.* Negated: a host bug. (h) *One host.* Plural:
   two hosts serving `mmcif`. (i) *A host answers at once.* Negated: still loading.
3. **The negated world.** (a) `ref`, with the stripped test naming it; (b) `changes`, `reattach`, identity by key;
   (c) the wrapping adapter and its `fallback`; (d) `field` on `values`; (e) `unitRows`, a pure function, the host
   writing; (f) keys fixed by `structureNodes` from identifiers, no adapter choosing them; (g) dynamic nodes validated as
   any; (h) one adapter a format, conflicts an error, the hosts composed in one adapter; (i) `open` awaits, cancellable.
4. **Held both.** Integrated: (b) to (i) in the Decision. Not one mechanism, by design: (e) outward is a function, not
   an adapter. Standing: (a) a host-only `ref` against portability (Gaps 2); reopen without region invalidation.

## Consequences

- **`seqx-core`**: `data/protocol.ts` (`OpenedSource.changes?`, `RegionRequest.field?`; both optional, every runtime
  unchanged); `data/store.ts` (`reattach`, node content replaced by key); `context/` (`withAdapter`, one adapter a
  format, conflicts listed; the resolver and the overlay bind by format with no order); `lang/` (`values.field`, `source.ref`
  if Gaps 2: the 0.5.0 bump through ML1); `io/structure` (new entry: `StructureDescription`, `structureNodes`);
  `structure/unitRows`. Tests: one conformance suite every adapter passes (the protocol's, plus: equal nodes from a
  description by MC's path and by a hand-built one; keys unchanged when a chain is removed; `reattach` keeps identity;
  `changes` coalesced on a fake clock; a host's malformed node an issue); properties (`structureNodes` keys independent
  of row order; a description round-tripped through nodes and `structureFacts` gives the same unit facts).
- **`seqx-viewer`**: the store model subscribes to `changes`, reopens in place a frame later. Tests: an update keeps the
  selection by key and drops parts on a removed unit with one issue.
- **MC** (brief changes): CIF categories → `StructureDescription` → `structureNodes`; its conformance suite split (the
  CIF reading against recorded counts; `structureNodes` on descriptions); `values` by `field` (B-factor, pLDDT).
- **Storybook**: the fake adapter beside the fake viewer; a `Host data` story (the fake serving 2OCJ, a pLDDT track from
  it, "drop chain" followed by the sequences) with its headless twin on a fake clock, which asserts no asset read and
  equal nodes against the default adapter; the Extensibility page lists the data adapters and `withAdapter`.
- **Designs when agreed**: [sources.md](sources.md) ("a runtime per format" read as the default adapter; a host's
  adapter replaces it); [protocol.md](protocol.md) (the two optional members); [data-store.md](data-store.md) and
  [dynamic-nodes.md](dynamic-nodes.md) (a reopen replaces content by key; EC-176's conflict rule for one reading only);
  [ports.md](ports.md) (registry entries: data adapters by format; the "seq data ports" family named);
  [molecules.md](molecules.md) (Reading: through `StructureDescription`; the Tension on the mmCIF runtime closed);
  [extensibility.md](extensibility.md) (formats row: an adapter per format, a host's included).
- **`_spec` sections that change**: core's Data (Protocol, Store), Context (Registry entries), Language (`values`,
  `source`), a new Seq io section; viewer's StoreModel.

## Gaps for the author

1. **Names.** Keep `FormatRuntime` and the `format-runtime` kind in code, the designs calling them the default adapters
   of the data ports, rather than renaming to `DataAdapter`? Recommended **yes, keep** (no churn; "ideally nothing but
   a name" met in the designs).
2. **Data with no file.** A `ref` param on `source`, read only by host adapters, with `uri` beside it where one exists
   and `ref` alone allowed but reported (`seq.source.host_only`, info)? Recommended **yes**.
3. **Host updates.** `OpenedSource.changes` and a store `reattach` that keeps node identity by key and replaces content
   (EC-176's rule kept for two answers of one reading)? Recommended **yes**.
4. **Several values in one source.** `values.field` and `RegionRequest.field`, rather than one source node per channel?
   Recommended **yes**.
5. **One adapter a format.** Two runtimes claiming one format become a registration error (today the first wins), with
   `withAdapter` to replace a default? Recommended **yes**.
6. **No package.** `io/structure` an entry of `seqx-core`, a Mol\* adapter in a package of its own only when built,
   rather than a `seqx-io` package now? Recommended **yes**.

## Decided by the author (2026-11-20)

- Gaps 1, 5, 6: "Accept all three": the code keeps `FormatRuntime` (the designs say the default adapter); one adapter
  per format, `withAdapter` replaces a default; no `seqx-io` package, the structure description in `seqx-core`'s
  `io/structure`.
- Gap 2: `source.ref` beside a `uri` or alone, alone flagged as not portable (info issue, EC-185).
- Gap 3: host updates by `changes` and `reattach`, every tile dropped and read again; region-level invalidation a
  later refinement (a tension).
- Gap 4: `values.field`, in the 0.5.0 bump.
