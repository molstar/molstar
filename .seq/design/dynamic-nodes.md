# Dynamic nodes: data that arrives as structure, not only as tables

## Intent

The author, 2026-11-17, on dynamic nodes being deferred: "wouldn't this be needed to have a good streaming protocol
though? i would almost lean to do it now, because it might be really hard to do later."

A source read by region does not only send rows. A transcripts file sends transcripts, each with its own coordinate
system (its residues), a codon map onto the genome, and a source of its protein's domains; a reads file sends coverage
and mismatches as their own annotations; a session file sends tracks. If the protocol, the store and the resolved model
only know tables, adding structure later changes all three. So the protocol carries nodes from the start, and the
first slice exercises them.

## Precedents

- The previous implementation: manifests and patches adding data nodes (`{kind, key, params}`, identity by source and
  key), checked as document nodes, refs that may wait for a node not yet held, a lifetime held by what uses a node,
  item coordinate systems made on demand from an item's `map`, nested sources requested on need with a depth limit,
  track templates (`for_each`) showing dynamic nodes (EC-176 to EC-180; [sources/seqstar](../sources/seqstar/README.md)).
- JBrowse: adapters that return features with subfeatures and their own coordinates; synteny and alignments adapters
  yielding more than rows ([sources/precedents](../sources/precedents/README.md)).

## Decision

- **An answer carries tables and nodes**: the source protocol's answer for a region is its tables plus the data nodes
  the region brings (sequences, units, maps, annotations, nested sources); a source's coverage answer may list nodes
  known up front ([sources.md](sources.md)).
- **Identity is (source, key)**: the same key in a later answer is the same node; a node that leaves and returns keeps
  its identity and its coordinate system's id (EC-176). Address: `/data/<source>/(node :key <key>)`, or a kind's own
  form (`(item :id T1)`, `(sequence :name …)`), never a spec path ([addressing.md](addressing.md)).
- **Checked by the same description**: a dynamic node is validated as a spec node of its kind, by the one language
  description; its refs name spec nodes, held dynamic nodes or nodes of the same answer; anything else is an issue and
  the node is left out. The same key with other params keeps the held node, with an issue. Refs are absolute addresses
  (the language's `ref` takes no relative form), so a runtime writes absolute ones; a relative form would be a core
  repr change, not wanted.
- **Where they live**: in the data store, as entries of their source ([data-store.md](data-store.md)); never in the spec,
  its history or its patches (EC-175).
- **The overlay is a tier of its own** ([composition.md](composition.md)): the static model and the resolved view come
  from the spec (pure, cached by hash) and never hold a dynamic node; the overlay holds the dynamic nodes in the store,
  versioned with it, and the draw lists read both. Systems, maps and
  annotations from the overlay join the model as the static ones do; paths through a dynamic map are composed when the
  overlay's version changes, for the systems it touched, never on hover.
- **Naming before arrival**: a spec may name a dynamic node before it is held; what names it waits (a layer pending, its
  track not shown, noted once per source), never "unplaced".
- **Item coordinate systems on demand**: an item whose rows carry a `map` (a transcript's CDS) gets its own coordinate
  system when something names it or a template selects it, made when a region holding the item arrives (EC-177). The
  system's length is `floor(coding / 3)` (a stop counts as a residue); it stays so until letters arrive.
- **Nested sources on need**: requested only when something needs their data, never past a depth limit (EC-179); wanted
  by region for the shown instances, not only the template's root source.
- **Lifetime**: a node stays while something holds it (a region in the store that carried it, a placement shown, the
  selection, a highlight, a ref, another held node naming it) and is evicted quietly with what it alone named. A node
  is also held while shown (a placed instance, a node whose item system a shown layer uses, a node a held node names):
  `nodes.holdShown` (default true). Why: a shown track must not vanish when its region's tiles are evicted.
  **Shown means in view** (decided, the author, 2026-11-19: "Hold only while in view"; found by CS: a template's
  instance held its node, so every gene's exon track, once shown, stayed, 29 empty tracks after five places): a
  template's instance holds its node only while the node's region meets the view and its margin; out of it the node
  is evicted as any other and its track goes. A track never vanishes while its node's region is on screen, which is
  what `holdShown` is for.
- **Track templates are in version 1**: a template (`for_each` over a source's nodes, by kind, with order and filter) is
  one entry of the resolved view, expanded at layout into a track per matching held node, in key order; an arrival
  re-expands only its template (EC-178). The resolved view stays pure: it holds the template, not its expansion.
  - `for_each.node_kind` (not `kind`, a reserved param key of the description language) is narrowed to `ranges` and
    `values`.
  - An instance's label is the template's label plus its top ancestor's key; nodes have no label param (open: add one
    to the description).
  - Expansion is capped by `maxTemplateTracks` (default 200), with an issue past it.
  - open: today every template re-expands on each overlay version, cheap at this size; the design is that only the
    arrival's own does (EC-178). Recommendation: per-template re-expansion as a followup, with the overlay keyed on
    what it reads of the model.

## Consequences

- The first slice's by-region source is synthetic and asynchronous, and sends a transcript as a dynamic node (its item
  coordinate system, its codon map onto the region, a nested source of domains), shown by a track template; FASTA with
  its index follows ([version-1.md](version-1.md)).
- The store's and the protocol's tests cover arrival, conflict, waiting refs, eviction and return with the same id, and
  a template expanding and contracting, with a fake clock.
